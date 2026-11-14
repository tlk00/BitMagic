"""Interruptible mono depot ambience using standard local audio players."""
import array
import os
import select
import shutil
import signal
import subprocess
import sys
import tempfile
import threading
import time
import wave
from pathlib import Path

AUDIO_RETRY_SECONDS = 30
STREAM_STALL_SECONDS = 10


def stream_command(player, rate):
    """Return the Linux player's command for a mono, signed 16-bit PCM stream."""
    if player == 'pw-play':
        return [player, '--raw', '--rate=' + str(rate), '--channels=1', '--format=s16', '-']
    if player == 'paplay':
        return [player, '--raw', '--rate=' + str(rate), '--channels=1', '--format=s16le']
    if player == 'aplay':
        return [player, '-q', '-t', 'raw', '-f', 'S16_LE', '-r', str(rate), '-c', '1']
    raise ValueError('unsupported streaming player: ' + player)


class StreamPlayer:
    """An agent-owned Linux audio player accepting PCM over a nonblocking pipe."""

    def __init__(self, agent, player, rate, stop):
        """Start player and register its process so Ctrl-C also stops it."""
        self.agent = agent
        self.stop = stop
        self.error_log = tempfile.TemporaryFile(dir=getattr(agent, 'temp_dir', getattr(agent, 'home', None)))
        self.process = None
        try:
            with agent.mutex:
                if stop.is_set() or agent.stopping.is_set():
                    raise RuntimeError('Audio stopped')
                self.process = subprocess.Popen(stream_command(player, rate), stdin=subprocess.PIPE,
                                                stdout=subprocess.DEVNULL, stderr=self.error_log,
                                                start_new_session=True, bufsize=0)
                os.set_blocking(self.process.stdin.fileno(), False)
                agent.children.add(self.process)
        except BaseException:
            if self.process:
                self.process.kill()
                self.process.wait()
            self.error_log.close()
            raise

    def write(self, data):
        """Write all PCM while remaining responsive to shutdown and player failure."""
        remaining = memoryview(data)
        last_progress = time.monotonic()
        while remaining and not self.stop.is_set() and not self.agent.stopping.is_set():
            if self.process.poll() is not None:
                raise RuntimeError(self.error_detail())
            _, writable, _ = select.select([], [self.process.stdin], [], .2)
            if writable:
                try:
                    count = os.write(self.process.stdin.fileno(), remaining)
                except BlockingIOError:
                    continue
                except BrokenPipeError:
                    raise RuntimeError(self.error_detail()) from None
                remaining = remaining[count:]
                last_progress = time.monotonic()
            elif time.monotonic() - last_progress > STREAM_STALL_SECONDS:
                raise RuntimeError('Audio player stopped consuming data')
        if remaining:
            raise RuntimeError('Audio stopped')

    def error_detail(self):
        """Return a short player error without exposing an unbounded log."""
        self.error_log.seek(0)
        detail = self.error_log.read(2048).decode(errors='replace').strip()
        return 'Audio player exited' + (': ' + detail if detail else '')

    def close(self):
        """Drain queued audio, then reap and unregister the player process."""
        proc = self.process
        if proc is None:
            return
        try:
            if proc.stdin:
                proc.stdin.close()
            if (self.stop.is_set() or self.agent.stopping.is_set()) and proc.poll() is None:
                self.agent.signal_group(proc, signal.SIGTERM)
            try:
                proc.wait(timeout=3)
            except subprocess.TimeoutExpired:
                self.agent.signal_group(proc, signal.SIGKILL)
                proc.wait()
        finally:
            with self.agent.mutex:
                self.agent.children.discard(proc)
            self.error_log.close()
            self.process = None


def read_loop(mode):
    """Read a supported mono ambience loop, returning its sample rate and samples."""
    if mode not in ('depot', 'train', 'train-normal'):
        raise ValueError('unknown ambience: ' + mode)
    filenames = {'depot': 'depot.wav', 'train': 'train-ambient.wav',
                 'train-normal': 'train-normal-ambient.wav'}
    asset = Path(__file__).resolve().parent / 'assets' / filenames[mode]
    with wave.open(str(asset), 'rb') as stream:
        rate = stream.getframerate()
        if stream.getnchannels() != 1 or stream.getsampwidth() != 2:
            raise ValueError('depot audio must be 16-bit mono PCM')
        samples = array.array('h', stream.readframes(stream.getnframes()))
    if sys.byteorder != 'little':
        samples.byteswap()
    if not samples:
        raise ValueError('empty ambience audio')
    return rate, samples


def mix_preview_background(path, volume, mode='train-normal'):
    """Overlay a steady moving-train bed on a mono 44.1 kHz preview WAV."""
    rate, background = read_loop(mode)
    with wave.open(str(path), 'rb') as stream:
        if (stream.getnchannels(), stream.getsampwidth(), stream.getframerate()) != (1, 2, 44100):
            raise ValueError('preview audio must be 44.1 kHz mono PCM')
        data = array.array('h', stream.readframes(stream.getnframes()))
    if sys.byteorder != 'little':
        data.byteswap()
    gain = .8 * volume / .35
    for index, value in enumerate(data):
        bed = background[int(index * rate / 44100) % len(background)]
        data[index] = max(-32768, min(32767, round(value + bed * gain)))
    if sys.byteorder != 'little':
        data.byteswap()
    with wave.open(str(path), 'wb') as stream:
        stream.setparams((1, 2, 44100, 0, 'NONE', 'not compressed'))
        stream.writeframes(data.tobytes())


def render_preview(path, volume, mode='depot', seconds=None):
    """Write one scaled ambience loop to path and return its duration in seconds."""
    rate, samples = read_loop(mode)
    if seconds is not None:
        samples = samples[:int(seconds * rate)]
        fade = min(len(samples), rate)
        for i in range(fade):
            samples[len(samples) - fade + i] = int(samples[len(samples) - fade + i] * (fade - i - 1) / fade)
    gain = volume / .35 * 1.15
    scaled = array.array('h', (max(-32768, min(32767, int(value * gain))) for value in samples))
    if sys.byteorder != 'little':
        scaled.byteswap()
    with wave.open(str(path), 'wb') as stream:
        stream.setparams((1, 2, rate, 0, 'NONE', 'not compressed'))
        stream.writeframes(scaled.tobytes())
    return len(scaled) / rate


class Depot:
    """Play a looping recording and fade it around foreground sounds.

    The agent owns all playback children. Linux uses a continuous PCM stream;
    macOS retains the portable one-second WAV playback fallback.
    """

    def __init__(self, agent):
        """Bind the worker to an agent's audio configuration and child registry."""
        self.agent = agent
        self.pause = threading.Event()
        self.duck_level = .25
        self.quiet = threading.Event()
        self.unavailable = threading.Event()
        self.stop = threading.Event()
        self.thread = None
        self.mode = 'depot'

    def set_mode(self, mode):
        """Select an ambience with a quiet transition."""
        if mode not in ('depot', 'train', 'train-normal'):
            raise ValueError('unknown ambience: ' + mode)
        self.mode = mode

    def start(self):
        """Start the optional background worker without delaying the supervisor."""
        self.thread = threading.Thread(target=self.run, daemon=True)
        self.thread.start()

    def mute(self, level=.25):
        """Fade to the requested level and wait until the background is reduced."""
        self.duck_level = level
        self.pause.set()
        while self.thread and self.thread.is_alive() and not self.agent.stopping.is_set():
            if self.unavailable.is_set() or self.quiet.wait(.2):
                break

    def unmute(self):
        """Allow the worker to fade back in after foreground playback."""
        self.quiet.clear()
        self.pause.clear()

    def close(self):
        """Stop the worker and allow its current one-second clip to finish."""
        self.stop.set()
        if self.thread:
            self.thread.join(timeout=4)

    def run(self):
        """Loop the bundled recording, fading out in two seconds and in over three."""
        try:
            current = 'depot'
            rate, samples = read_loop(current)
            if sys.platform.startswith('linux'):
                self.run_stream(current, rate, samples)
                return
            players = [p for p in ('afplay', 'pw-play', 'paplay', 'aplay') if shutil.which(p)]
            gain, position = 0.0, 0
            recovering = False
            with tempfile.TemporaryDirectory(prefix='depot-', dir=getattr(self.agent, 'temp_dir', self.agent.home)) as directory:
                chunk = Path(directory) / 'ambient.wav'
                while not self.stop.is_set():
                    if self.agent.stopping.is_set():
                        self.stop.wait(.1)
                        continue
                    if not players:
                        self.unavailable.set()
                        self.quiet.set()
                        if not recovering:
                            self.agent.say('Depot ambience unavailable; retrying in {} seconds'.format(AUDIO_RETRY_SECONDS))
                            recovering = True
                        self.stop.wait(AUDIO_RETRY_SECONDS)
                        players = [p for p in ('afplay', 'pw-play', 'paplay', 'aplay') if shutil.which(p)]
                        continue
                    if self.mode != current and gain <= 0:
                        current = self.mode
                        rate, samples = read_loop(current)
                        position = 0
                    if self.pause.is_set() and self.mode == current and gain <= self.duck_level + 1e-6:
                        self.quiet.set()
                        self.unavailable.set()
                    elif not self.pause.is_set():
                        self.quiet.clear()
                    target = 0.0 if self.mode != current else (self.duck_level if self.pause.is_set() else 1.0)
                    end_gain = max(target, gain - .5) if gain > target else min(target, gain + 1 / 3)
                    volume = self.agent.args.sound_volume / .35 * 1.15
                    data = array.array('h')
                    edge = max(1, round(rate * .004))
                    for i in range(rate):
                        amplitude = (gain + (end_gain - gain) * i / rate) * volume
                        # External players reopen the device per chunk; taper
                        # both ends to avoid clicks from abrupt sample jumps.
                        amplitude *= min(1.0, i / edge, (rate - 1 - i) / edge)
                        value = int(samples[position] * amplitude)
                        data.append(max(-32768, min(32767, value)))
                        position = (position + 1) % len(samples)
                    if sys.byteorder != 'little':
                        data.byteswap()
                    with wave.open(str(chunk), 'wb') as stream:
                        stream.setparams((1, 2, rate, 0, 'NONE', 'not compressed'))
                        stream.writeframes(data.tobytes())
                    played = False
                    last_error = None
                    for player in players:
                        if self.agent.stopping.is_set() or self.stop.is_set():
                            break
                        try:
                            self.agent.command([player, str(chunk)], timeout=4)
                            played = True
                            break
                        except (OSError, RuntimeError) as error:
                            last_error = error
                            continue
                    if not played:
                        if self.agent.stopping.is_set() or self.stop.is_set():
                            gain = 0.0
                            continue
                        self.agent.say('Depot ambience unavailable; retrying in {} seconds{}'.format(
                            AUDIO_RETRY_SECONDS, ': ' + str(last_error) if last_error else ''))
                        self.quiet.set()
                        recovering = True
                        gain = 0.0
                        self.stop.wait(AUDIO_RETRY_SECONDS)
                        players = [p for p in ('afplay', 'pw-play', 'paplay', 'aplay') if shutil.which(p)]
                        continue
                    if recovering:
                        self.agent.say('Depot ambience restored')
                        recovering = False
                        self.unavailable.clear()
                    gain = end_gain
        except (OSError, ValueError, EOFError, wave.Error) as error:
            self.agent.say('Depot ambience unavailable: {}; continuing without background audio'.format(error))
        finally:
            self.quiet.set()

    def run_stream(self, current, rate, samples):
        """Stream both backgrounds through one Linux player between announcements."""
        players = [p for p in ('pw-play', 'paplay', 'aplay') if shutil.which(p)]
        failed = set()
        stream = None
        recovering = False
        gain, position = 0.0, 0
        next_write = time.monotonic()
        try:
            while not self.stop.is_set():
                if self.agent.stopping.is_set():
                    if stream:
                        stream.close()
                        stream = None
                    gain = 0.0
                    self.stop.wait(.1)
                    continue
                if self.mode != current and gain <= 0:
                    current = self.mode
                    new_rate, samples = read_loop(current)
                    if new_rate != rate and stream:
                        stream.close()
                        stream = None
                    rate, position = new_rate, 0
                if self.pause.is_set() and self.mode == current and gain <= self.duck_level + 1e-6:
                    self.quiet.set()
                elif not self.pause.is_set():
                    self.quiet.clear()
                target = 0.0 if self.mode != current else (self.duck_level if self.pause.is_set() else 1.0)
                count = max(1, rate // 4)
                elapsed = count / rate
                end_gain = (max(target, gain - elapsed / 2) if gain > target
                            else min(target, gain + elapsed / 3))
                volume = self.agent.args.sound_volume / .35 * 1.15
                data = array.array('h')
                for i in range(count):
                    amplitude = (gain + (end_gain - gain) * i / count) * volume
                    value = int(samples[(position + i) % len(samples)] * amplitude)
                    data.append(max(-32768, min(32767, value)))
                if sys.byteorder != 'little':
                    data.byteswap()
                if stream is None:
                    player = next((p for p in players if p not in failed), None)
                    if player is None:
                        self.agent.say('Depot ambience unavailable; retrying in {} seconds'.format(AUDIO_RETRY_SECONDS))
                        self.unavailable.set()
                        self.quiet.set()
                        recovering = True
                        gain = 0.0
                        self.stop.wait(AUDIO_RETRY_SECONDS)
                        players = [p for p in ('pw-play', 'paplay', 'aplay') if shutil.which(p)]
                        failed.clear()
                        continue
                    try:
                        stream = StreamPlayer(self.agent, player, rate, self.stop)
                        next_write = time.monotonic()
                    except (OSError, RuntimeError) as error:
                        failed.add(player)
                        self.agent.say('Background player {} failed: {}'.format(player, error))
                        continue
                try:
                    stream.write(data.tobytes())
                except (OSError, RuntimeError) as error:
                    stream.close()
                    stream = None
                    gain = 0.0
                    if not self.stop.is_set() and not self.agent.stopping.is_set():
                        failed.add(player)
                        self.agent.say('Background player {} failed: {}'.format(player, error))
                    continue
                if recovering:
                    self.agent.say('Depot ambience restored')
                    recovering = False
                    self.unavailable.clear()
                gain = end_gain
                position = (position + count) % len(samples)
                next_write = max(next_write + elapsed, time.monotonic())
                while not self.stop.is_set() and not self.agent.stopping.is_set():
                    delay = next_write - time.monotonic()
                    if delay <= 0:
                        break
                    self.stop.wait(min(.1, delay))
        finally:
            if stream:
                stream.close()

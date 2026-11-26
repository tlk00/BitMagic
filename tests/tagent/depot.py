"""Interruptible mono depot ambience using standard local audio players."""
import array
import math
import os
import random
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
AMBIENCE_FILES = {
    'depot': 'depot.wav',
    'train': 'train-ambient.wav',
    'train-normal': 'train-normal-ambient.wav',
    'train-stress64': 'stress64-carriage.wav',
}


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
    if mode not in AMBIENCE_FILES:
        raise ValueError('unknown ambience: ' + mode)
    return read_mono(Path(__file__).resolve().parent / 'assets' / AMBIENCE_FILES[mode])


def read_mono(path):
    """Load the rate and native-endian samples of a nonempty 16-bit mono WAV."""
    with wave.open(str(path), 'rb') as stream:
        rate = stream.getframerate()
        if stream.getnchannels() != 1 or stream.getsampwidth() != 2:
            raise ValueError('background audio must be 16-bit mono PCM: ' + str(path))
        samples = array.array('h', stream.readframes(stream.getnframes()))
    if sys.byteorder != 'little':
        samples.byteswap()
    if not samples:
        raise ValueError('empty background audio: ' + str(path))
    return rate, samples


class Garratt:
    """Crossfade a one-time departure and occasional Garratt runs into the wagon."""

    def __init__(self, rate, recording):
        """Keep the Garratt PCM at the wagon sample rate for one stress64 stage."""
        self.rate = rate
        self.recording = recording
        self.departed = False
        self.next_at = 0
        self.position = self.length = 0
        self.fade_in = self.fade_out = 1

    def start_if_due(self, now, busy):
        """Start a due recording when foreground audio is clear. Return whether it started."""
        if busy or self.position < self.length or now < self.next_at:
            return False
        self.position = 0
        self.length = len(self.recording) if self.departed else min(len(self.recording), self.rate * 40)
        self.fade_in = self.rate * (12 if self.departed else 3)
        self.fade_out = self.rate * (14 if self.departed else 9)
        self.departed = True
        self.next_at = now + random.uniform(30 * 60, 50 * 60)
        return True

    def mix(self, samples, position, count):
        """Return the next count wagon samples with a constant-power Garratt crossfade."""
        for index in range(count):
            wagon = samples[(position + index) % len(samples)]
            offset = self.position + index
            if offset < self.length:
                weight = min(1, offset / self.fade_in, (self.length - 1 - offset) / self.fade_out)
                wagon = wagon * math.sqrt(1 - weight) + self.recording[offset] * math.sqrt(weight)
            yield wagon


def read_garratt(rate):
    """Read the bundled Garratt recording at the wagon rate."""
    path = Path(__file__).resolve().parent / 'assets' / 'stress64-garratt.wav'
    recorded_rate, recording = read_mono(path)
    if recorded_rate != rate:
        raise ValueError('Garratt and wagon audio rates differ')
    return recording


def ambient_chunk(samples, position, count, gain, end_gain, volume, garratt=None, edge=0):
    """Render count samples with a gain ramp and optional train mix and edge tapers."""
    if garratt and garratt.position < garratt.length:
        values = garratt.mix(samples, position, count)
    else:
        values = (samples[(position + index) % len(samples)] for index in range(count))
    data = array.array('h')
    for index, value in enumerate(values):
        amplitude = (gain + (end_gain - gain) * index / count) * volume
        if edge:
            # Standalone players reopen the audio device per chunk.
            amplitude *= min(1, index / edge, (count - 1 - index) / edge)
        data.append(max(-32768, min(32767, int(value * amplitude))))
    if sys.byteorder != 'little':
        data.byteswap()
    return data


def render_stress64_preview(path, volume, voice_path=None):
    """Write a wagon, departure and shortened return demo with optional conductor audio."""
    rate, wagon = read_loop('train-stress64')
    garratt = Garratt(rate, read_garratt(rate)[:52 * rate])
    voice = array.array('h')
    if voice_path is not None:
        voice_rate, voice = read_mono(voice_path)
        if voice_rate != 2 * rate:
            raise ValueError('stress64 demo voice must be at twice the wagon sample rate')
        voice = array.array('h', ((left + right) // 2 for left, right in zip(voice[::2], voice[1::2])))
    voice_start = 45 * rate
    voice_end = voice_start + len(voice)
    return_at = max(57, math.ceil(voice_end / rate) + 4)
    seconds = return_at + 58
    end = seconds * rate - 1
    gain = volume / .35 * 1.15
    garratt.start_if_due(0, False)
    with wave.open(str(path), 'wb') as output:
        output.setparams((1, 2, rate, 0, 'NONE', 'not compressed'))
        for second in range(seconds):
            if second == return_at:
                garratt.start_if_due(garratt.next_at, False)
            position = second * rate
            data = array.array('h')
            for index, sample in enumerate(garratt.mix(wagon, position, rate), position):
                duck = 0
                if voice and voice_start - rate < index < voice_end + 3 * rate:
                    duck = min(1, (index - voice_start + rate) / rate,
                               (voice_end + 3 * rate - index) / (3 * rate))
                spoken = voice[index - voice_start] if voice_start <= index < voice_end else 0
                fade = min(1, (end - index) / (4 * rate))
                value = int((sample * gain * (1 - .75 * duck) + spoken) * fade)
                data.append(max(-32768, min(32767, value)))
            garratt.position = min(garratt.length, garratt.position + rate)
            if sys.byteorder != 'little':
                data.byteswap()
            output.writeframesraw(data.tobytes())
    return seconds


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
        self.request = ('depot', 0)

    def set_mode(self, mode):
        """Select an ambience with a quiet transition."""
        if mode not in AMBIENCE_FILES:
            raise ValueError('unknown ambience: ' + mode)
        current, generation = self.request
        if mode != current:
            self.request = (mode, generation + 1)

    def load_loop(self, request):
        """Load the requested background and optional Garratt recording for this stage."""
        mode, _ = request
        rate, samples = read_loop(mode)
        garratt = Garratt(rate, read_garratt(rate)) if mode == 'train-stress64' else None
        return rate, samples, garratt

    def start_garratt(self, garratt, current):
        """Defer a due Garratt departure or passage until the foreground queue clears."""
        if garratt is None or current != self.request:
            return
        queued = self.agent.audio_thread
        busy = self.pause.is_set() or (queued is not None and queued.is_alive())
        departure = not garratt.departed
        if garratt.start_if_due(time.monotonic(), busy):
            self.agent.say('STRESS64 AMBIENCE | Garratt ' + ('departure' if departure else 'return'))

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
            current = self.request
            rate, samples, garratt = self.load_loop(current)
            if sys.platform.startswith('linux'):
                self.run_stream(current, rate, samples, garratt)
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
                    if self.request != current and gain <= 0:
                        current = self.request
                        rate, samples, garratt = self.load_loop(current)
                        position = 0
                    self.start_garratt(garratt, current)
                    if self.pause.is_set() and self.request == current and gain <= self.duck_level + 1e-6:
                        self.quiet.set()
                    elif not self.pause.is_set():
                        self.quiet.clear()
                    target = 0.0 if self.request != current else (self.duck_level if self.pause.is_set() else 1.0)
                    end_gain = max(target, gain - .5) if gain > target else min(target, gain + 1 / 3)
                    volume = self.agent.args.sound_volume / .35 * 1.15
                    edge = max(1, round(rate * .004))
                    data = ambient_chunk(samples, position, rate, gain, end_gain, volume, garratt, edge)
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
                        self.unavailable.set()
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
                    position = (position + rate) % len(samples)
                    if garratt:
                        garratt.position = min(garratt.length, garratt.position + rate)
        except (OSError, ValueError, EOFError, wave.Error) as error:
            self.agent.say('Depot ambience unavailable: {}; continuing without background audio'.format(error))
        finally:
            self.quiet.set()

    def run_stream(self, current, rate, samples, garratt):
        """Stream the current loop and optional train mix through one Linux player."""
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
                if self.request != current and gain <= 0:
                    current = self.request
                    new_rate, samples, garratt = self.load_loop(current)
                    if new_rate != rate and stream:
                        stream.close()
                        stream = None
                    rate, position = new_rate, 0
                self.start_garratt(garratt, current)
                if self.pause.is_set() and self.request == current and gain <= self.duck_level + 1e-6:
                    self.quiet.set()
                elif not self.pause.is_set():
                    self.quiet.clear()
                target = 0.0 if self.request != current else (self.duck_level if self.pause.is_set() else 1.0)
                count = max(1, rate // 4)
                elapsed = count / rate
                end_gain = (max(target, gain - elapsed / 2) if gain > target
                            else min(target, gain + elapsed / 3))
                volume = self.agent.args.sound_volume / .35 * 1.15
                data = ambient_chunk(samples, position, count, gain, end_gain, volume, garratt)
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
                if garratt:
                    garratt.position = min(garratt.length, garratt.position + count)
                next_write = max(next_write + elapsed, time.monotonic())
                while not self.stop.is_set() and not self.agent.stopping.is_set():
                    delay = next_write - time.monotonic()
                    if delay <= 0:
                        break
                    self.stop.wait(min(.1, delay))
        finally:
            if stream:
                stream.close()

"""Interruptible mono depot ambience using standard local audio players."""
import array
import shutil
import sys
import tempfile
import threading
import wave
from pathlib import Path


def read_loop(mode):
    """Read a supported mono ambience loop, returning its sample rate and samples."""
    if mode not in ('depot', 'train'):
        raise ValueError('unknown ambience: ' + mode)
    asset = Path(__file__).resolve().parent / 'assets' / (
        'depot.wav' if mode == 'depot' else 'train-ambient.wav')
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


def render_preview(path, volume, mode='depot'):
    """Write one scaled ambience loop to path and return its duration in seconds."""
    rate, samples = read_loop(mode)
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

    The agent owns all playback children; the worker uses one-second WAV chunks
    so fade requests are portable across the supported external players.
    """

    def __init__(self, agent):
        """Bind the worker to an agent's audio configuration and child registry."""
        self.agent = agent
        self.pause = threading.Event()
        self.quiet = threading.Event()
        self.stop = threading.Event()
        self.thread = None
        self.mode = 'depot'

    def set_mode(self, mode):
        """Select depot or moving-train ambience with a quiet transition."""
        if mode not in ('depot', 'train'):
            raise ValueError('unknown ambience: ' + mode)
        self.mode = mode

    def start(self):
        """Start the optional background worker without delaying the supervisor."""
        self.thread = threading.Thread(target=self.run, daemon=True)
        self.thread.start()

    def mute(self):
        """Fade out and wait briefly until the background player is silent."""
        self.pause.set()
        while self.thread and self.thread.is_alive() and not self.agent.stopping.is_set():
            if self.quiet.wait(.2):
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
            players = [p for p in ('afplay', 'pw-play', 'paplay', 'aplay') if shutil.which(p)]
            if not players:
                return
            gain, position = 0.0, 0
            with tempfile.TemporaryDirectory(prefix='depot-', dir=self.agent.home) as directory:
                chunk = Path(directory) / 'ambient.wav'
                while not self.stop.is_set():
                    if self.agent.stopping.is_set():
                        self.stop.wait(.1)
                        continue
                    if self.mode != current and gain <= 0:
                        current = self.mode
                        rate, samples = read_loop(current)
                        position = 0
                    if self.pause.is_set() and gain <= 0:
                        self.quiet.set()
                        self.stop.wait(.1)
                        continue
                    self.quiet.clear()
                    target = 0.0 if self.pause.is_set() or self.mode != current else 1.0
                    end_gain = max(target, gain - .5) if target == 0 else min(target, gain + 1 / 3)
                    volume = self.agent.args.sound_volume / .35 * 1.15
                    data = array.array('h')
                    for i in range(rate):
                        amplitude = (gain + (end_gain - gain) * i / rate) * volume
                        value = int(samples[position] * amplitude)
                        data.append(max(-32768, min(32767, value)))
                        position = (position + 1) % len(samples)
                    if sys.byteorder != 'little':
                        data.byteswap()
                    with wave.open(str(chunk), 'wb') as stream:
                        stream.setparams((1, 2, rate, 0, 'NONE', 'not compressed'))
                        stream.writeframes(data.tobytes())
                    played = False
                    for player in players:
                        if self.agent.stopping.is_set() or self.stop.is_set():
                            break
                        try:
                            self.agent.command([player, str(chunk)], timeout=4)
                            played = True
                            break
                        except (OSError, RuntimeError):
                            continue
                    if not played:
                        if self.agent.stopping.is_set() or self.stop.is_set():
                            gain = 0.0
                            continue
                        self.agent.say('Depot ambience unavailable; continuing without background audio')
                        return
                    gain = end_gain
        except (OSError, ValueError, EOFError, wave.Error) as error:
            self.agent.say('Depot ambience unavailable: {}; continuing without background audio'.format(error))
        finally:
            self.quiet.set()

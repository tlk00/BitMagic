"""PCM synthesis and nonblocking local playback checks; no audio device needed."""
import array
import tempfile
import unittest
import wave
from pathlib import Path
from unittest.mock import patch

import soundtrack
import tagent


class SoundtrackTests(unittest.TestCase):
    """Validate both themed sound banks and their supervisor integration."""

    def test_banks(self):
        """Every event is valid mono PCM, distinct per theme, with bounded peaks."""
        with tempfile.TemporaryDirectory() as directory:
            for event in soundtrack.EVENTS:
                banks = []
                for theme in soundtrack.THEMES:
                    path = Path(directory) / (theme + '.wav')
                    duration = soundtrack.render(path, theme, event, 0.35, 3, 8)
                    if theme == 'steampunk' and event == 'depart':
                        self.assertEqual(duration, 23.0)
                    if theme == 'steampunk' and event == 'building':
                        self.assertGreater(duration, soundtrack.morse_intervals(soundtrack.TELEGRAM)[-1][1])
                    with wave.open(str(path), 'rb') as stream:
                        self.assertEqual(stream.getparams()[:3], (1, 2, 44100))
                        self.assertAlmostEqual(stream.getnframes() / 44100, duration, places=4)
                        data = stream.readframes(stream.getnframes())
                    values = array.array('h', data)
                    self.assertLess(max(abs(value) for value in values), 16000)
                    self.assertGreater(max(values), 0)
                    banks.append(data)
                # SOS is intentionally the same timing with different horn pitch.
                self.assertNotEqual(*banks)

    def test_telegram_timing(self):
        """Keying preserves dash/dot and word gaps, and encodes the full message."""
        self.assertEqual(soundtrack.morse_intervals('N E', unit=1, lead=0),
                         [(0, 3), (4, 5), (12, 13)])
        self.assertEqual(soundtrack.TELEGRAM,
                         'New BitMagic updates, new BitMagic updates')
        keyed = soundtrack.morse_intervals(soundtrack.TELEGRAM)
        self.assertAlmostEqual(keyed[-1][1] + .8, 38.375)
        self.assertAlmostEqual(soundtrack.morse_intervals('E')[0][1] - .3, .105)
        with self.assertRaises(ValueError):
            soundtrack.morse_intervals('!')

    def test_key_track_has_only_contact_impacts(self):
        """A held dash is silent between its two dry impacts, independent of steam."""
        track = soundtrack.telegraph_track([(0.1, 0.4)], 0.6)
        rate = soundtrack.RATE
        self.assertTrue(any(track[int(.1 * rate):int(.13 * rate)]))
        self.assertFalse(any(track[int(.15 * rate):int(.39 * rate)]))
        self.assertTrue(any(track[int(.4 * rate):int(.43 * rate)]))
        self.assertEqual(track, soundtrack.telegraph_track([(0.1, 0.4)], 0.6))

    def test_steam_crescendo_without_key(self):
        """The independent steam stem grows smoothly across the message."""
        import math
        duration = soundtrack.morse_intervals(soundtrack.TELEGRAM)[-1][1] + .8
        silent_key = array.array('f', [0]) * int(duration * soundtrack.RATE)
        with tempfile.TemporaryDirectory() as directory, \
             patch.object(soundtrack, 'telegraph_track', return_value=silent_key):
            path = Path(directory) / 'steam.wav'
            soundtrack.render(path, 'steampunk', 'building', 1)
            with wave.open(str(path), 'rb') as stream:
                values = array.array('h', stream.readframes(stream.getnframes()))
            def rms(start, end):
                """Return RMS amplitude in a mono sample window, in seconds."""
                window = values[int(start * soundtrack.RATE):int(end * soundtrack.RATE)]
                return math.sqrt(sum(value * value for value in window) / len(window))
            self.assertGreater(rms(duration * .8, duration * .9),
                               4 * rms(duration * .05, duration * .15))

    def test_occupancy_and_silence(self):
        """Pulse rate rises with worker occupancy; zero volume produces silence."""
        self.assertLess(soundtrack.pulse_rate(1, 8), soundtrack.pulse_rate(8, 8))
        self.assertEqual(soundtrack.pulse_rate(100, 8), soundtrack.pulse_rate(8, 8))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'silent.wav'
            soundtrack.render(path, 'steampunk', 'work', 0)
            with wave.open(str(path), 'rb') as stream:
                self.assertFalse(any(stream.readframes(stream.getnframes())))

    def test_demo_sequence(self):
        """Demo bypasses checkout validation and exports the full theme sequence."""
        with tempfile.TemporaryDirectory() as directory, patch('sys.argv', [
                'tagent', '--sound-demo', '--sound-theme', 'dieselpunk',
                '--jobs', '8', '--state-dir', directory]):
            agent = tagent.Agent(tagent.arguments())
            with patch.object(agent, 'validate') as validate, \
                 patch.object(tagent.shutil, 'which', return_value='/fake/afplay'), \
                 patch.object(agent, 'command', return_value='') as command:
                self.assertEqual(agent.run(), 0)
                validate.assert_not_called()
                self.assertEqual(command.call_count, 8)
                self.assertTrue((Path(directory) / 'sound-dieselpunk-work-1.wav').exists())
                self.assertTrue((Path(directory) / 'sound-dieselpunk-work-8.wav').exists())
                with wave.open(str(Path(directory) / 'demo-dieselpunk.wav')) as stream:
                    self.assertEqual(stream.getnchannels(), 1)
                    self.assertGreater(stream.getnframes() / soundtrack.RATE, 20)

    def test_steam_demo_includes_mono_dispatch(self):
        """Demo plays each variant explicitly and combines mono clips correctly."""
        with tempfile.TemporaryDirectory() as directory, patch('sys.argv', [
                'tagent', '--sound-demo', '--sound-theme', 'steampunk',
                '--jobs', '8', '--state-dir', directory]):
            agent = tagent.Agent(tagent.arguments())
            def short_clip(path, theme, event, volume, workers, jobs):
                """Stand in for long synthesized events using a short mono clip."""
                with wave.open(str(path), 'wb') as stream:
                    stream.setparams((1, 2, soundtrack.RATE, 0, 'NONE', 'not compressed'))
                    stream.writeframes(b'\0' * (soundtrack.RATE // 10 * 2))
                return .1
            with patch.object(soundtrack, 'render', side_effect=short_clip), \
                 patch.object(tagent.shutil, 'which', return_value='/fake/afplay'), \
                 patch.object(agent, 'command', return_value='') as command, \
                 patch.object(tagent.random, 'choice', side_effect=AssertionError('random demo')):
                self.assertEqual(agent.run(), 0)
                self.assertEqual(command.call_count, 10)
                with wave.open(str(Path(directory) / 'demo-steampunk.wav')) as stream:
                    self.assertEqual(stream.getnchannels(), 1)
                    self.assertGreater(stream.getnframes() / stream.getframerate(), 5)
                with wave.open(str(Path(directory) / 'sound-steampunk-dispatch.wav')) as stream:
                    self.assertEqual(stream.getnchannels(), 1)

    def test_steam_train_schedule(self):
        """Departure starts the five-minute timer; other signals cannot reset it."""
        with tempfile.TemporaryDirectory() as directory, patch('sys.argv', [
                'tagent', '--sound-theme', 'steampunk', '--state-dir', directory]):
            agent = tagent.Agent(tagent.arguments())
            clock = [100.0]
            played = []
            with patch.object(tagent.time, 'monotonic', side_effect=lambda: clock[0]), \
                 patch.object(tagent.random, 'choice', side_effect=['work', 'dispatch']) as coin, \
                 patch.object(agent, 'play_sound', side_effect=lambda event, workers: played.append(event)):
                agent.sound('work')
                self.assertIsNone(agent.audio_thread)  # No train outside a test cycle.
                agent.start_test_sound()
                agent.audio_thread.join()
                self.assertEqual(played, ['depart'])
                clock[0] = 399
                agent.sound('work')
                self.assertEqual(played, ['depart'])
                agent.sound('build')
                agent.audio_thread.join()
                clock[0] = 400
                agent.sound('work')
                agent.audio_thread.join()
                self.assertEqual(played, ['depart', 'build', 'work'])
                clock[0] = 699
                agent.sound('work')
                self.assertEqual(played.count('work'), 1)
                clock[0] = 700
                agent.sound('work')
                agent.audio_thread.join()
                self.assertEqual(played, ['depart', 'build', 'work', 'dispatch'])
                self.assertEqual(coin.call_count, 2)
                agent.args.sound_theme = 'dieselpunk'
                clock[0] = 800
                agent.sound('work')
                agent.audio_thread.join()
                self.assertEqual(played[-1], 'work')
                self.assertEqual(coin.call_count, 2)

    def test_dispatch_asset_mono_and_volume(self):
        """The approved portable clip remains mono and obeys the volume setting."""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'dispatch.wav'
            source_path = Path(soundtrack.__file__).parent / 'assets' / 'dispatch.wav'
            with wave.open(str(source_path), 'rb') as source:
                reference = source.readframes(source.getnframes())
            duration = soundtrack.render_dispatch(path, .5)
            with wave.open(str(path), 'rb') as output:
                self.assertEqual(output.getparams()[:3], (1, 2, 44100))
                self.assertEqual(output.readframes(output.getnframes()), reference)
                self.assertGreater(duration, 4)
            soundtrack.render_dispatch(path, 0)
            with wave.open(str(path), 'rb') as output:
                self.assertFalse(any(output.readframes(output.getnframes())))

    def test_dispatch_preview_requires_steampunk(self):
        """The dispatcher voice is unavailable in diesel-punk mode."""
        with patch('sys.argv', ['tagent', '--sound-preview', 'dispatch',
                                '--sound-theme', 'dieselpunk']), \
             patch('sys.stderr'):
            with self.assertRaises(SystemExit) as error:
                tagent.arguments()
            self.assertEqual(error.exception.code, 2)

    def test_preview_without_checkout_or_audio_device(self):
        """Preview bypasses Git validation and dispatches local player asynchronously."""
        with tempfile.TemporaryDirectory() as directory, patch('sys.argv', [
                'tagent', '--sound-preview', 'work', '--sound-theme', 'steampunk',
                '--state-dir', directory]):
            agent = tagent.Agent(tagent.arguments())
            with patch.object(agent, 'validate') as validate, \
                 patch.object(tagent.shutil, 'which', return_value='/fake/afplay'), \
                 patch.object(agent, 'command', return_value='') as command:
                self.assertEqual(agent.run(), 0)
                agent.audio_thread.join(timeout=10)
                self.assertFalse(agent.audio_thread.is_alive())
                validate.assert_not_called()
                self.assertEqual(command.call_args.args[0][0], 'afplay')
                # Playback must allow the entire slow train to pass, not cut it
                # off at the former eight-second operational timeout.
                with wave.open(command.call_args.args[0][1], 'rb') as stream:
                    duration = stream.getnframes() / stream.getframerate()
                self.assertGreater(command.call_args.kwargs['timeout'], duration)


if __name__ == '__main__':
    unittest.main()

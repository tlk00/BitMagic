"""PCM synthesis and nonblocking local playback checks; no audio device needed."""
import array
import tempfile
import threading
import sys
import unittest
import wave
from pathlib import Path
from unittest.mock import patch

import soundtrack
import tagent


class SoundtrackTests(unittest.TestCase):
    """Validate the steampunk sound bank and their supervisor integration."""

    def setUp(self):
        """Keep host Piper installations from changing unrelated fallback tests."""
        self.piper_environment = patch.dict('os.environ', {
            'PIPER_PYTHON': '/nonexistent/tagent-test-piper-python',
            'PIPER_MODEL': '/nonexistent/tagent-test-piper.onnx'})
        self.piper_environment.start()
        self.addCleanup(self.piper_environment.stop)

    def test_station_announcement_waits_and_ducks_bed(self):
        """Station speech starts at eight seconds and reduces the bed by 20 percent."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bed, voice, result = (root / name for name in ('bed.wav', 'voice.wav', 'out.wav'))
            for path, duration, value in ((bed, 12, 1000), (voice, 2, 2000)):
                with wave.open(str(path), 'wb') as stream:
                    stream.setparams((1, 2, 44100, 0, 'NONE', 'not compressed'))
                    stream.writeframes(array.array('h', [value] * (duration * 44100)).tobytes())
            soundtrack.render_station_visit(result, voice, .35 / 1.15, bed)
            with wave.open(str(result), 'rb') as stream:
                data = array.array('h', stream.readframes(stream.getnframes()))
            self.assertEqual(data[6 * 44100], 1000)
            self.assertEqual(data[8 * 44100], 2800)
            self.assertEqual(data[10 * 44100], 800)
            self.assertGreater(data[11 * 44100], 800)

    def test_banks(self):
        """Every event is valid mono PCM with bounded peaks."""
        with tempfile.TemporaryDirectory() as directory:
            for event in soundtrack.EVENTS:
                banks = []
                for theme in soundtrack.THEMES:
                    path = Path(directory) / (theme + '.wav')
                    duration = soundtrack.render(path, theme, event, 0.35, 3, 8)
                    if theme == 'steampunk' and event == 'depart':
                        self.assertAlmostEqual(duration, 56.1, places=4)
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

    def test_station_music_crossfades_into_departure(self):
        """The first 30% of the train overlaps the end of the longer march."""
        with tempfile.TemporaryDirectory() as directory:
            plain = Path(directory) / 'departure.wav'
            combined = Path(directory) / 'station-and-departure.wav'
            soundtrack.render(plain, 'steampunk', 'depart', .35, 2, 8, station_music=False)
            duration = soundtrack.render(combined, 'steampunk', 'depart', .35, 2, 8)
            self.assertAlmostEqual(duration, 56.1, places=4)
            with wave.open(str(plain), 'rb') as source:
                source.setpos(round(6.9 * soundtrack.RATE))
                tail = array.array('h', source.readframes(source.getnframes()))
                tail = array.array('h', (round(value * .9) for value in tail)).tobytes()
            with wave.open(str(combined), 'rb') as source:
                self.assertEqual(source.getnchannels(), 1)
                self.assertTrue(any(source.readframes(33 * soundtrack.RATE)))
                source.setpos(40 * soundtrack.RATE)
                self.assertEqual(source.readframes(source.getnframes()), tail)

    def test_train_stop_preview_is_mono_and_faded(self):
        """The locomotive stop has quiet boundaries and honors output volume."""
        import tempfile
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'stop.wav'
            duration = soundtrack.render_train_stop(path, .35)
            self.assertGreater(duration, 5)
            self.assertLess(duration, 7)
            with wave.open(str(path)) as stream:
                self.assertEqual((stream.getnchannels(), stream.getframerate()), (1, 44100))
                samples = array.array('h', stream.readframes(stream.getnframes()))
            self.assertLess(abs(samples[0]), 3)
            self.assertLess(abs(samples[-1]), 3)
            self.assertGreater(max(abs(sample) for sample in samples), 500)
            self.assertLessEqual(max(abs(sample) for sample in samples), 12000)

    def test_station_and_train_fade_out(self):
        """March and both train effects taper smoothly towards silence."""
        import math
        def rms(values):
            """Measure RMS in a nonempty mono sample window."""
            return math.sqrt(sum(value * value for value in values) / len(values))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'sound.wav'
            soundtrack.render(path, 'steampunk', 'depart', .35, station_music=True)
            with wave.open(str(path)) as stream:
                samples = array.array('h', stream.readframes(stream.getnframes()))
            def window(start, end):
                return samples[int(start * soundtrack.RATE):int(end * soundtrack.RATE)]
            self.assertLess(rms(window(55.8, 56.0)), rms(window(50, 51)) * .2)
            soundtrack.render(path, 'steampunk', 'work', .35)
            with wave.open(str(path)) as stream:
                samples = array.array('h', stream.readframes(stream.getnframes()))
            self.assertLess(rms(window(10.2, 10.4)), rms(window(8, 9)) * .1)

    def test_phase_announcement_precedes_station_and_departure(self):
        """Both recorded announcements precede the same station/train sequence."""
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory) / 'base.wav'
            soundtrack.render(base, 'steampunk', 'depart', .35, 2, 8)
            with wave.open(str(base), 'rb') as source:
                expected = source.readframes(source.getnframes())
            intros = []
            for phase in ('perf', 'stress'):
                path = Path(directory) / (phase + '.wav')
                duration = soundtrack.render(path, 'steampunk', 'depart', .35, 2, 8,
                                             announcement=phase)
                with wave.open(str(Path(soundtrack.__file__).parent / 'assets' /
                                   ('station-' + phase + '.wav')), 'rb') as asset:
                    intro_frames = round(asset.getnframes() * soundtrack.RATE / asset.getframerate())
                with wave.open(str(path), 'rb') as source:
                    self.assertEqual(source.getnchannels(), 1)
                    intros.append(source.readframes(intro_frames))
                    self.assertFalse(any(source.readframes(round(.3 * soundtrack.RATE))))
                    self.assertEqual(source.readframes(source.getnframes()), expected)
                self.assertGreater(duration, 53)
            self.assertNotEqual(*intros)

    def test_phase_audio_queue_does_not_block_and_prioritizes_failure(self):
        """Queued phase sound is nonblocking; failures discard pending chatter."""
        import threading
        with tempfile.TemporaryDirectory() as directory, patch('sys.argv', [
                'tagent', '--sound-theme', 'steampunk', '--state-dir', directory]):
            agent = tagent.Agent(tagent.arguments())
            started = threading.Event()
            release = threading.Event()
            played = []
            def play(event, workers):
                """Hold the current clip so the following phase must queue."""
                played.append(event)
                if event == 'building':
                    started.set()
                    release.wait(5)
            try:
                with patch.object(agent, 'play_sound', side_effect=play):
                    agent.sound('building')
                    self.assertTrue(started.wait(2))
                    agent.sound('depart-perf')
                    agent.sound('fail')
                    release.set()
                    agent.audio_thread.join(timeout=5)
                    self.assertFalse(agent.audio_thread.is_alive())
                    self.assertEqual(played, ['building', 'fail'])
            finally:
                release.set()

    def test_stale_phase_chatter_is_dropped_and_repeats_coalesce(self):
        """Pending old-phase reports never play, and repeated reports share one slot."""
        with tempfile.TemporaryDirectory() as directory, patch('sys.argv', [
                'tagent', '--state-dir', directory]):
            agent = tagent.Agent(tagent.arguments())
            started, release = threading.Event(), threading.Event()
            played = []
            def play(event, workers):
                played.append(event)
                if event == 'building':
                    started.set()
                    release.wait(3)
            try:
                with patch.object(agent, 'play_sound', side_effect=play):
                    agent.sound('building')
                    self.assertTrue(started.wait(2))
                    agent.sound('station-perf')
                    pending = agent.audio_thread
                    for _ in range(100):
                        agent.sound('station-perf')
                    self.assertIs(agent.audio_thread, pending)
                    with agent.audio_mutex:
                        agent.audio_phase_generation += 1
                    agent.sound('station-stress')
                    release.set()
                    self.assertTrue(agent.wait_for_audio())
                    self.assertEqual(played, ['building', 'station-stress'])
            finally:
                release.set()

    def test_build_status_expires_when_build_finishes(self):
        """A queued old build report cannot follow the final build summary."""
        with tempfile.TemporaryDirectory() as directory, patch('sys.argv', [
                'tagent', '--state-dir', directory]):
            agent = tagent.Agent(tagent.arguments())
            started, release = threading.Event(), threading.Event()
            played = []
            def play(event, workers):
                played.append(event)
                if event == 'building':
                    started.set()
                    release.wait(3)
            try:
                with patch.object(agent, 'play_sound', side_effect=play):
                    agent.sound('building')
                    self.assertTrue(started.wait(2))
                    agent.sound('build-warning', workers=1)
                    agent.finish_build(True)
                    release.set()
                    self.assertTrue(agent.wait_for_audio())
                    self.assertEqual(played, ['building', 'build-final-ok'])
            finally:
                release.set()

    def test_repeated_failure_preserves_distress_exchange(self):
        """A second failure does not cancel an SOS already queued or playing."""
        with tempfile.TemporaryDirectory() as directory, patch('sys.argv', [
                'tagent', '--state-dir', directory]):
            agent = tagent.Agent(tagent.arguments())
            started, release = threading.Event(), threading.Event()
            played = []
            def play(event, workers):
                played.append(event)
                if event == 'fail':
                    started.set()
                    release.wait(3)
                return True
            try:
                with patch.object(agent, 'play_sound', side_effect=play):
                    agent.sound('fail')
                    self.assertTrue(started.wait(2))
                    distress = agent.audio_thread
                    agent.sound('fail')
                    self.assertIs(agent.audio_thread, distress)
                    release.set()
                    self.assertTrue(agent.wait_for_audio())
                    self.assertEqual(played, ['fail', 'distress-engineer', 'distress-dispatch'])
            finally:
                release.set()

    def test_finish_audio_drains_active_exchange_and_rejects_new_chatter(self):
        """Shutdown waits for foreground work while closing admission atomically."""
        with tempfile.TemporaryDirectory() as directory, patch('sys.argv', [
                'tagent', '--state-dir', directory]):
            agent = tagent.Agent(tagent.arguments())
            started, release, finished = threading.Event(), threading.Event(), threading.Event()
            played = []
            def play(event, workers):
                played.append(event)
                started.set()
                release.wait(3)
            def finish():
                agent.finish_audio()
                finished.set()
            try:
                with patch.object(agent, 'play_sound', side_effect=play):
                    agent.sound('build')
                    self.assertTrue(started.wait(2))
                    thread = threading.Thread(target=finish)
                    thread.start()
                    for _ in range(100):
                        with agent.audio_mutex:
                            if agent.audio_closing:
                                break
                        threading.Event().wait(.01)
                    self.assertTrue(agent.audio_closing)
                    agent.sound('repo-unchanged')
                    self.assertFalse(finished.is_set())
                    release.set()
                    thread.join(2)
                    self.assertTrue(finished.is_set())
                    self.assertEqual(played, ['build'])
            finally:
                release.set()

    def test_finish_audio_cancels_delayed_joker_reply(self):
        """Shutdown closes the queue before a delayed reply can enqueue audio."""
        with tempfile.TemporaryDirectory() as directory, patch('sys.argv', [
                'tagent', '--state-dir', directory]):
            agent = tagent.Agent(tagent.arguments())
            self.assertTrue(agent.finish_audio())
            with patch.object(tagent.random, 'uniform', return_value=0), \
                 patch.object(agent, 'play_after') as play:
                agent.queue_joker_reply(agent.audio_epoch, agent.audio_phase_generation)
                agent.sound('build')
                play.assert_not_called()
                self.assertIsNone(agent.audio_thread)

    def test_foreground_sound_mutes_depot(self):
        """Foreground playback starts after the depot fades and then restores it."""
        with tempfile.TemporaryDirectory() as directory, patch('sys.argv', [
                'tagent', '--sound-theme', 'steampunk', '--state-dir', directory]):
            agent = tagent.Agent(tagent.arguments())
            order = []
            class Background:
                """Track the ordering of foreground and background playback."""
                def mute(self):
                    order.append('muted')
                def unmute(self):
                    order.append('resumed')
            agent.depot = Background()
            with patch.object(agent, 'play_sound', side_effect=lambda *_: order.append('sound')):
                agent.sound('build-clear')
                agent.audio_thread.join(2)
            self.assertEqual(order, ['muted', 'sound', 'resumed'])

    def test_station_reports_repeat_at_random_intervals(self):
        """Each phase uses a bounded timer and waits for foreground audio."""
        with tempfile.TemporaryDirectory() as directory, patch('sys.argv', [
                'tagent', '--state-dir', directory]):
            agent = tagent.Agent(tagent.arguments())
            agent.audio_available = True
            with patch.object(tagent.time, 'monotonic', return_value=100), \
                 patch.object(tagent.random, 'uniform', return_value=750) as interval, \
                 patch.object(agent, 'sound') as sound, \
                 patch.object(agent, 'start_conductor_download'):
                agent.start_test_sound('perf')
                self.assertEqual(agent.next_station_report, 850)
                agent.report_station_status()
                self.assertEqual(sound.call_count, 2)
                with patch.object(tagent.time, 'monotonic', return_value=850):
                    with patch.object(agent, 'audio_thread') as active_audio:
                        active_audio.is_alive.return_value = True
                        agent.report_station_status()
                        self.assertEqual(sound.call_count, 2)
                        self.assertEqual(agent.next_station_report, 850)
                    agent.report_station_status()
                    sound.assert_called_with('station-perf')
                    self.assertEqual(agent.next_station_report, 1600)
                    agent.report_station_status()
                    self.assertEqual(sound.call_count, 3)
                agent.start_test_sound('stress')
                self.assertEqual(agent.station_phase, 'stress')
                with patch.object(tagent.time, 'monotonic', return_value=850):
                    agent.report_station_status()
                    sound.assert_called_with('station-stress')
                self.assertTrue(all(call.args == (600, 1200) for call in interval.call_args_list))

    def test_station_announcements_use_conductor_and_express_for_both_fast_phases(self):
        """Destinations and fast-mode workload messages match each testing phase."""
        with tempfile.TemporaryDirectory() as directory, patch('sys.argv', [
                'tagent', '--state-dir', directory, '--fast']):
            agent = tagent.Agent(tagent.arguments())
            with patch.object(agent, 'synthesize_announcement', return_value=1) as speech:
                agent.render_station_report(Path(directory) / 'out.wav', False)
                self.assertIn('Performance Test Station', speech.call_args.args[1])
                self.assertIn('express service', speech.call_args.args[1])
                self.assertIn('reduced repetitions and fixtures', speech.call_args.args[1])
                self.assertTrue(speech.call_args.kwargs['conductor'])
                agent.render_station_report(Path(directory) / 'out.wav', True)
                self.assertIn('Stress Test Station', speech.call_args.args[1])
                self.assertIn('express service, with reduced checks', speech.call_args.args[1])
                agent.args.fast = False
                agent.render_station_report(Path(directory) / 'out.wav', True)
                self.assertNotIn('express', speech.call_args.args[1])
                self.assertIn('regular testing schedule', speech.call_args.args[1])
                agent.args.pedantic = True
                agent.render_station_report(Path(directory) / 'out.wav', True)
                self.assertIn('expect significant delays', speech.call_args.args[1])

    def test_station_destination_choices_avoid_stop_and_early_arrival(self):
        """Each periodic destination option avoids alarm and premature arrival."""
        with tempfile.TemporaryDirectory() as directory, patch('sys.argv', ['tagent', '--state-dir', directory]):
            agent = tagent.Agent(tagent.arguments())
            for stress in (False, True):
                for index in range(4):
                    with patch.object(tagent.random, 'choice', side_effect=lambda choices: choices[index]), \
                         patch.object(agent, 'synthesize_announcement', return_value=1) as speech:
                        agent.render_station_report(Path(directory) / 'out.wav', stress)
                        message = speech.call_args.args[1].lower()
                        self.assertNotIn('next stop', message)
                        self.assertNotIn('approaching', message)
                        self.assertIn('station', message)

    def test_phase_change_keeps_train_timer(self):
        """Both phases announce while the periodic five-minute clock is preserved."""
        with tempfile.TemporaryDirectory() as directory, patch('sys.argv', [
                'tagent', '--sound-theme', 'steampunk', '--state-dir', directory]):
            agent = tagent.Agent(tagent.arguments())
            with patch.object(tagent.time, 'monotonic', return_value=100), \
                 patch.object(agent, 'sound') as sound:
                agent.start_test_sound('perf')
                agent.start_test_sound('stress')
                self.assertEqual(agent.next_train_at, 400)
                self.assertEqual([call.args[0] for call in sound.call_args_list],
                                 ['depart-perf', 'station-perf', 'depart-stress', 'station-stress'])

    def test_long_audio_wait_is_interruptible_and_not_cut_at_55_seconds(self):
        """A slow station/departure clip can finish; individual waits stay short."""
        with tempfile.TemporaryDirectory() as directory, patch('sys.argv', [
                'tagent', '--state-dir', directory]):
            agent = tagent.Agent(tagent.arguments())
            clock = [0.0]
            waits = []
            class LongClip:
                """Fake an audio thread that completes after 60 seconds."""
                def is_alive(self):
                    """Report whether simulated audio is still playing."""
                    return clock[0] < 60
                def join(self, timeout):
                    """Advance the simulated clock by the requested wait."""
                    waits.append(timeout)
                    clock[0] += timeout
            agent.audio_thread = LongClip()
            with patch.object(tagent.time, 'monotonic', side_effect=lambda: clock[0]):
                self.assertTrue(agent.wait_for_audio())
            self.assertEqual(clock[0], 60)
            self.assertLessEqual(max(waits), 1)

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

    def test_steam_demo_includes_mono_dispatch(self):
        """Demo plays each variant explicitly and combines mono clips correctly."""
        with tempfile.TemporaryDirectory() as directory, patch('sys.argv', [
                'tagent', '--sound-demo', '--sound-theme', 'steampunk',
                '--jobs', '8', '--state-dir', directory]):
            agent = tagent.Agent(tagent.arguments())
            agent.args.tts = 'off'
            def short_clip(path, theme, event, volume, workers, jobs, **kwargs):
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
                self.assertEqual(command.call_count, 19)
                with wave.open(str(Path(directory) / 'demo-steampunk.wav')) as stream:
                    self.assertEqual(stream.getnchannels(), 1)
                    self.assertGreater(stream.getnframes() / stream.getframerate(), 5)
                with wave.open(str(Path(directory) / 'sound-steampunk-dispatch.wav')) as stream:
                    self.assertEqual(stream.getnchannels(), 1)

    def test_bare_sound_preview_plays_all_assets(self):
        """A bare preview includes both backgrounds and produces a correctly timed demo."""
        with tempfile.TemporaryDirectory() as directory, patch('sys.argv', [
                'tagent', '--sound-preview', '--tts', 'off', '--state-dir', directory]):
            agent = tagent.Agent(tagent.arguments())
            events = []
            def short_clip(path, theme, event, volume, workers, jobs, **kwargs):
                with wave.open(str(path), 'wb') as stream:
                    stream.setparams((1, 2, soundtrack.RATE, 0, 'NONE', 'not compressed'))
                    stream.writeframes(b'\0' * (soundtrack.RATE // 10 * 2))
                return .1
            def background(path, volume, mode='depot'):
                events.append(mode)
                with wave.open(str(path), 'wb') as stream:
                    stream.setparams((1, 2, 22050, 0, 'NONE', 'not compressed'))
                    stream.writeframes(b'\0' * (22050 // 5 * 2))
                return .2
            with patch.object(soundtrack, 'render', side_effect=short_clip), \
                 patch.object(tagent.depot, 'render_preview', side_effect=background), \
                 patch.object(tagent.shutil, 'which', return_value='/fake/afplay'), \
                 patch.object(agent, 'command', return_value='') as command:
                self.assertEqual(agent.run(), 0)
                self.assertEqual(command.call_count, 21)
                self.assertEqual(events, ['depot', 'train'])
            with wave.open(str(Path(directory) / 'demo-steampunk.wav')) as stream:
                self.assertEqual(stream.getnchannels(), 1)
                # Two 0.2-second backgrounds keep their duration after upsampling.
                self.assertGreater(stream.getnframes() / stream.getframerate(), 7.5)

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
                self.assertEqual(played, ['depart-perf', 'station-perf'])
                clock[0] = 399
                agent.sound('work')
                self.assertEqual(played, ['depart-perf', 'station-perf'])
                agent.sound('build')
                agent.audio_thread.join()
                clock[0] = 400
                agent.sound('work')
                agent.audio_thread.join()
                self.assertEqual(played, ['depart-perf', 'station-perf', 'build', 'work'])
                clock[0] = 699
                agent.sound('work')
                self.assertEqual(played.count('work'), 1)
                clock[0] = 700
                agent.sound('work')
                agent.audio_thread.join()
                self.assertEqual(played, ['depart-perf', 'station-perf', 'build', 'work', 'dispatch'])
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

    def test_build_warning_counter_incremental(self):
        """Count real complete diagnostic lines once across logs, including colors."""
        with tempfile.TemporaryDirectory() as directory:
            perf = Path(directory) / 'perf.log'
            stress = Path(directory) / 'stress.log'
            perf.write_bytes(b'file.cpp:5:8: warning: one\n2 warnings generated.\n'
                             b'line | some warning: in source\nfile.cpp:9: war')
            counter = tagent.WarningCounter()
            self.assertEqual(counter.observe(perf), 1)
            self.assertEqual(counter.observe(perf), 1)
            with perf.open('ab') as stream:
                stream.write(b'ning: completed\nfile.cpp:1: \x1b[31mwarning:\x1b[0m three\n')
            self.assertEqual(counter.observe(perf), 3)
            stress.write_bytes(b'ld: warning: linker issue\nunit.cpp:7: fatal error: bad')
            self.assertEqual(counter.observe(stress), 4)
            self.assertEqual(counter.errors, 0)
            self.assertEqual(counter.observe(stress, final=True), 4)
            self.assertEqual(counter.errors, 1)
            counter.observe(stress, final=True)
            self.assertEqual(counter.errors, 1)

    def test_build_status_reports_actual_count_and_interval(self):
        """Dispatch uses accumulated counts; a summary is never double counted."""
        with tempfile.TemporaryDirectory() as directory, patch('sys.argv', [
                'tagent', '--sound-theme', 'steampunk', '--state-dir', directory]):
            agent = tagent.Agent(tagent.arguments())
            log = Path(directory) / 'build-perf.log'
            log.write_text('foo.cpp:1: warning: first\n1 warning generated.\n')
            agent.next_build_report = 0
            with patch.object(agent, 'sound') as sound, \
                 patch.object(tagent.time, 'monotonic', return_value=100):
                agent.build_progress(log)
                sound.assert_called_once_with('build-warning', workers=(1, 0, 0, 0, 0))
                agent.build_progress(log)
                self.assertEqual(sound.call_count, 1)
            self.assertEqual(agent.build_warnings.count, 1)
            self.assertIn('No compiler warnings so far', tagent.build_report_text(0))
            self.assertIn('1 compiler warning so far', tagent.build_report_text(1))
            self.assertIn('3 compiler warnings so far', tagent.build_report_text(3))
            self.assertIn('Sample builds: 2 warnings so far',
                          tagent.build_report_text(3, [1, 2, 3, 4, 2]))
            self.assertIn('Inverted list builds: 4 ready, 6 pending',
                          tagent.build_report_text(3, [1, 2, 3, 4, 2, 4, 6]))

    def test_build_report_tts_and_recorded_fallback(self):
        """Female voice receives an exact count; unavailable TTS uses safe fallback."""
        with tempfile.TemporaryDirectory() as directory, patch('sys.argv', [
                'tagent', '--sound-theme', 'steampunk', '--state-dir', directory]):
            agent = tagent.Agent(tagent.arguments())
            agent.tts_engine = '/fake/espeak-ng'
            agent.args.tts = 'auto'
            output = Path(directory) / 'report.wav'
            spoken = []
            def synth(argv, **kwargs):
                """Produce a short mono TTS WAV and capture voice options/text."""
                spoken.append(argv)
                voice = argv[argv.index('-w') + 1]
                with wave.open(voice, 'wb') as stream:
                    stream.setparams((1, 2, 22050, 0, 'NONE', 'not compressed'))
                    stream.writeframes(b'\x10\x10' * 4410)
            with patch.object(agent, 'command', side_effect=synth):
                agent.render_build_report(output, (3, 2, 4, 5, 7))
            self.assertIn(spoken[0][spoken[0].index('-v') + 1], ('en-us+f3', 'en-us+f4', 'en-us+m3'))
            self.assertIn('3 compiler warnings so far', spoken[0][-1])
            self.assertIn('Performance builds: 2 ready, 4 pending', spoken[0][-1])
            self.assertIn('Stress builds: 5 ready, 7 pending', spoken[0][-1])
            with patch.object(agent, 'command', side_effect=synth):
                agent.render_build_report(output, (3, 2, 4, 5, 7, 1, 4, 6))
            self.assertIn('Inverted list builds: 4 ready, 6 pending', spoken[-1][-1])
            with wave.open(str(output)) as stream:
                self.assertEqual(stream.getnchannels(), 1)
            agent.tts_engine = None
            with patch.object(soundtrack, 'mix_dynamic_dispatch', wraps=soundtrack.mix_dynamic_dispatch) as mix:
                agent.render_build_report(output, 0)
                self.assertEqual(Path(mix.call_args.args[0]).name, 'build-clear.wav')
                agent.render_build_report(output, 3)
                self.assertEqual(Path(mix.call_args.args[0]).name, 'build-warning.wav')

    def test_voice_defaults_and_explicit_overrides(self):
        """Voice models follow the log directory unless explicitly overridden."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            with patch.dict('os.environ', {'PIPER_MODEL': '', 'PIPER_CONDUCTOR_MODEL': '',
                                          'PIPER_ENGINEER_MODEL': ''}), \
                 patch('sys.argv', ['tagent', '--state-dir', str(root)]):
                agent = tagent.Agent(tagent.arguments())
            self.assertEqual(agent.piper_model, root / 'piper-voices/en_US-lessac-medium.onnx')
            self.assertEqual(agent.piper_conductor_model, root / 'piper-voices/en_GB-southern_english_female-low.onnx')
            self.assertEqual(agent.joker_model, root / 'piper-voices/en_US-danny-low.onnx')
            self.assertEqual(agent.temp_dir, root / 'tmp')
            self.assertEqual(agent.command([sys.executable, '-c',
                'import os; print(os.environ["TMPDIR"])']).strip(), str(agent.temp_dir))
            with patch('sys.argv', ['tagent', '--state-dir', str(root),
                    '--piper-model', str(root / 'custom.onnx'),
                    '--piper-engineer-model', str(root / 'engineer.onnx')]):
                custom = tagent.Agent(tagent.arguments())
            self.assertEqual(custom.piper_model, root / 'custom.onnx')
            self.assertEqual(custom.joker_model, root / 'engineer.onnx')

    def test_voice_downloads_use_log_directory(self):
        """All three downloaded voices and their configs land beneath LOG_DIR."""
        with tempfile.TemporaryDirectory() as directory, \
             patch.dict('os.environ', {'PIPER_MODEL': '', 'PIPER_CONDUCTOR_MODEL': '',
                                       'PIPER_ENGINEER_MODEL': ''}), \
             patch('sys.argv', ['tagent', '--state-dir', directory, '--piper-python', sys.executable]):
            agent = tagent.Agent(tagent.arguments())
            def command(argv, **kwargs):
                if 'piper.download_voices' in argv:
                    target = Path(argv[argv.index('--data-dir') + 1])
                    self.assertEqual(target.parent, agent.voices_dir)
                    (target / (argv[-1] + '.onnx')).write_bytes(b'model')
                    (target / (argv[-1] + '.onnx.json')).write_text('{}')
                else:
                    with wave.open(argv[argv.index('-f') + 1], 'wb') as output:
                        output.setparams((1, 2, 22050, 0, 'NONE', 'not compressed'))
                        output.writeframes(b'\0\0' * 32)
            with patch.object(agent, 'ensure_piper_engine'), patch.object(agent, 'command', side_effect=command):
                for model, kwargs in ((agent.piper_model, {}),
                                      (agent.piper_conductor_model, {'conductor': True}),
                                      (agent.joker_model, {'engineer': True})):
                    agent.download_piper_voice(model.stem, model, **kwargs)
                    self.assertTrue(model.is_file())
                    self.assertTrue(Path(str(model) + '.json').is_file())
            self.assertEqual(len(list(agent.voices_dir.iterdir())), 6)

    def test_interrupt_preview_waits_for_engineer_voice_with_joker_disabled(self):
        """Emergency preview waits for Danny even when jokes are disabled."""
        with tempfile.TemporaryDirectory() as directory, patch('sys.argv', [
                'tagent', '--state-dir', directory, '--sound-preview', 'interrupt', '--no-joker']):
            agent = tagent.Agent(tagent.arguments())
            ready = threading.Event()
            def download(voice, model, engineer=False):
                self.assertEqual(voice, 'en_US-danny-low')
                self.assertTrue(engineer)
                ready.wait(.05)
                ready.set()
            def sound(*args, **kwargs):
                self.assertTrue(ready.is_set(), 'preview started before engineer setup')
            with patch.object(agent, 'report_audio'), \
                 patch.object(agent, 'start_piper_download'), \
                 patch.object(agent, 'start_conductor_download'), \
                 patch.object(agent, 'download_piper_voice', side_effect=download), \
                 patch.object(agent, 'sound', side_effect=sound):
                self.assertEqual(agent.run(), 0)

    def test_installed_engineer_model_is_selected_for_interrupt(self):
        """Installed Danny is selected at startup and the actual engine is reported."""
        with tempfile.TemporaryDirectory() as directory:
            model = Path(directory) / 'en_US-danny-low.onnx'
            model.write_bytes(b'model')
            Path(str(model) + '.json').write_text('{}')
            primary = Path(directory) / 'en_US-lessac-medium.onnx'
            with patch('sys.argv', ['tagent', '--state-dir', directory,
                                    '--piper-python', sys.executable, '--piper-model', str(primary),
                                    '--piper-engineer-model', str(model)]):
                agent = tagent.Agent(tagent.arguments())
            self.assertTrue(agent.joker_available)
            with patch.object(agent, 'command') as command, \
                 patch.object(soundtrack, 'mix_dynamic_dispatch', return_value=1):
                agent.play_sound('interrupt-engineer', 1)
            self.assertIn(str(model), command.call_args_list[0].args[0])
            self.assertIn('TTS: Engineer using Piper', (Path(directory) / 'agent.log').read_text())

    def test_interrupt_exchange_stops_children_first(self):
        """Interrupt speech follows child cleanup and uses separate crew voices."""
        with tempfile.TemporaryDirectory() as directory, patch('sys.argv', [
                'tagent', '--state-dir', directory]):
            agent = tagent.Agent(tagent.arguments())
            agent.audio_available = True
            agent.stopping.set()
            events = []
            with patch.object(agent, 'stop_children', side_effect=lambda: events.append('stopped')), \
                 patch.object(agent, 'play_sound', side_effect=lambda event, workers: events.append(event)):
                agent.announce_interruption()
            self.assertEqual(events, ['stopped', 'interrupt-engineer', 'interrupt-dispatch'])
            self.assertTrue(agent.stopping.is_set())
            with patch.object(agent, 'synthesize_announcement', return_value=None) as speech:
                agent.play_sound('interrupt-engineer', 1)
                self.assertIn('Emergency brakes applied', speech.call_args.args[1])
                self.assertTrue(speech.call_args.kwargs['engineer'])
                agent.play_sound('interrupt-dispatch', 1)
                self.assertIn('Brakes acknowledged', speech.call_args.args[1])
                self.assertNotIn('engineer', speech.call_args.kwargs)

    def test_distress_sos_and_radio_exchange(self):
        """SOS precedes engineer and dispatch; disabled TTS keeps the SOS."""
        with tempfile.TemporaryDirectory() as directory, patch('sys.argv', [
                'tagent', '--state-dir', directory]):
            agent = tagent.Agent(tagent.arguments())
            agent.audio_available = True
            with patch.object(agent, 'play_sound', return_value=True) as play:
                agent.play_after(None, 'fail', 1, agent.audio_epoch)
                self.assertEqual([call.args[0] for call in play.call_args_list],
                                 ['fail', 'distress-engineer', 'distress-dispatch'])
                play.reset_mock()
                agent.args.tts = 'off'
                agent.play_after(None, 'fail', 1, agent.audio_epoch)
                self.assertEqual([call.args[0] for call in play.call_args_list], ['fail'])
            agent.args.tts = 'auto'
            with patch.object(agent, 'synthesize_announcement', return_value=None) as synth:
                agent.play_sound('distress-engineer', 1)
                self.assertIn('Mayday, mayday, mayday', synth.call_args.args[1])
                self.assertTrue(synth.call_args.kwargs['engineer'])
                agent.play_sound('distress-dispatch', 1)
                self.assertIn('Emergency response team is on the way', synth.call_args.args[1])
                self.assertNotIn('engineer', synth.call_args.kwargs)
                self.assertTrue(synth.call_args.kwargs['distress'])
            agent.args.sound_volume = .4
            with patch.object(soundtrack, 'render', return_value=1) as render, \
                 patch.object(tagent.shutil, 'which', return_value='/fake/player'), \
                 patch.object(agent, 'command'):
                agent.play_sound('fail', 1)
                self.assertAlmostEqual(render.call_args.args[3], .5)
                agent.args.sound_volume = 1.0
                agent.play_sound('fail', 1)
                self.assertEqual(render.call_args.args[3], 1.0)

    def test_joker_reply_probability_and_shuffled_deck(self):
        """A quarter of jokes have no reply; replies exhaust the deck before reuse."""
        with tempfile.TemporaryDirectory() as directory, patch('sys.argv', ['tagent', '--state-dir', directory]):
            agent = tagent.Agent(tagent.arguments())
            with patch.object(agent, 'play_sound', return_value=True), \
                 patch.object(tagent.threading, 'Thread') as thread, \
                 patch.object(tagent.random, 'random', return_value=.24):
                agent.play_after(None, 'joker', 'joke', agent.audio_epoch)
                thread.assert_not_called()
            with patch.object(agent, 'play_sound', return_value=True), \
                 patch.object(tagent.threading, 'Thread') as thread, \
                 patch.object(tagent.random, 'random', return_value=.25):
                agent.play_after(None, 'joker', 'joke', agent.audio_epoch)
                thread.assert_called_once()
            agent.audio_thread = None
            replies = []
            with patch.object(tagent.random, 'uniform', return_value=0), \
                 patch.object(agent, 'play_after') as play:
                for _ in tagent.DISPATCH_JOKER_REPLIES:
                    agent.queue_joker_reply(agent.audio_epoch)
                    agent.audio_thread.join(timeout=2)
                    replies.append(play.call_args.args[2])
            self.assertEqual(set(replies), set(tagent.DISPATCH_JOKER_REPLIES))
            self.assertEqual(len(replies), len(set(replies)))

    def test_joker_reply_scheduling_and_cancellation(self):
        """Dispatch responds after the pause unless a failure supersedes it."""
        with tempfile.TemporaryDirectory() as directory, patch('sys.argv', [
                'tagent', '--state-dir', directory]):
            agent = tagent.Agent(tagent.arguments())
            with patch.object(tagent.random, 'uniform', return_value=0) as delay, \
                 patch.object(agent, 'play_after') as play:
                agent.queue_joker_reply(agent.audio_epoch)
                agent.audio_thread.join(timeout=2)
                delay.assert_called_once_with(3, 7)
                self.assertEqual(play.call_args.args[1], 'joker-reply')
                self.assertIn(play.call_args.args[2], tagent.DISPATCH_JOKER_REPLIES)
                play.reset_mock()
                agent.audio_epoch += 1
                agent.queue_joker_reply(agent.audio_epoch - 1)
                play.assert_not_called()

    def test_repo_report_tts_and_fallback(self):
        """Remote outcomes have distinct speech and work without an installed TTS engine."""
        with tempfile.TemporaryDirectory() as directory, patch('sys.argv', [
                'tagent', '--state-dir', directory]):
            agent = tagent.Agent(tagent.arguments())
            agent.tts_engine = '/fake/espeak-ng'
            agent.args.tts = 'auto'
            output = Path(directory) / 'report.wav'
            with patch.object(agent, 'command') as command, \
                 patch.object(soundtrack, 'mix_dynamic_dispatch', return_value=1):
                agent.render_repo_report(output, True)
                self.assertIn('updates detected', command.call_args.args[0][-1])
                agent.render_repo_report(output, False)
                self.assertIn('No new updates', command.call_args.args[0][-1])
                agent.render_repo_report(output, True, redirect=True)
                self.assertIn('Freight train diverted to the build depot', command.call_args.args[0][-1])
                agent.render_repo_report(output, True, restart=True)
                self.assertIn('Dispatch to BitMagic. Updates received. Restarting TAgent', command.call_args.args[0][-1])
            agent.tts_engine = None
            with patch.object(soundtrack, 'render', return_value=1) as render:
                agent.render_repo_report(output, True)
                self.assertEqual(render.call_args.args[2], 'build')
                agent.render_repo_report(output, False)
                self.assertEqual(render.call_args.args[2], 'idle')

    def test_piper_preferred_and_falls_back_to_espeak(self):
        """The configured Piper environment is preferred; errors fall back to eSpeak."""
        with tempfile.TemporaryDirectory() as directory:
            model = Path(directory) / 'voice.onnx'
            model.write_bytes(b'model')
            Path(str(model) + '.json').write_text('{}')
            with patch('sys.argv', ['tagent', '--state-dir', directory, '--tts', 'auto',
                                     '--piper-python', sys.executable, '--piper-model', str(model)]):
                agent = tagent.Agent(tagent.arguments())
            self.assertTrue(agent.piper_available)
            agent.piper_conductor_model = model
            agent.piper_conductor_available = True
            agent.tts_engine = '/fake/espeak-ng'
            calls = []
            def synth(argv, **kwargs):
                """Produce valid mono audio and record the engine used."""
                calls.append(argv)
                flag = '-f' if '-f' in argv else '-w'
                with wave.open(argv[argv.index(flag) + 1], 'wb') as stream:
                    stream.setparams((1, 2, 22050, 0, 'NONE', 'not compressed'))
                    stream.writeframes(b'\x10\x10' * 2205)
            output = Path(directory) / 'report.wav'
            with patch.object(agent, 'command', side_effect=synth):
                agent.render_perf_progress(output, ('perf_release_64', 12))
            self.assertEqual(calls[0][:3], [str(agent.piper_python), '-m', 'piper'])
            self.assertEqual(calls[0][calls[0].index('-f') - 1], str(model))
            self.assertIn('perf 64 completed 12 tests', calls[0][-1])
            calls.clear()
            def broken_piper(argv, **kwargs):
                if argv[0] == str(agent.piper_python):
                    raise tagent.CommandFailure('Piper failed')
                return synth(argv, **kwargs)
            with patch.object(agent, 'command', side_effect=broken_piper):
                agent.render_perf_progress(output, ('perf_release_64', 12))
            self.assertTrue(agent.piper_conductor_failed)
            self.assertFalse(agent.piper_failed)
            self.assertEqual(calls[0][0], '/fake/espeak-ng')
            self.assertIn('falling back to eSpeak', (Path(directory) / 'agent.log').read_text())
            with wave.open(str(output)) as stream:
                self.assertEqual((stream.getnchannels(), stream.getframerate()), (1, 44100))
            agent.piper_conductor_failed = False
            calls.clear()
            def corrupt_piper(argv, **kwargs):
                if argv[0] == str(agent.piper_python):
                    Path(argv[argv.index('-f') + 1]).write_bytes(b'invalid WAV')
                else:
                    synth(argv, **kwargs)
            with patch.object(agent, 'command', side_effect=corrupt_piper):
                agent.render_perf_progress(output, ('perf_release_64', 12))
            self.assertTrue(agent.piper_conductor_failed)
            self.assertEqual(calls[0][0], '/fake/espeak-ng')
            agent.piper_conductor_failed = False
            agent.tts_engine = None
            expected = Path(directory) / 'recorded.wav'
            soundtrack.render_dispatch(expected, agent.args.sound_volume)
            with patch.object(agent, 'command', side_effect=tagent.CommandFailure('Piper failed')):
                agent.render_dispatch(output, 4)
            self.assertEqual(output.read_bytes(), expected.read_bytes())

    def test_speech_gain_preserves_background(self):
        """Lowering voice gain changes speech by 15% without changing the bed."""
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'voice.wav'
            with wave.open(str(source), 'wb') as stream:
                stream.setparams((1, 2, 22050, 0, 'NONE', 'not compressed'))
                stream.writeframes(array.array('h', [10000] * 22050).tobytes())
            samples = []
            for gain in (0, 1, .85):
                path = Path(directory) / (str(gain) + '.wav')
                soundtrack.mix_dynamic_dispatch(source, path, .35, 1, 1, voice_gain=gain)
                with wave.open(str(path), 'rb') as stream:
                    samples.append(array.array('h', stream.readframes(stream.getnframes())))
            index = 44100 // 2
            full = samples[1][index] - samples[0][index]
            quiet = samples[2][index] - samples[0][index]
            self.assertAlmostEqual(quiet / full, .85, delta=.001)
            self.assertEqual(samples[0][0], samples[1][0])

    def test_piper_requires_model_configuration(self):
        """An incomplete voice installation is skipped without preventing startup."""
        with tempfile.TemporaryDirectory() as directory:
            model = Path(directory) / 'voice.onnx'
            model.write_bytes(b'model')
            with patch('sys.argv', ['tagent', '--state-dir', directory,
                                     '--piper-python', sys.executable, '--piper-model', str(model)]):
                agent = tagent.Agent(tagent.arguments())
            self.assertFalse(agent.piper_available)

    def test_piper_conductor_and_dispatch_use_separate_models(self):
        """Perf progress uses the UK model while dispatch retains the US model."""
        with tempfile.TemporaryDirectory() as directory:
            primary = Path(directory) / 'en_US-lessac-medium.onnx'
            conductor = Path(directory) / 'en_GB-southern_english_female-low.onnx'
            for model in (primary, conductor):
                model.write_bytes(b'model')
                Path(str(model) + '.json').write_text('{}')
            with patch('sys.argv', ['tagent', '--state-dir', directory, '--tts', 'auto',
                                    '--piper-python', sys.executable, '--piper-model', str(primary),
                                    '--piper-conductor-model', str(conductor)]):
                agent = tagent.Agent(tagent.arguments())
            with patch.object(agent, 'command') as command, \
                 patch.object(soundtrack, 'mix_dynamic_dispatch', return_value=1):
                agent.render_dispatch(Path(directory) / 'out.wav', 6)
                self.assertIn(str(primary), command.call_args.args[0])
                agent.render_perf_progress(Path(directory) / 'out.wav', ('perf_release_64', 12))
                self.assertIn(str(conductor), command.call_args.args[0])
                self.assertIn('TTS: Conductor using Piper', (Path(directory) / 'agent.log').read_text())
                agent.piper_conductor_available = False
                agent.tts_engine = '/fake/espeak-ng'
                agent.render_perf_progress(Path(directory) / 'out.wav', ('perf_release_64', 12))
                self.assertEqual(command.call_args.args[0][0], '/fake/espeak-ng')
                self.assertIn('TTS: Conductor using eSpeak NG fallback', (Path(directory) / 'agent.log').read_text())

    def test_missing_conductor_voice_downloads_separately(self):
        """The conductor model is downloaded without replacing the dispatch voice."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            primary = root / 'en_US-lessac-medium.onnx'
            conductor = root / 'en_GB-southern_english_female-low.onnx'
            primary.write_bytes(b'primary')
            Path(str(primary) + '.json').write_text('{}')
            with patch('sys.argv', ['tagent', '--state-dir', directory, '--tts', 'auto',
                                    '--piper-python', sys.executable, '--piper-model', str(primary),
                                    '--piper-conductor-model', str(conductor)]):
                agent = tagent.Agent(tagent.arguments())
            def command(argv, **kwargs):
                if '-c' in argv:
                    return ''
                if 'piper.download_voices' in argv:
                    staging = Path(argv[argv.index('--data-dir') + 1])
                    self.assertEqual(argv[-1], conductor.stem)
                    (staging / conductor.name).write_bytes(b'conductor')
                    (staging / (conductor.name + '.json')).write_text('{}')
                else:
                    with wave.open(argv[argv.index('-f') + 1], 'wb') as stream:
                        stream.setparams((1, 2, 22050, 0, 'NONE', 'not compressed'))
                        stream.writeframes(b'\0\0' * 2205)
            with patch.object(agent, 'command', side_effect=command):
                agent.start_conductor_download()
                agent.conductor_download_thread.join(5)
            self.assertTrue(agent.piper_conductor_available)
            self.assertEqual(conductor.read_bytes(), b'conductor')
            self.assertEqual(primary.read_bytes(), b'primary')

    def test_conductor_setup_starts_before_builds(self):
        """Normal startup requests both voice preparations before validation/builds."""
        with tempfile.TemporaryDirectory() as directory, patch('sys.argv', [
                'tagent', '--state-dir', directory, '--tts', 'auto', '--no-joker']):
            agent = tagent.Agent(tagent.arguments())
            with patch.object(agent, 'report_audio'), \
                 patch.object(agent, 'start_piper_download') as primary, \
                 patch.object(agent, 'start_conductor_download') as conductor, \
                 patch.object(agent, 'download_piper_voice'), \
                 patch.object(agent, 'validate', side_effect=tagent.Failure('stop')):
                with self.assertRaisesRegex(tagent.Failure, 'stop'):
                    agent.run()
            primary.assert_called_once_with()
            conductor.assert_called_once_with()

    def test_conductor_download_retries_after_failure(self):
        """A failed startup download can be retried when the perf phase begins."""
        with tempfile.TemporaryDirectory() as directory, patch('sys.argv', [
                'tagent', '--state-dir', directory, '--tts', 'auto',
                '--piper-python', sys.executable,
                '--piper-conductor-model', str(Path(directory) / 'en_GB-southern_english_female-low.onnx')]):
            agent = tagent.Agent(tagent.arguments())
            with patch.object(agent, 'download_piper_voice') as download:
                agent.start_conductor_download()
                agent.conductor_download_thread.join(5)
                agent.start_conductor_download()
                agent.conductor_download_thread.join(5)
                self.assertEqual(download.call_count, 2)

    def test_missing_piper_voice_downloads_and_validates_in_background(self):
        """A missing voice is staged and verified without blocking the caller."""
        import threading
        with tempfile.TemporaryDirectory() as directory:
            model = Path(directory) / 'voices' / 'en_US-lessac-medium.onnx'
            with patch('sys.argv', ['tagent', '--state-dir', directory, '--tts', 'auto',
                                     '--piper-python', sys.executable, '--piper-model', str(model)]):
                agent = tagent.Agent(tagent.arguments())
            started, release = threading.Event(), threading.Event()
            def command(argv, **kwargs):
                """Emulate the downloader and Piper's WAV validation probe."""
                if '-c' in argv:
                    return ''
                if 'piper.download_voices' in argv:
                    started.set()
                    release.wait(5)
                    staging = Path(argv[argv.index('--data-dir') + 1])
                    (staging / model.name).write_bytes(b'complete model')
                    (staging / (model.name + '.json')).write_text('{}')
                else:
                    with wave.open(argv[argv.index('-f') + 1], 'wb') as stream:
                        stream.setparams((1, 2, 22050, 0, 'NONE', 'not compressed'))
                        stream.writeframes(b'\0\0' * 2205)
            try:
                with patch.object(agent, 'command', side_effect=command):
                    agent.start_piper_download()
                    self.assertTrue(started.wait(2))
                    self.assertFalse(agent.piper_available)
                    release.set()
                    agent.piper_download_thread.join(5)
                    self.assertFalse(agent.piper_download_thread.is_alive())
                self.assertTrue(agent.piper_available)
                self.assertEqual(model.read_bytes(), b'complete model')
                self.assertTrue(Path(str(model) + '.json').is_file())
            finally:
                release.set()

    def test_piper_voice_download_failure_preserves_fallback(self):
        """A partial download is discarded and eSpeak remains available."""
        with tempfile.TemporaryDirectory() as directory:
            model = Path(directory) / 'voices' / 'en_US-lessac-medium.onnx'
            with patch('sys.argv', ['tagent', '--state-dir', directory, '--tts', 'auto',
                                     '--piper-python', sys.executable, '--piper-model', str(model)]):
                agent = tagent.Agent(tagent.arguments())
            agent.tts_engine = '/fake/espeak-ng'
            def failed_download(argv, **kwargs):
                if '-c' in argv:
                    return ''
                staging = Path(argv[argv.index('--data-dir') + 1])
                (staging / model.name).write_bytes(b'partial')
                raise tagent.CommandFailure('offline')
            with patch.object(agent, 'command', side_effect=failed_download):
                agent.start_piper_download()
                agent.piper_download_thread.join(5)
            self.assertFalse(agent.piper_available)
            self.assertFalse(model.exists())
            self.assertEqual(agent.tts_engine, '/fake/espeak-ng')
            self.assertIn('continuing with TTS fallback', (Path(directory) / 'agent.log').read_text())

    def test_piper_voice_download_can_be_disabled(self):
        """A user can opt out of downloading missing voices."""
        with tempfile.TemporaryDirectory() as directory, patch('sys.argv', [
                'tagent', '--state-dir', directory, '--piper-python', sys.executable,
                '--piper-model', str(Path(directory) / 'en_US-lessac-medium.onnx'),
                '--no-piper-download']):
            agent = tagent.Agent(tagent.arguments())
            agent.start_piper_download()
            self.assertIsNone(agent.piper_download_thread)

    def test_piper_bootstrap_creates_user_venv_and_installs_voice(self):
        """A fresh target is provisioned in its own venv without system installs."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / 'piper-demo'
            venv = root / '.venv'
            python = venv / 'bin' / 'python'
            model = root / 'en_US-lessac-medium.onnx'
            with patch('sys.argv', ['tagent', '--state-dir', directory, '--tts', 'auto',
                                     '--piper-python', str(python), '--piper-model', str(model)]):
                agent = tagent.Agent(tagent.arguments())
            calls, installed = [], [False]
            def command(argv, **kwargs):
                calls.append(argv)
                if 'venv' in argv:
                    python.parent.mkdir(parents=True)
                    python.write_text('#!/bin/sh\nexit 0\n')
                    python.chmod(0o755)
                    (venv / 'pyvenv.cfg').write_text('home = /tmp\n')
                elif '-c' in argv:
                    if not installed[0]:
                        raise tagent.CommandFailure('module missing')
                elif 'pip' in argv:
                    installed[0] = True
                elif 'piper.download_voices' in argv:
                    staging = Path(argv[argv.index('--data-dir') + 1])
                    (staging / model.name).write_bytes(b'model')
                    (staging / (model.name + '.json')).write_text('{}')
                else:
                    with wave.open(argv[argv.index('-f') + 1], 'wb') as stream:
                        stream.setparams((1, 2, 22050, 0, 'NONE', 'not compressed'))
                        stream.writeframes(b'\0\0' * 2205)
            with patch.object(agent, 'command', side_effect=command):
                agent.start_piper_download()
                agent.piper_download_thread.join(5)
            self.assertTrue(agent.piper_available)
            self.assertTrue(any('venv' in argv for argv in calls))
            install = next(argv for argv in calls if 'pip' in argv)
            self.assertEqual(install[0], str(python))
            self.assertIn('piper-tts', install)

    def test_piper_bootstrap_does_not_install_into_system_python(self):
        """A configured interpreter outside a venv is never modified."""
        with tempfile.TemporaryDirectory() as directory:
            python = Path(directory) / 'system' / 'bin' / 'python'
            python.parent.mkdir(parents=True)
            python.write_text('#!/bin/sh\nexit 0\n')
            python.chmod(0o755)
            with patch('sys.argv', ['tagent', '--state-dir', directory, '--piper-python', str(python)]):
                agent = tagent.Agent(tagent.arguments())
            with patch.object(agent, 'command', side_effect=tagent.CommandFailure('missing')) as command:
                with self.assertRaisesRegex(tagent.Failure, 'virtual environment'):
                    agent.ensure_piper_engine()
            self.assertEqual(command.call_count, 1)

    def test_piper_bootstrap_install_failure_uses_fallback(self):
        """A failed package installation cannot block or disable eSpeak fallback."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / 'piper-demo'
            python = root / '.venv' / 'bin' / 'python'
            model = root / 'en_US-lessac-medium.onnx'
            with patch('sys.argv', ['tagent', '--state-dir', directory, '--tts', 'auto',
                                     '--piper-python', str(python), '--piper-model', str(model)]):
                agent = tagent.Agent(tagent.arguments())
            agent.tts_engine = '/fake/espeak-ng'
            def command(argv, **kwargs):
                if 'venv' in argv:
                    python.parent.mkdir(parents=True)
                    python.write_text('#!/bin/sh\nexit 0\n')
                    python.chmod(0o755)
                    (root / '.venv' / 'pyvenv.cfg').write_text('home = /tmp\n')
                else:
                    raise tagent.CommandFailure('offline or missing module')
            with patch.object(agent, 'command', side_effect=command):
                agent.start_piper_download()
                agent.piper_download_thread.join(5)
            self.assertFalse(agent.piper_available)
            self.assertEqual(agent.tts_engine, '/fake/espeak-ng')
            self.assertFalse(model.exists())
            self.assertIn('continuing with TTS fallback', (Path(directory) / 'agent.log').read_text())

    def test_perf_progress_female_tts_and_fallback(self):
        """The conductor speaks the chosen variant/count or uses a neutral sound."""
        with tempfile.TemporaryDirectory() as directory, patch('sys.argv', [
                'tagent', '--state-dir', directory]):
            agent = tagent.Agent(tagent.arguments())
            agent.args.tts = 'auto'
            agent.args.dispatcher_voice = 'en-us+m3,en-us+f4'
            agent.tts_engine = '/fake/espeak-ng'
            with patch.object(agent, 'command') as command, \
                 patch.object(soundtrack, 'mix_dynamic_dispatch', return_value=1):
                agent.render_perf_progress(Path(directory) / 'out.wav', ('perf_release_avx2', 63))
                argv = command.call_args.args[0]
                self.assertEqual(argv[argv.index('-v') + 1], 'en-us+f4')
                self.assertIn('perf avx2 completed 63 tests', argv[-1])
                self.assertIn('Approximately 21 percent of test elements complete', argv[-1])
                self.assertEqual(tagent.perf_estimated_percent(299), 99)
                self.assertEqual(tagent.perf_estimated_percent(350), 99)
                agent.render_perf_progress(Path(directory) / 'out.wav', ('perf_release_avx2', 350))
                self.assertIn('perf avx2 train is delayed. 350 tests completed', command.call_args.args[0][-1])
                agent.render_perf_progress(Path(directory) / 'out.wav', ('perf_release_64', 12))
                self.assertIn('perf 64 completed 12 tests', command.call_args.args[0][-1])
            agent.tts_engine = None
            with patch.object(soundtrack, 'render', return_value=1) as render:
                agent.render_perf_progress(Path(directory) / 'out.wav', ('perf_release', 2))
                self.assertEqual(render.call_args.args[2], 'start')

    def test_build_final_tts_and_exact_recorded_numbers(self):
        """Success/failure and both diagnostic totals survive online and offline paths."""
        import json
        with tempfile.TemporaryDirectory() as directory, patch('sys.argv', [
                'tagent', '--sound-theme', 'steampunk', '--state-dir', directory]):
            agent = tagent.Agent(tagent.arguments())
            agent.tts_engine = '/fake/espeak-ng'
            output = Path(directory) / 'summary.wav'
            spoken = []
            def synth(argv, **kwargs):
                """Create a short TTS result and record its exact announcement."""
                spoken.append(argv)
                with wave.open(argv[argv.index('-w') + 1], 'wb') as stream:
                    stream.setparams((1, 2, 22050, 0, 'NONE', 'not compressed'))
                    stream.writeframes(b'\x10\x10' * 4410)
            with patch.object(agent, 'command', side_effect=synth):
                agent.render_build_summary(output, (12, 34, False))
            self.assertIn('Build failed. 12 compiler errors. 34 compiler warnings.', spoken[0][-1])
            self.assertIn(spoken[0][spoken[0].index('-v') + 1], ('en-us+f3', 'en-us+f4', 'en-us+m3'))
            agent.tts_engine = None
            collected = []
            def capture(voice, path, volume, workers, budget, background, radio=False):
                """Capture assembled offline speech before its temporary file is removed."""
                with wave.open(str(voice), 'rb') as stream:
                    collected.append(stream.readframes(stream.getnframes()))
                return 1
            with patch.object(soundtrack, 'mix_dynamic_dispatch', side_effect=capture):
                agent.render_build_summary(output, (12, 34, False))
            assets = Path(soundtrack.__file__).parent / 'assets'
            meta = json.loads((assets / 'build-summary-voice.json').read_text())
            with wave.open(str(assets / 'build-summary-voice.wav'), 'rb') as stream:
                position = 0
                for word in ('failed', 'errors', 'one', 'two', 'warnings', 'three', 'four'):
                    start, length = meta[word]
                    stream.setpos(start)
                    fragment = stream.readframes(length)
                    position = collected[0].find(fragment, position)
                    self.assertGreaterEqual(position, 0, word)
                    position += len(fragment)

    def test_dynamic_dispatch_uses_observed_count(self):
        """A real fake TTS executable produces mono speech and receives the actual count."""
        with tempfile.TemporaryDirectory() as directory, patch('sys.argv', [
                'tagent', '--sound-theme', 'steampunk', '--state-dir', directory]):
            agent = tagent.Agent(tagent.arguments())
            engine = Path(directory) / 'fake-espeak'
            message = Path(directory) / 'spoken.txt'
            engine.write_text('#!' + sys.executable + "\nimport sys,wave\n" +
                'open(' + repr(str(message)) + ',"w").write(sys.argv[-1])\n' +
                "path=sys.argv[sys.argv.index('-w')+1]\n" +
                "with wave.open(path,'wb') as w:\n w.setparams((1,2,22050,0,'NONE','not compressed')); w.writeframes(b'\\x10\\x10'*4410)\n")
            engine.chmod(0o755)
            agent.tts_engine = str(engine)
            output = Path(directory) / 'out.wav'
            agent.render_dispatch(output, 6)
            self.assertEqual(message.read_text(), 'Dispatch to BitMagic. Freight rolling. 6 test crews at work.')
            with wave.open(str(output), 'rb') as stream:
                self.assertEqual(stream.getparams()[:3], (1, 2, 44100))
                self.assertGreater(stream.getnframes(), 44100)
            agent.render_dispatch(output, 1)
            self.assertIn('1 test crew at work', message.read_text())

    def test_tts_missing_disabled_or_failed_uses_recording(self):
        """Unavailable/disabled/broken TTS always falls back to the exact bundled WAV."""
        with tempfile.TemporaryDirectory() as directory, patch('sys.argv', [
                'tagent', '--sound-theme', 'steampunk', '--state-dir', directory]):
            agent = tagent.Agent(tagent.arguments())
            agent.args.sound_volume = .5
            output = Path(directory) / 'actual.wav'
            expected = Path(directory) / 'expected.wav'
            soundtrack.render_dispatch(expected, .5)
            agent.tts_engine = None
            agent.render_dispatch(output, 4)
            self.assertEqual(output.read_bytes(), expected.read_bytes())
            self.assertIn('unavailable', (Path(directory) / 'agent.log').read_text())
            agent.args.tts = 'off'
            agent.tts_engine = '/nonexistent/espeak-ng'
            agent.render_dispatch(output, 4)
            self.assertEqual(output.read_bytes(), expected.read_bytes())
            agent.args.tts = 'auto'
            agent.tts_warned = False
            broken = Path(directory) / 'broken-tts'
            broken.write_text('#!/bin/sh\nexit 7\n')
            broken.chmod(0o755)
            agent.tts_engine = str(broken)
            agent.render_dispatch(output, 4)
            self.assertEqual(output.read_bytes(), expected.read_bytes())
            self.assertIn('failed', (Path(directory) / 'agent.log').read_text())
            empty = Path(directory) / 'empty-tts'
            empty.write_text("#!/bin/sh\nwhile [ \"$1\" != -w ]; do shift; done\nshift\n: > \"$1\"\n")
            empty.chmod(0o755)
            agent.tts_engine = str(empty)
            agent.render_dispatch(output, 4)
            self.assertEqual(output.read_bytes(), expected.read_bytes())

    def test_missing_audio_is_nonfatal_and_reported(self):
        """No player and no writable console leaves the agent running with a notice."""
        with tempfile.TemporaryDirectory() as directory, patch('sys.argv', [
                'tagent', '--sound-theme', 'steampunk', '--tts', 'off', '--state-dir', directory]):
            agent = tagent.Agent(tagent.arguments())
            agent.prepare_generated_audio()
            import builtins
            original_open = builtins.open
            def no_console(path, *args, **kwargs):
                """Deny only the console; leave WAV and log I/O operational."""
                if str(path) == '/dev/console':
                    raise PermissionError('no console')
                return original_open(path, *args, **kwargs)
            with patch.object(tagent.shutil, 'which', return_value=None), \
                 patch.object(tagent.os, 'access', return_value=True), \
                 patch('builtins.open', side_effect=no_console):
                agent.report_audio()
                agent.play_sound('dispatch', 1)
            log = (Path(directory) / 'agent.log').read_text()
            self.assertIn('No audio player found', log)
            self.assertIn('audio unavailable', log)
            self.assertTrue(agent.audio_available)

    def test_player_failure_disables_sound_but_not_agent(self):
        """A present player with no working device falls back and allows a later retry."""
        import builtins
        with tempfile.TemporaryDirectory() as directory, patch('sys.argv', [
                'tagent', '--sound-theme', 'steampunk', '--tts', 'off', '--state-dir', directory]):
            agent = tagent.Agent(tagent.arguments())
            agent.prepare_generated_audio()
            original_open = builtins.open
            def no_console(path, *args, **kwargs):
                """Reject only the console device, preserving file access."""
                if str(path) == '/dev/console':
                    raise PermissionError('no console')
                return original_open(path, *args, **kwargs)
            with patch.object(tagent.shutil, 'which', side_effect=lambda name: '/fake/afplay' if name == 'afplay' else None), \
                 patch.object(agent, 'command', side_effect=tagent.CommandFailure('no device')), \
                 patch('builtins.open', side_effect=no_console):
                agent.play_sound('dispatch', 1)
            self.assertTrue(agent.audio_available)
            log = (Path(directory) / 'agent.log').read_text()
            self.assertIn('No working audio player', log)
            self.assertIn('audio unavailable', log)
            self.assertTrue(agent.audio_available)

    def test_launcher_detects_and_uses_local_tts(self):
        """The actual shell launcher detects eSpeak and previews six active crews."""
        import os
        import subprocess
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            binary = root / 'bin'
            binary.mkdir()
            engine = binary / 'espeak-ng'
            engine.write_text('#!' + sys.executable + "\nimport os,sys,wave\n" +
                "open(os.environ['TTS_RECORD'],'w').write(sys.argv[-1])\n" +
                "path=sys.argv[sys.argv.index('-w')+1]\n" +
                "with wave.open(path,'wb') as w:\n w.setparams((1,2,22050,0,'NONE','not compressed')); w.writeframes(b'\\x10\\x10'*4410)\n")
            engine.chmod(0o755)
            player = binary / 'afplay'
            player.write_text('#!/bin/sh\nexit 0\n')
            player.chmod(0o755)
            env = dict(os.environ, PATH=str(binary) + os.pathsep + os.environ['PATH'],
                       PYTHON=sys.executable, TTS_MODE='auto', SOUND='yes',
                       SOUND_VOLUME='0.35', STARTUP='build', TTS_RECORD=str(root / 'spoken'))
            script = Path(__file__).with_name('run_tagent.sh')
            result = subprocess.run(['bash', str(script), '--sound-theme', 'steampunk',
                                     '--jobs', '6', '--state-dir', str(root / 'state'),
                                     '--sound-preview', 'dispatch'], env=env, capture_output=True,
                                    text=True, timeout=20)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn('TTS: eSpeak NG available', result.stdout)
            self.assertIn('6 test crews at work', (root / 'spoken').read_text())
            with wave.open(str(root / 'state/sound-steampunk-dispatch.wav')) as stream:
                self.assertEqual(stream.getnchannels(), 1)
                self.assertEqual(stream.getframerate(), 44100)

    def test_no_audio_preview_still_exports_wav(self):
        """Headless preview writes the fallback WAV while reporting no sound device."""
        with tempfile.TemporaryDirectory() as directory, patch('sys.argv', [
                'tagent', '--sound-theme', 'steampunk', '--tts', 'off',
                '--sound-preview', 'dispatch', '--state-dir', directory]):
            agent = tagent.Agent(tagent.arguments())
            with patch.object(tagent.shutil, 'which', return_value=None), \
                 patch.object(tagent.os, 'access', return_value=False):
                self.assertEqual(agent.run(), 0)
                agent.audio_thread.join(timeout=10)
                self.assertFalse(agent.audio_thread.is_alive())
                self.assertTrue(agent.audio_available)
                self.assertTrue((Path(directory) / 'generated-audio' / 'sound-steampunk-dispatch.wav').exists())
                self.assertIn('playback will be retried', (Path(directory) / 'agent.log').read_text())

    def test_dispatch_preview_requires_steampunk(self):
        """Removed dieselpunk theme is rejected."""
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

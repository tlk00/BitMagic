"""Checks for the optional, interruptible background audio worker."""
import array
import subprocess
import tempfile
import threading
import time
import unittest
import wave
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import depot
import generate_depot_ambient


class DepotTests(unittest.TestCase):
    """Exercise fades and failed players without a physical sound device."""

    def test_sherman_passage_has_no_scripted_gaps(self):
        """A steady source remains audible throughout the generated depot loop."""
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'source.wav'
            result = Path(directory) / 'result.wav'
            with wave.open(str(source), 'wb') as stream:
                stream.setparams((1, 2, 22050, 0, 'NONE', 'not compressed'))
                stream.writeframes(b'\0\0' * (670 * 22050))
                stream.writeframes(b'\x00\x10' * (65 * 22050))
            generate_depot_ambient.render(source, result)
            with wave.open(str(result), 'rb') as stream:
                self.assertEqual(stream.getnchannels(), 1)
                for second in range(60):
                    chunk = array.array('h', stream.readframes(22050))
                    self.assertGreater(sum(abs(value) for value in chunk) / len(chunk), 700,
                                       'unexpected quiet at second {}'.format(second))

    def test_preview_renders_mono_loop(self):
        """The standalone preview is a scaled minute of portable PCM audio."""
        with tempfile.TemporaryDirectory() as directory:
            for mode, seconds in (('depot', 60), ('train', 300)):
                path = Path(directory) / (mode + '.wav')
                self.assertEqual(depot.render_preview(path, .35, mode), seconds)
                with wave.open(str(path)) as stream:
                    self.assertEqual((stream.getnchannels(), stream.getsampwidth()), (1, 2))
                    self.assertEqual(stream.getnframes(), seconds * stream.getframerate())

    def test_switch_modes_and_resume_after_test_interruption(self):
        """The worker changes loops and survives a temporary test-stop signal."""
        with tempfile.TemporaryDirectory() as directory:
            played = threading.Event()
            loaded = threading.Event()
            original = depot.read_loop
            def read_loop(mode):
                """Observe the train loop being selected."""
                result = original(mode)
                if mode == 'train':
                    loaded.set()
                return result
            def command(argv, timeout):
                """Simulate short player calls."""
                played.set()
                time.sleep(.02)
            agent = SimpleNamespace(home=Path(directory), args=SimpleNamespace(sound_volume=.35),
                                    stopping=threading.Event(), command=command, say=lambda _: None)
            worker = depot.Depot(agent)
            with patch.object(depot.shutil, 'which', side_effect=lambda name: name == 'afplay'), \
                 patch.object(depot, 'read_loop', side_effect=read_loop):
                try:
                    worker.start()
                    self.assertTrue(played.wait(2))
                    worker.set_mode('train')
                    self.assertTrue(loaded.wait(3))
                    agent.stopping.set()
                    time.sleep(.15)
                    self.assertTrue(worker.thread.is_alive())
                    played.clear()
                    agent.stopping.clear()
                    self.assertTrue(played.wait(2))
                finally:
                    worker.close()

    def test_fades_to_silence_and_resumes(self):
        """A foreground clip waits until the background player has stopped."""
        with tempfile.TemporaryDirectory() as directory:
            plays = []
            edges = []
            played = threading.Event()
            def command(argv, timeout):
                """Observe each complete background chunk as a fake player."""
                self.assertTrue(Path(argv[1]).exists())
                with wave.open(argv[1], 'rb') as stream:
                    samples = array.array('h', stream.readframes(stream.getnframes()))
                    edges.append((samples[0], samples[-1]))
                plays.append(time.monotonic())
                played.set()
                time.sleep(.02)
            agent = SimpleNamespace(home=Path(directory), args=SimpleNamespace(sound_volume=.35),
                                    stopping=threading.Event(), command=command, say=lambda _: None)
            worker = depot.Depot(agent)
            with patch.object(depot.shutil, 'which', side_effect=lambda name: name == 'afplay'):
                try:
                    worker.start()
                    self.assertTrue(played.wait(2))
                    worker.mute()
                    self.assertTrue(worker.quiet.is_set())
                    count = len(plays)
                    time.sleep(.15)
                    self.assertEqual(len(plays), count)
                    played.clear()
                    worker.unmute()
                    self.assertTrue(played.wait(2))
                finally:
                    worker.close()
            self.assertFalse(worker.thread.is_alive())
            self.assertTrue(edges)
            self.assertTrue(all(first == 0 and last == 0 for first, last in edges))

    def test_player_failure_releases_foreground(self):
        """A broken audio device ends ambience and cannot block an announcement."""
        with tempfile.TemporaryDirectory() as directory:
            messages = []
            def command(argv, timeout):
                """Simulate an unavailable audio device."""
                raise RuntimeError('no device')
            agent = SimpleNamespace(home=Path(directory), args=SimpleNamespace(sound_volume=.35),
                                    stopping=threading.Event(), command=command, say=messages.append)
            worker = depot.Depot(agent)
            with patch.object(depot.shutil, 'which', side_effect=lambda name: name == 'afplay'):
                worker.start()
                worker.thread.join(3)
                worker.mute()
            self.assertFalse(worker.thread.is_alive())
            self.assertTrue(worker.quiet.is_set())
            self.assertTrue(any('continuing without background audio' in text for text in messages))

    def test_repeated_player_processes_keep_running(self):
        """A sequence of real player processes continues past three seconds."""
        with tempfile.TemporaryDirectory() as directory:
            script = Path(directory) / 'fake-player'
            script.write_text('#!/bin/sh\nsleep 1\n')
            script.chmod(0o755)
            calls = []
            fifth = threading.Event()
            def command(argv, timeout):
                """Run a one-second fake player and observe successive calls."""
                calls.append(time.monotonic())
                if len(calls) >= 5:
                    fifth.set()
                subprocess.run([str(script), argv[1]], check=True, timeout=timeout)
            agent = SimpleNamespace(home=Path(directory), args=SimpleNamespace(sound_volume=.35),
                                    stopping=threading.Event(), command=command, say=lambda _: None)
            worker = depot.Depot(agent)
            with patch.object(depot.shutil, 'which', side_effect=lambda name: name == 'afplay'):
                try:
                    worker.start()
                    self.assertTrue(fifth.wait(8), 'player stopped before the fifth audio chunk')
                    self.assertGreaterEqual(calls[-1] - calls[0], 3.5)
                finally:
                    worker.close()
            self.assertFalse(worker.thread.is_alive())

"""Check the engineer's cadence and full-day variety independently of playback."""
import unittest
from unittest import mock
import joker


class JokerTests(unittest.TestCase):
    def test_full_day_without_repeats(self):
        schedule = joker.Schedule()
        self.assertGreaterEqual(len(set(schedule.lines)), 144)
        heard = []
        with mock.patch('joker.time.monotonic', return_value=100):
            schedule.start()
            self.assertGreaterEqual(schedule.deadline, 700)
            self.assertLessEqual(schedule.deadline, 1300)
            self.assertIsNone(schedule.due())
            for _ in schedule.lines:
                schedule.deadline = 0
                heard.append(schedule.due().removeprefix('Train engineer Joe to BitMagic. '))
                self.assertIsNone(schedule.due())
        self.assertEqual(len(heard), len(set(heard)))

    def test_new_phase_preserves_deck(self):
        schedule = joker.Schedule()
        schedule.deadline = 0
        schedule.due()
        remaining = list(schedule.deck)
        schedule.start()
        self.assertEqual(schedule.deck, remaining)

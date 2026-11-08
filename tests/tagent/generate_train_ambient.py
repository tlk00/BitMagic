#!/usr/bin/env python3
"""Mix CC0 train recordings into a quiet, looping mono ambience.

Inputs must be 16-bit mono WAV at 22050 Hz; source attribution and acquisition
are documented in README.md. No external packages are needed for mixing.
"""
import argparse
import array
import math
import sys
import wave
from pathlib import Path

RATE = 22050
OVERLAP_SECONDS = 12
TARGET_RMS = 0.014


def read_segment(path, start, duration):
    """Read and normalize a selected segment of a 22050 Hz, 16-bit mono WAV."""
    with wave.open(str(path), 'rb') as stream:
        if (stream.getnchannels(), stream.getsampwidth(), stream.getframerate()) != (1, 2, RATE):
            raise ValueError('input must be 22050 Hz, 16-bit mono WAV: ' + str(path))
        stream.setpos(start * RATE)
        samples = array.array('h', stream.readframes(duration * RATE))
    if sys.byteorder != 'little':
        samples.byteswap()
    if len(samples) != duration * RATE:
        raise ValueError('input is shorter than the selected segment: ' + str(path))
    rms = math.sqrt(sum(value * value for value in samples) / len(samples)) / 32768
    if rms == 0:
        raise ValueError('input segment is silent: ' + str(path))
    return samples, TARGET_RMS / rms


def render(first_path, second_path, third_path, output):
    """Crossfade three normalized inputs and add occasional BITMAGIC Morse."""
    first, first_gain = read_segment(first_path, 20, 200)
    second, second_gain = read_segment(second_path, 20, 36)
    third, third_gain = read_segment(third_path, 2, 52)
    third_again, third_again_gain = read_segment(third_path, 0, 52)
    first_length, other_length = 200 * RATE, 36 * RATE
    overlap = OVERLAP_SECONDS * RATE
    freight_overlap = 4 * RATE
    freight_segment = 52 * RATE
    freight_length = freight_segment * 2 - freight_overlap
    second_start = first_length - overlap
    third_start = second_start + other_length - overlap
    length = third_start + freight_length
    def freight_at(index):
        """Repeat the freight recording with a four-second internal crossfade."""
        second_freight_start = freight_segment - freight_overlap
        if index < second_freight_start:
            return third[index] * third_gain
        if index < freight_segment:
            fraction = (index - second_freight_start) / freight_overlap
            return (third[index] * third_gain * (1 - fraction) +
                    third_again[index - second_freight_start] * third_again_gain * fraction)
        return third_again[index - second_freight_start] * third_again_gain

    def mix_at(index):
        """Return one normalized sample with a linear overlap of the sources."""
        if index < second_start:
            return first[index] * first_gain
        if index < first_length:
            fraction = (index - second_start) / overlap
            return (first[index] * first_gain * (1 - fraction) +
                    second[index - second_start] * second_gain * fraction)
        if index < third_start:
            return second[index - second_start] * second_gain
        if index < second_start + other_length:
            fraction = (index - third_start) / overlap
            return (second[index - second_start] * second_gain * (1 - fraction) +
                    freight_at(index - third_start) * fraction)
        return freight_at(index - third_start)

    data = array.array('h')
    for index in range(overlap, length):
        value = mix_at(index)
        if index >= length - overlap:
            fraction = (index - (length - overlap)) / overlap
            value = value * (1 - fraction) + mix_at(index - (length - overlap)) * fraction
        data.append(max(-32768, min(32767, int(value))))

    morse = {'B': '-...', 'I': '..', 'T': '-', 'M': '--', 'A': '.-',
             'G': '--.', 'C': '-.-.'}
    for beginning in (20.0, 105.0, 190.0):
        cursor = beginning
        for letter in 'BITMAGIC':
            for symbol in morse[letter]:
                duration = .11 if symbol == '.' else .33
                start, frames = int(cursor * RATE), int(duration * RATE)
                for i in range(frames):
                    elapsed = i / RATE
                    envelope = min(1, elapsed / .012, (duration - elapsed) / .012)
                    value = data[start + i] + int(32767 * .012 * envelope *
                                                  math.sin(2 * math.pi * 660 * elapsed))
                    data[start + i] = max(-32768, min(32767, value))
                cursor += duration + .11
            cursor += .22
    if sys.byteorder != 'little':
        data.byteswap()
    with wave.open(str(output), 'wb') as stream:
        stream.setparams((1, 2, RATE, 0, 'NONE', 'not compressed'))
        stream.writeframes(data.tobytes())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('wheels', type=Path, help='first converted source WAV')
    parser.add_argument('between_wagons', type=Path, help='second converted source WAV')
    parser.add_argument('freight', type=Path, help='third converted source WAV')
    parser.add_argument('--output', type=Path,
                        default=Path(__file__).resolve().parent / 'assets' / 'train-ambient.wav')
    args = parser.parse_args()
    render(args.wheels, args.between_wagons, args.freight, args.output)

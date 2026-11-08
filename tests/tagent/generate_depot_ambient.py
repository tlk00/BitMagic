#!/usr/bin/env python3
"""Mix the CC BY-SA Sherman recording into mono depot ambience.

Input must be a 22050 Hz, 16-bit mono WAV converted from the source listed in
README.md. The generated asset is distributed under the same CC BY-SA license.
"""
import argparse
import array
import math
import random
import sys
import wave
from pathlib import Path

RATE = 22050
SECONDS = 60


def render(source_path, output_path):
    """Mix a continuous Sherman passage and soft workshop sounds without steam."""
    rng = random.Random(8273)
    samples = array.array('f', [0.0]) * (RATE * SECONDS)
    low = 0.0
    for i in range(len(samples)):
        t = i / RATE
        low = .994 * low + .006 * rng.uniform(-1, 1)
        fade = min(1, t / 3, (SECONDS - t) / 3)
        samples[i] = fade * (.009 * math.sin(2 * math.pi * 55 * t) +
                             .004 * math.sin(2 * math.pi * 83 * t) + .12 * low)

    # Short, damped metal taps from a distant workbench.
    for base, taps, frequency in ((5.8, (0, .46), 390), (15.1, (0, .37, .83), 480),
                                  (27.8, (0, .51), 320), (39.6, (0, .32, .75), 420),
                                  (52.3, (0, .55), 370)):
        for number, offset in enumerate(taps):
            start = int((base + offset) * RATE)
            for i in range(int(.25 * RATE)):
                t = i / RATE
                body = math.sin(2 * math.pi * frequency * t) + .25 * math.sin(
                    2 * math.pi * frequency * 1.43 * t)
                samples[start + i] += .025 * (1 if number == 0 else .72) * math.exp(-t * 30) * body

    with wave.open(str(source_path), 'rb') as source:
        if (source.getnchannels(), source.getsampwidth(), source.getframerate()) != (1, 2, RATE):
            raise ValueError('source must be 22050 Hz, 16-bit mono WAV')
        # This steady passage avoids the approach/departure envelopes in the
        # earlier excerpts. Blend its last five seconds into its beginning.
        source.setpos(670 * RATE)
        passage = array.array('h', source.readframes(65 * RATE))
        if sys.byteorder != 'little':
            passage.byteswap()
        if len(passage) != 65 * RATE:
            raise ValueError('source is shorter than the requested passage')
        rms = math.sqrt(sum(value * value for value in passage) / len(passage)) / 32768
        gain = .85 * min(1.0, .05 / rms) if rms else 0.0
        overlap = 5 * RATE
        for i in range(SECONDS * RATE):
            value = passage[i + overlap]
            if i >= (SECONDS * RATE - overlap):
                offset = i - (SECONDS * RATE - overlap)
                fraction = offset / overlap
                value = value * (1 - fraction) + passage[offset] * fraction
            samples[i] += gain * value / 32768
    pcm = array.array('h', (max(-32768, min(32767, int(value * 32767))) for value in samples))
    if sys.byteorder != 'little':
        pcm.byteswap()
    with wave.open(str(output_path), 'wb') as output:
        output.setparams((1, 2, RATE, 0, 'NONE', 'not compressed'))
        output.writeframes(pcm.tobytes())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('sherman', type=Path, help='converted Sherman source WAV')
    parser.add_argument('--output', type=Path,
                        default=Path(__file__).resolve().parent / 'assets' / 'depot.wav')
    args = parser.parse_args()
    render(args.sherman, args.output)

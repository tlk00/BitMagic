"""Steampunk/dieselpunk soundtrack with synthesized sounds and a bundled voice.

44.1 kHz mono PCM synthesis uses seeded noise, harmonics and soft envelopes. The
engine pulse rate represents worker occupancy, not measured CPU utilization.
"""
import array
import math
from pathlib import Path
import random
import sys
import wave

RATE = 44100
EVENTS = ('start', 'building', 'build', 'depart', 'work', 'success', 'fail', 'idle')
THEMES = ('steampunk', 'dieselpunk')


TELEGRAM = 'New BitMagic updates, new BitMagic updates'
MORSE = dict(zip('ABCDEFGHIJKLMNOPQRSTUVWXYZ', (
    '.-', '-...', '-.-.', '-..', '.', '..-.', '--.', '....', '..', '.---',
    '-.-', '.-..', '--', '-.', '---', '.--.', '--.-', '.-.', '...', '-',
    '..-', '...-', '.--', '-..-', '-.--', '--..')))
MORSE[','] = '--..--'


def render_dispatch(path, volume=0.35):
    """Write the bundled mono radio message at the requested volume.

    path is the destination WAV; volume must be finite and in [0,1]. The source
    is the approved voice/steam/wheel mix at reference volume 0.5. Returns clip
    duration in seconds. No audio device or runtime speech service is needed.
    """
    if not math.isfinite(volume) or not 0 <= volume <= 1:
        raise ValueError('Volume must be between 0 and 1')
    asset = Path(__file__).resolve().parent / 'assets' / 'dispatch.wav'
    with wave.open(str(asset), 'rb') as source:
        if (source.getnchannels(), source.getsampwidth(), source.getframerate()) != (1, 2, RATE):
            raise ValueError('Invalid dispatch audio format')
        frames = source.getnframes()
        mono = array.array('h', source.readframes(frames))
    if sys.byteorder != 'little':
        mono.byteswap()
    scaled = array.array('h')
    scale = volume / 0.5
    for sample in mono:
        scaled.append(max(-32768, min(32767, round(sample * scale))))
    if sys.byteorder != 'little':
        scaled.byteswap()
    with wave.open(str(path), 'wb') as destination:
        destination.setparams((1, 2, RATE, 0, 'NONE', 'not compressed'))
        destination.writeframes(scaled.tobytes())
    return frames / RATE


def morse_intervals(text, unit=0.105, lead=0.3):
    """Return key-down (start,end) seconds for text, using standard Morse spacing.

    unit is the positive dot duration; lead is initial silence in seconds.
    Unsupported symbols raise ValueError. Dashes last three units; element,
    letter and word gaps last one, three and seven units respectively.
    """
    if not math.isfinite(unit) or unit <= 0 or not math.isfinite(lead) or lead < 0:
        raise ValueError('Invalid Morse timing')
    intervals = []
    cursor = lead
    for word in text.upper().split():
        for letter in word:
            if letter not in MORSE:
                raise ValueError('Unsupported Morse character: ' + letter)
            for mark in MORSE[letter]:
                length = unit * (3 if mark == '-' else 1)
                intervals.append((cursor, cursor + length))
                cursor += length + unit
            cursor += 2 * unit
        cursor += 4 * unit
    return intervals


def telegraph_track(intervals, duration):
    """Return an independent mono float track of dry metal key impacts.

    intervals contains Morse key-down (start,end) seconds; duration bounds the
    allocated track. Only contact closure and release sound, for 50 ms each.
    Nothing sounds during a held dash. An independent random source generates
    contact transients, so the steam noise cannot become the keyed carrier.
    """
    track = array.array('f', [0]) * int(duration * RATE)
    key_rng = random.Random(8191)
    for start, end in intervals:
        for edge, accent, pitch in ((start, 0.95, 1.0), (end, 0.72, 1.23)):
            first = int(edge * RATE)
            for offset in range(int(0.05 * RATE)):
                index = first + offset
                if not 0 <= index < len(track):
                    continue
                elapsed = offset / RATE
                # A sharp click into a small steel stop, with very short
                # inharmonic resonances; no sustained tone or keyed noise bed.
                attack = 1 - math.exp(-elapsed * 14000)
                metal = (
                    0.46 * math.sin(2 * math.pi * 2016 * pitch * elapsed) * math.exp(-elapsed * 110)
                    + 0.26 * math.sin(2 * math.pi * 3685 * pitch * elapsed) * math.exp(-elapsed * 190)
                    + 0.11 * math.sin(2 * math.pi * 5352 * pitch * elapsed) * math.exp(-elapsed * 380)
                    + 0.04 * math.sin(2 * math.pi * 650 * elapsed) * math.exp(-elapsed * 450))
                contact = 0.25 * key_rng.uniform(-1, 1) * math.exp(-elapsed * 1600)
                track[index] += accent * attack * math.tanh(1.5 * (metal + contact))
    return track


def pulse_rate(workers, budget):
    """Return engine pulses/second for active workers and a positive job budget."""
    return 1.5 + 6.5 * min(1, max(0, workers) / max(1, budget))


def wheel_clatter(t, rate, noise):
    """Return one stylized rail-joint sample with two pairs of axle impacts.

    t is seconds, rate is carriage-pattern repetitions/second, and noise is a
    filtered random sample in [-1,1]. The quieter pickup and accented second
    strike make a ta-DAM rhythm. This cadence is independent of cylinder exhaust;
    it is an audible stylization rather than a rolling-stock physics model.
    """
    phase = (t * rate) % 1
    value = 0.0
    for offset, accent in ((0.0, 0.72), (0.13, 1.0), (0.48, 0.65), (0.61, 0.94)):
        elapsed = ((phase - offset) % 1) / rate
        # Very short broadband attack, then a damped body knock instead of a bell.
        attack = (1 - math.exp(-elapsed * 3500)) * math.exp(-elapsed * 145)
        body = (1 - math.exp(-elapsed * 1800)) * math.exp(-elapsed * 32)
        value += accent * (0.65 * noise * attack + body * (
            0.48 * math.sin(2 * math.pi * 82 * elapsed)
            + 0.18 * math.sin(2 * math.pi * 205 * elapsed)))
    return value


def render(path, theme, event, volume=0.35, workers=1, budget=1):
    """Write a bounded mono WAV to path; return duration in seconds.

    theme/event select a known sound. volume is amplitude in [0,1]. workers and
    budget set engine cadence. No audio device is accessed by this function.
    """
    if theme not in THEMES or event not in EVENTS:
        raise ValueError('Unknown sound theme or event')
    if not math.isfinite(volume) or not 0 <= volume <= 1:
        raise ValueError('Volume must be between 0 and 1')
    steam = theme == 'steampunk'
    duration = {'start': 1.1, 'building': 1.6, 'build': 1.4, 'work': 3.0,
                'depart': 3.0, 'success': 2.3, 'fail': 6.0, 'idle': 1.5}[event]
    if steam and event in ('work', 'depart'):
        duration = 23.0 if event == 'depart' else 10.5
    telegram = morse_intervals(TELEGRAM) if steam and event == 'building' else []
    if telegram:
        duration = telegram[-1][1] + 0.8
    key_track = telegraph_track(telegram, duration) if telegram else None
    rng = random.Random(1729)
    samples = array.array('h')
    low_noise = 0.0
    steam_noise = 0.0
    cadence = pulse_rate(workers, budget)
    if steam:
        cadence *= 0.7
    rail_rate = 0.50 + 0.35 * min(1, max(0, workers) / max(1, budget))
    # True Morse SOS: three dots, three dashes, three dots, repeated three times.
    sos = []
    cursor = 0.0
    for repetition in range(3):
        for group, length in enumerate((1, 3, 1)):
            for symbol in range(3):
                sos.append((cursor, cursor + length * 0.055))
                cursor += (length + 1) * 0.055
            cursor += 0.11
        cursor += 0.22
    for index in range(int(duration * RATE)):
        t = index / RATE
        noise = rng.uniform(-1, 1)
        low_noise = 0.94 * low_noise + 0.06 * noise
        steam_noise = 0.65 * steam_noise + 0.35 * noise
        fade = min(1, t / 0.025, (duration - t) / 0.08)
        value = 0.0
        if event in ('work', 'depart'):
            phase = (t * cadence) % 1
            if steam:
                # Paired cylinder exhaust: an accented CHUFF followed by a
                # softer chuff. Filtered noise gives steam body rather than clicks.
                # Outside observer: a gradual approach, locomotive first, then a
                # longer freight-car tail. Integrating a mild Doppler ratio keeps
                # the source clock continuous through the closest point.
                crossing = duration * 0.36
                transition = 0.45
                source_time = t - 0.06 * transition * (
                    math.log(math.cosh((t - crossing) / transition))
                    - math.log(math.cosh(-crossing / transition)))
                departing = event == 'depart'
                if departing:
                    # Start at the platform, then accelerate: integrate speed to
                    # obtain a continuous axle/exhaust clock, rather than jumping
                    # the frequency at successive animation frames.
                    # Double the scene length without lowering the pitch: the
                    # event/acceleration timeline runs at half speed, while tonal
                    # oscillators and noise remain at the original sample rate.
                    scene_time = t / 2
                    moving = max(0, scene_time - 1.6)
                    source_time = 2 * (0.27 * moving + 0.10 * moving * moving)
                    cadence = 2.5
                pair_phase = (source_time * cadence / 2) % 1
                chuff = 0.0
                for offset, accent in ((0.0, 1.0), (0.48, 0.82)):
                    local = pair_phase - offset
                    if local >= 0:
                        chuff += accent * (1 - math.exp(-local * 100)) * math.exp(-local * 10)
                engine_distance = (t - crossing) / (duration * 0.19)
                engine_gain = 0.12 + 1.15 / (1 + engine_distance * engine_distance)
                if departing:
                    engine_gain = 1.1
                value = engine_gain * chuff * (2.4 * steam_noise + 0.65 * noise)
                value += engine_gain * chuff * 0.35 * math.sin(2 * math.pi * 95 * source_time)
                # Freight bogies trail the locomotive: slower, low paired thuds,
                # with broad rolling rumble rather than a high metallic jingle.
                wagon_distance = (t - duration * 0.57) / (duration * 0.24)
                wagon_gain = 0.06 + 0.94 / (1 + wagon_distance * wagon_distance)
                if departing:
                    wagon_gain = min(1, max(0, (scene_time - 1.6) / 2))
                value += wagon_gain * wheel_clatter(source_time, rail_rate, steam_noise)
                value += wagon_gain * (0.3 * low_noise +
                    0.065 * math.sin(2 * math.pi * 37 * source_time) +
                    0.04 * math.sin(2 * math.pi * 54 * source_time))
                if departing:
                    # Steam valves hiss under the prolonged bell/whistle scene.
                    value += steam_noise * (0.10 + 0.18 * math.exp(
                        -((scene_time - 2.5) / 1.2) ** 2))
                    # Platform bell overlaps the steam whistle and first chuffs.
                    for strike_at in (0.1, 1.0, 1.9, 2.8):
                        elapsed = t - strike_at * 2
                        if elapsed >= 0:
                            envelope = (1 - math.exp(-elapsed * 120)) * math.exp(-elapsed * 2.5)
                            value += envelope * (
                                0.52 * math.sin(2 * math.pi * 698.5 * elapsed)
                                + 0.20 * math.sin(2 * math.pi * 1928 * elapsed))
                    whistle_time = scene_time - 0.8
                    if 0 < whistle_time < 2.6:
                        envelope = min(1, whistle_time / 0.25, (2.6 - whistle_time) / 0.4)
                        value += envelope * (0.06 * noise + sum(
                            0.19 * math.sin(2 * math.pi * hz * (whistle_time * 2) +
                                            0.5 * math.sin(2 * math.pi * 4.5 * t))
                            for hz in (330, 412.5, 495)))
            else:
                # Low crankshaft harmonics, cylinder knocks and filtered rumble.
                value = 0.26 * math.sin(2 * math.pi * 48 * t)
                value += 0.15 * math.sin(2 * math.pi * 96 * t)
                value += (0.5 * low_noise + 0.22 * noise) * math.exp(-phase * 15)
        elif event == 'building':
            if steam:
                # Slow steam crescendo only; Morse never gates its amplitude.
                progress = t / duration
                swell = progress * progress * (3 - 2 * progress)
                value = steam_noise * (0.04 + 0.5 * swell)
            else:
                value = 0.3 * math.sin(2 * math.pi * (30 * t + 24 * t * t)) + 0.2 * low_noise
        elif event == 'fail':
            if any(start <= t < end for start, end in sos):
                value = 0.48 * math.sin(2 * math.pi * (880 if steam else 660) * t)
        elif event in ('success', 'start'):
            # Steam whistle chord or a diesel air horn, with gentle vibrato.
            fundamental = 440 if steam else 146.83
            vibrato = 0.6 * math.sin(2 * math.pi * 5 * t)
            gate = 1 if event == 'success' else max(0, math.sin(math.pi * t / duration))
            value = gate * sum(0.22 * math.sin(2 * math.pi * fundamental * ratio * t + vibrato)
                               for ratio in ((1, 1.25, 1.5) if steam else (1, 1.5, 2)))
            value += 0.07 * noise if steam else 0.04 * low_noise
        elif event == 'build':
            # Four celebratory mechanical bell/ignition notes: cha-cha.
            note = int(t / 0.3)
            local = t % 0.3
            hz = (523.25, 659.25, 783.99, 1046.5, 783.99)[min(note, 4)]
            if not steam:
                hz *= 0.5
            value = math.exp(-local * 12) * (0.4 * math.sin(2 * math.pi * hz * t)
                    + 0.16 * math.sin(2 * math.pi * hz * 2.76 * t))
        elif event == 'idle':
            value = 0.32 * math.sin(2 * math.pi * (330 * t - 75 * t * t / duration))
            value *= math.sin(math.pi * t / duration)
            value += 0.13 * noise * (1 - t / duration) if steam else 0.06 * low_noise
        value = math.tanh(value) * fade * volume
        # Small stereo spread; mono playback remains compatible.
        if steam and event in ('work', 'depart'):
            # The locomotive crosses first; the spatial centre of the carriages
            # follows later. Equal-power panning also folds down cleanly to mono.
            pan = 0.85 * math.tanh((t - duration * 0.48) / (duration * 0.22))
            left = math.sqrt((1 - pan) / 2)
            right = math.sqrt((1 + pan) / 2)
            if event == 'depart':
                pan = 0.1 + 0.75 * (t / duration)
                left = math.sqrt((1 - pan) / 2)
                right = math.sqrt((1 + pan) / 2)
                value *= 1 - 0.85 * (t / duration) ** 1.3
            else:
                value *= math.sin(math.pi * t / duration) ** 0.45
        else:
            pan = 0.08 * math.sin(2 * math.pi * 0.7 * t)
            left = 1 - pan
            right = 1 + pan
        if key_track is not None:
            # Mix independent stems after shaping the steam, so a key impact
            # cannot compress or modulate the background. Slightly place the
            # key to the left while the steam remains centred.
            key = key_track[index] * fade * volume
            left_sample = value * 0.55 + key * 0.80
            right_sample = value * 0.55 + key * 0.55
        else:
            left_sample = value * left
            right_sample = value * right
        # Fold both positions into a single mono sample for all output devices.
        samples.append(int(32767 * (left_sample + right_sample) / 2))
    if sys.byteorder != 'little':
        samples.byteswap()
    with wave.open(str(path), 'wb') as stream:
        stream.setparams((1, 2, RATE, 0, 'NONE', 'not compressed'))
        stream.writeframes(samples.tobytes())
    return duration

#!/usr/bin/env python3
"""
12-stage Circle inverse FFT (interpolation) over Mersenne31.

Fixed configuration:
    log_n = 12
    N     = 4096

The transform follows the Plonky3-style Circle FFT structure:
  * M31 field p = 2^31 - 1
  * standard circle domain
  * 12 DIF interpolation stages
  * stage 0 uses y-derived twiddles
  * later stages use x-derived twiddles
  * all twiddles are batch-inverted before the inverse FFT

All 4096 input evaluation values are printed before the transform.
"""

from __future__ import annotations

import random
import time
from dataclasses import dataclass
from typing import List

P = (1 << 31) - 1
LOG_N = 12
N = 1 << LOG_N
FULL_GENERATOR = (311_014_874, 1_584_694_829)


def add_mod(a: int, b: int) -> int:
    s = a + b
    if s >= P:
        s -= P
    return s


def sub_mod(a: int, b: int) -> int:
    d = a - b
    if d < 0:
        d += P
    return d


def mul_mod(a: int, b: int) -> int:
    return (a * b) % P


@dataclass(frozen=True)
class Point:
    x: int
    y: int

    def add(self, other: "Point") -> "Point":
        return Point(
            (self.x * other.x - self.y * other.y) % P,
            (self.x * other.y + self.y * other.x) % P,
        )

    def double(self) -> "Point":
        return Point(
            (2 * self.x * self.x - 1) % P,
            (2 * self.x * self.y) % P,
        )


def generator(bits: int) -> Point:
    g = Point(*FULL_GENERATOR)
    for _ in range(31 - bits):
        g = g.double()
    return g


def bit_reverse(x: int, bits: int) -> int:
    y = 0
    for _ in range(bits):
        y = (y << 1) | (x & 1)
        x >>= 1
    return y


def bit_reverse_list(values: List[Point]) -> List[Point]:
    n = len(values)
    bits = n.bit_length() - 1
    return [values[bit_reverse(i, bits)] for i in range(n)]


def standard_coset0(log_n: int) -> List[Point]:
    half = 1 << (log_n - 1)
    shift = generator(log_n + 1)
    subgroup_generator = generator(log_n - 1)

    points = []
    cur = shift
    for _ in range(half):
        points.append(cur)
        cur = cur.add(subgroup_generator)

    return points


def compute_twiddles(log_n: int) -> List[List[int]]:
    points = bit_reverse_list(standard_coset0(log_n))

    twiddles = []
    twiddles.append([pt.y for pt in points])
    twiddles.append([points[i].x for i in range(0, len(points), 2)])

    for _ in range(log_n - 2):
        prev = twiddles[-1]
        cur = [
            (2 * prev[i] * prev[i] - 1) % P
            for i in range(0, len(prev), 2)
        ]
        twiddles.append(cur)

    return twiddles


def batch_inverse(xs: List[int]) -> List[int]:
    n = len(xs)
    prefix = [1] * n

    acc = 1
    for i, x in enumerate(xs):
        if x == 0:
            raise ZeroDivisionError("zero twiddle")
        prefix[i] = acc
        acc = (acc * x) % P

    acc_inv = pow(acc, P - 2, P)

    out = [0] * n
    for i in range(n - 1, -1, -1):
        out[i] = (acc_inv * prefix[i]) % P
        acc_inv = (acc_inv * xs[i]) % P

    return out


def compute_inverse_twiddles(log_n: int) -> List[List[int]]:
    return [batch_inverse(layer) for layer in compute_twiddles(log_n)]


def circle_ifft_in_place(values: List[int], inv_twiddles: List[List[int]]) -> List[float]:
    n = len(values)
    stage_times = []

    for inv_ts in inv_twiddles:
        t0 = time.perf_counter()

        num_blocks = len(inv_ts)
        block_size = n // num_blocks
        half = block_size >> 1

        for block, tw in enumerate(inv_ts):
            base = block * block_size
            upper = base + half

            for j in range(half):
                i0 = base + j
                i1 = upper + j

                a = values[i0]
                b = values[i1]

                values[i0] = add_mod(a, b)
                values[i1] = mul_mod(sub_mod(a, b), tw)

        stage_times.append(time.perf_counter() - t0)

    inv_n = pow(n, P - 2, P)
    for i in range(n):
        values[i] = mul_mod(values[i], inv_n)

    return stage_times


def main() -> None:
    print("=" * 72)
    print("12-stage Circle inverse FFT / interpolation")
    print(f"Field              : M31 = {P}")
    print(f"log_n              : {LOG_N}")
    print(f"N                  : {N}")
    print(f"Stages             : {LOG_N}")
    print(f"Butterflies/stage  : {N // 2}")
    print(f"Total butterflies  : {(N // 2) * LOG_N}")
    print("=" * 72)

    rng = random.Random(1)
    input_values = [rng.randrange(P) for _ in range(N)]

    print()
    print("ALL INPUT VALUES")
    print("-" * 72)
    for i, value in enumerate(input_values):
        print(f"input[{i:4d}] = {value}")

    print()
    print("Preparing inverse twiddles ...")
    t0 = time.perf_counter()
    inv_twiddles = compute_inverse_twiddles(LOG_N)
    t1 = time.perf_counter()

    print(f"Inverse twiddle generation : {(t1 - t0) * 1e3:.3f} ms")
    print(
        "Twiddles per stage        : "
        + ", ".join(str(len(x)) for x in inv_twiddles)
    )

    coefficients = input_values.copy()

    print()
    print("Running 12-stage Circle inverse FFT ...")
    t0 = time.perf_counter()
    stage_times = circle_ifft_in_place(coefficients, inv_twiddles)
    t1 = time.perf_counter()

    print(f"Total IFFT latency         : {(t1 - t0) * 1e3:.3f} ms")
    print(f"Total IFFT latency         : {(t1 - t0) * 1e6:.3f} us")

    print()
    print("PER-STAGE LATENCY")
    print("-" * 72)
    for stage, dt in enumerate(stage_times):
        print(
            f"stage {stage:2d}: "
            f"twiddles={len(inv_twiddles[stage]):4d}, "
            f"butterflies={N // 2:4d}, "
            f"time={dt * 1e3:9.3f} ms"
        )

    print()
    print("ALL OUTPUT COEFFICIENTS")
    print("-" * 72)
    for i, value in enumerate(coefficients):
        print(f"output[{i:4d}] = {value}")


if __name__ == "__main__":
    main()

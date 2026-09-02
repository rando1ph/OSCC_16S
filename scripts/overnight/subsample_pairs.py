#!/usr/bin/env python3
"""Deterministic, pairing-preserving subsampler for paired FASTQ.

Reservoir sampling with a fixed seed. Reads both mates in lockstep and emits
at most --n pairs. Also reports the total number of pairs read.

Usage:
  subsample_pairs.py --r1 IN_R1 --r2 IN_R2 --out1 OUT_R1 --out2 OUT_R2
                     --n MAX_PAIRS [--seed SEED]

Writes "<out1>.meta.tsv" with raw_pairs and retained_pairs.
"""
import argparse
import gzip
import random
import sys


def open_maybe_gz(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rb")
    return open(path, "rb")


def read_fastq_pairs(r1, r2):
    with open_maybe_gz(r1) as f1, open_maybe_gz(r2) as f2:
        while True:
            rec1 = []
            rec2 = []
            for _ in range(4):
                line1 = f1.readline()
                line2 = f2.readline()
                if not line1 or not line2:
                    if line1 or line2:
                        raise ValueError("unbalanced paired FASTQ")
                    return
                rec1.append(line1)
                rec2.append(line2)
            yield b"".join(rec1), b"".join(rec2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--out1", required=True)
    ap.add_argument("--out2", required=True)
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--seed", type=int, default=20260901)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    reservoir1 = []
    reservoir2 = []
    total = 0

    for rec1, rec2 in read_fastq_pairs(args.r1, args.r2):
        total += 1
        if total <= args.n:
            reservoir1.append(rec1)
            reservoir2.append(rec2)
        else:
            j = rng.randrange(total)
            if j < args.n:
                reservoir1[j] = rec1
                reservoir2[j] = rec2

    retained = len(reservoir1)
    if retained == 0:
        # write empty (truncated) gzip files
        for path in (args.out1, args.out2):
            with gzip.open(path, "wb"):
                pass
    else:
        with gzip.open(args.out1, "wb") as o1, gzip.open(args.out2, "wb") as o2:
            for r1, r2 in zip(reservoir1, reservoir2):
                o1.write(r1)
                o2.write(r2)

    with open(args.out1 + ".meta.tsv", "w") as m:
        m.write("raw_pairs\tretained_pairs\n")
        m.write(f"{total}\t{retained}\n")

    print(f"raw_pairs={total} retained_pairs={retained} cap={args.n}")
    if total != retained:
        print("subsampled", flush=True)


if __name__ == "__main__":
    main()

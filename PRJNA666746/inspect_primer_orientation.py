import gzip
from collections import Counter

R1 = "raw/SRR12750669_1.fastq.gz"
R2 = "raw/SRR12750669_2.fastq.gz"

FORWARD = "CCTACGGGNGGCWGCAG"
REVERSE = "GACTACHVGGGTATCTAATCC"

IUPAC = {
    "A": set("A"),
    "C": set("C"),
    "G": set("G"),
    "T": set("T"),
    "N": set("ACGT"),
    "W": set("AT"),
    "H": set("ACT"),
    "V": set("ACG"),
}

def mismatches(seq, primer):
    seq = seq.upper()
    if len(seq) < len(primer):
        return 999

    n = 0
    for base, code in zip(seq, primer):
        if base not in IUPAC[code]:
            n += 1
    return n

def classify(seq):
    mf = mismatches(seq, FORWARD)
    mr = mismatches(seq, REVERSE)

    # 允许最多 2 个错配，只用于方向识别
    if mf <= 2 and mf < mr:
        return "F"
    if mr <= 2 and mr < mf:
        return "R"
    return "U"

r1_counts = Counter()
r2_counts = Counter()
pair_counts = Counter()

total = 0

with gzip.open(R1, "rt") as f1, gzip.open(R2, "rt") as f2:
    while True:
        h1 = f1.readline()
        h2 = f2.readline()

        if not h1 and not h2:
            break

        if not h1 or not h2:
            raise RuntimeError("R1/R2 have different read counts")

        s1 = f1.readline().strip()
        s2 = f2.readline().strip()

        f1.readline()
        f2.readline()
        f1.readline()
        f2.readline()

        c1 = classify(s1)
        c2 = classify(s2)

        r1_counts[c1] += 1
        r2_counts[c2] += 1
        pair_counts[c1 + c2] += 1
        total += 1

print("Total read pairs:", total)

print("\n=== R1 primer orientation ===")
for x in ["F", "R", "U"]:
    print(x, r1_counts[x], f"{r1_counts[x]/total*100:.2f}%")

print("\n=== R2 primer orientation ===")
for x in ["F", "R", "U"]:
    print(x, r2_counts[x], f"{r2_counts[x]/total*100:.2f}%")

print("\n=== Pair orientation ===")
for x in ["FR", "RF", "FF", "RR", "FU", "RU", "UF", "UR", "UU"]:
    print(x, pair_counts[x], f"{pair_counts[x]/total*100:.2f}%")

proper = pair_counts["FR"] + pair_counts["RF"]

print("\nOpposite-primer pairs:", proper)
print("Opposite-primer percentage:", f"{proper/total*100:.2f}%")

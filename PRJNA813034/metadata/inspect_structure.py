import csv
import re
from collections import Counter, defaultdict

src = "ena_run_metadata.tsv"

counts = Counter()
stage_group = Counter()

tumor_keys = {}
adjacent_keys = {}
precancer = []

with open(src, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for r in reader:
        title = r["experiment_title"]
        run = r["run_accession"]
        lib = r["library_name"]

        if "Precancer_" in title:
            m = re.search(r'Precancer_(\d+)', title)
            idx = int(m.group(1)) if m else None
            counts["Precancer"] += 1
            precancer.append((idx, run, lib, title))
            continue

        low = title.lower()

        m_stage = re.search(r't([1-4])_stage', low)
        m_idx = re.search(r'_(\d+)\s*$', title)

        if not m_stage or not m_idx:
            counts["Unparsed"] += 1
            print("UNPARSED:", run, title)
            continue

        stage = f"T{m_stage.group(1)}"
        idx = int(m_idx.group(1))
        key = (stage, idx)

        if "tissue adjacent to a tumor" in low:
            group = "Adjacent"
            counts[group] += 1
            stage_group[(stage, group)] += 1
            adjacent_keys[key] = (run, lib, title)

        elif "tumor tissue" in low:
            group = "Tumor"
            counts[group] += 1
            stage_group[(stage, group)] += 1
            tumor_keys[key] = (run, lib, title)

        else:
            counts["Other"] += 1
            print("OTHER:", run, title)

print("=== GROUP COUNTS ===")
for k in ["Tumor", "Adjacent", "Precancer", "Unparsed", "Other"]:
    print(f"{k}: {counts[k]}")

print("\n=== STAGE x GROUP ===")
for stage in ["T1", "T2", "T3", "T4"]:
    print(
        stage,
        "Tumor:", stage_group[(stage, "Tumor")],
        "Adjacent:", stage_group[(stage, "Adjacent")]
    )

matched = sorted(set(tumor_keys) & set(adjacent_keys))
adj_only = sorted(set(adjacent_keys) - set(tumor_keys))
tumor_only = sorted(set(tumor_keys) - set(adjacent_keys))

print("\n=== PAIR CHECK ===")
print("Matched Tumor-Adjacent keys:", len(matched))
print("Adjacent without matching Tumor:", len(adj_only))
print("Tumor without matching Adjacent:", len(tumor_only))

print("\n=== MATCHED KEYS ===")
for key in matched:
    print(f"{key[0]}_{key[1]:02d}")

print("\n=== ADJACENT ONLY ===")
for key in adj_only:
    print(key)

print("\n=== TUMOR ONLY COUNT BY STAGE ===")
c = Counter(k[0] for k in tumor_only)
for stage in ["T1", "T2", "T3", "T4"]:
    print(stage, c[stage])

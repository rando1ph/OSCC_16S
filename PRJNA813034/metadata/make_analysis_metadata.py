import csv
import re
from collections import Counter, defaultdict

src = "ena_run_metadata.tsv"
out = "analysis_metadata.tsv"

records = []

with open(src, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for r in reader:
        title = r["experiment_title"]
        low = title.lower()

        if "precancer_" in low:
            continue

        m_stage = re.search(r't([1-4])_stage', low)
        m_idx = re.search(r'_(\d+)\s*$', title)

        if not m_stage or not m_idx:
            continue

        stage = f"T{m_stage.group(1)}"
        idx = int(m_idx.group(1))

        # 主分析只保留有严格配对的 T1/T4 样本 1-10
        if stage not in {"T1", "T4"} or idx > 10:
            continue

        if "tissue adjacent to a tumor" in low:
            group = "Matched_Normal"
        elif "tumor tissue" in low:
            group = "Tumor"
        else:
            continue

        patient_id = f"{stage}_P{idx:02d}"

        records.append({
            "sample-id": r["run_accession"],
            "patient_id": patient_id,
            "group": group,
            "stage": stage,
            "sample_index": idx,
            "library_name": r["library_name"],
            "experiment_title": title
        })

records.sort(key=lambda x: (x["patient_id"], x["group"]))

with open(out, "w", newline="") as f:
    fields = [
        "sample-id",
        "patient_id",
        "group",
        "stage",
        "sample_index",
        "library_name",
        "experiment_title"
    ]
    writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(records)

groups = Counter(r["group"] for r in records)
patients = defaultdict(list)

for r in records:
    patients[r["patient_id"]].append(r["group"])

bad = {
    p: g for p, g in patients.items()
    if sorted(g) != ["Matched_Normal", "Tumor"]
}

print("Samples:", len(records))
print("Tumor:", groups["Tumor"])
print("Matched_Normal:", groups["Matched_Normal"])
print("Unique patients:", len(patients))
print("Incomplete/abnormal pairs:", len(bad))

print("Stage counts:")
for stage in ["T1", "T4"]:
    t = sum(r["stage"] == stage and r["group"] == "Tumor" for r in records)
    n = sum(r["stage"] == stage and r["group"] == "Matched_Normal" for r in records)
    print(stage, "Tumor:", t, "Matched_Normal:", n)

print("Written:", out)

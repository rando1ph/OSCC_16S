import csv
import re
from collections import Counter, defaultdict

src = "biosample_metadata.tsv"
out = "analysis_metadata.tsv"

records = []

with open(src, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for r in reader:
        isolate = r["isolate"]

        m = re.search(r'N(\d+)_', isolate)
        if not m:
            raise ValueError(f"Cannot extract patient ID: {isolate}")

        patient_id = f"P{int(m.group(1)):02d}"

        if "_OSCC" in isolate:
            group = "Tumor"
        elif "Matched Control" in isolate:
            group = "Matched_Normal"
        else:
            raise ValueError(f"Unknown group: {isolate}")

        records.append({
            "sample-id": r["run_accession"],
            "patient_id": patient_id,
            "group": group,
            "sample_alias": r["sample_alias"],
            "title": r["title"],
            "isolate": isolate,
            "tissue": r["tissue"],
            "age": r["age"],
            "sex": r["sex"],
        })

with open(out, "w", newline="") as f:
    fields = [
        "sample-id",
        "patient_id",
        "group",
        "sample_alias",
        "title",
        "isolate",
        "tissue",
        "age",
        "sex"
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

if bad:
    for p, g in bad.items():
        print("BAD:", p, g)

print("Written:", out)

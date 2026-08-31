import csv
import re
from collections import Counter, defaultdict

src = "ena_experiment_metadata.tsv"
out = "analysis_metadata.tsv"

records = []

with open(src, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for r in reader:
        alias = r["sample_alias"]
        library = r["library_name"]

        m = re.search(r'_S(\d+)$', alias)
        if not m:
            raise ValueError(f"Cannot extract patient ID from: {alias}")

        patient_id = f"P{int(m.group(1)):03d}"

        if library.startswith("Tumor_"):
            group = "Tumor"
        elif library.startswith("AN_"):
            group = "Matched_Normal"
        else:
            raise ValueError(f"Unknown library group: {library}")

        records.append({
            "sample-id": r["run_accession"],
            "patient_id": patient_id,
            "group": group,
            "sample_accession": r["sample_accession"],
            "sample_alias": alias,
            "experiment_accession": r["experiment_accession"],
            "library_name": library
        })

records.sort(key=lambda x: (x["patient_id"], x["group"]))

with open(out, "w", newline="") as f:
    fields = [
        "sample-id",
        "patient_id",
        "group",
        "sample_accession",
        "sample_alias",
        "experiment_accession",
        "library_name"
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

import csv
import os
from collections import defaultdict

base = os.path.expanduser("~/OSCC_16S/PRJNA813034")
raw = os.path.join(base, "raw")

analysis_runs = set()

with open("analysis_metadata.tsv", newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        analysis_runs.add(row["sample-id"])

files = defaultdict(dict)

with open("ena_run_metadata.tsv", newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        run = row["run_accession"]

        if run not in analysis_runs:
            continue

        urls = row["fastq_ftp"].split(";")

        if len(urls) != 2:
            raise ValueError(
                f"{run}: expected 2 FASTQ files, found {len(urls)}"
            )

        for url in urls:
            filename = url.rsplit("/", 1)[-1]

            if filename.endswith("_1.fastq.gz"):
                files[run]["forward"] = os.path.join(raw, filename)

            elif filename.endswith("_2.fastq.gz"):
                files[run]["reverse"] = os.path.join(raw, filename)

            else:
                raise ValueError(f"{run}: unexpected filename {filename}")

missing = []

for run in sorted(analysis_runs):
    if set(files[run]) != {"forward", "reverse"}:
        missing.append(run)

if missing:
    raise ValueError(f"Missing paired FASTQ mapping: {missing}")

with open("manifest.tsv", "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow([
        "sample-id",
        "forward-absolute-filepath",
        "reverse-absolute-filepath"
    ])

    for run in sorted(analysis_runs):
        writer.writerow([
            run,
            files[run]["forward"],
            files[run]["reverse"]
        ])

print("Analysis samples:", len(analysis_runs))
print("Manifest rows:", len(files))
print("Missing pairs:", len(missing))
print("Written: manifest.tsv")

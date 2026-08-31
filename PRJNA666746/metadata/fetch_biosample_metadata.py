import csv
import time
import urllib.request
import xml.etree.ElementTree as ET

input_file = "ena_run_samples.tsv"
output_file = "biosample_metadata.tsv"

rows = []
with open(input_file, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    rows = list(reader)

cache = {}

with open(output_file, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow([
        "run_accession",
        "sample_accession",
        "sample_alias",
        "title",
        "isolate",
        "tissue",
        "age",
        "sex"
    ])

    for i, row in enumerate(rows, 1):
        sample = row["sample_accession"]

        if sample not in cache:
            url = f"https://www.ebi.ac.uk/ena/browser/api/xml/{sample}"
            with urllib.request.urlopen(url, timeout=30) as r:
                xml_data = r.read()

            root = ET.fromstring(xml_data)

            title_node = root.find(".//TITLE")
            title = title_node.text.strip() if title_node is not None and title_node.text else ""

            attrs = {}
            for attr in root.findall(".//SAMPLE_ATTRIBUTE"):
                tag_node = attr.find("TAG")
                val_node = attr.find("VALUE")
                if tag_node is not None and val_node is not None:
                    tag = tag_node.text.strip() if tag_node.text else ""
                    val = val_node.text.strip() if val_node.text else ""
                    attrs[tag.lower()] = val

            cache[sample] = {
                "title": title,
                "isolate": attrs.get("isolate", ""),
                "tissue": attrs.get("tissue", ""),
                "age": attrs.get("age", ""),
                "sex": attrs.get("sex", "")
            }

            time.sleep(0.1)

        m = cache[sample]

        writer.writerow([
            row["run_accession"],
            sample,
            row["sample_alias"],
            m["title"],
            m["isolate"],
            m["tissue"],
            m["age"],
            m["sex"]
        ])

        print(f"[{i}/{len(rows)}] {sample}  {m['title']}  {m['isolate']}")

print(f"\nDone: {output_file}")

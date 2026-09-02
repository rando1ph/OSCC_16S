#!/usr/bin/env python3
"""Taxonomy filtering + genus collapse for one cohort.

Usage: taxonomy_processing.py <cohort_base_dir> <table_tsv> <taxonomy_tsv>
Outputs:
  <results>/genus_count.tsv       (samples x genera counts, genus-indexed)
  <results>/genus_relabund.tsv    (samples x genera relative abundance)
  <results>/taxonomy_filtered.tsv (feature -> taxonomy kept)
"""
import sys
import os
import re

import numpy as np
import pandas as pd


def genus_of(tax_string):
    if not isinstance(tax_string, str) or tax_string.strip() == "":
        return "Unclassified"
    parts = [p.strip() for p in tax_string.split(";") if p.strip()]
    for p in reversed(parts):
        if p.startswith("g__") and len(p) > 3:
            return p
    return "Unclassified"


def main():
    base, table_tsv, tax_tsv = sys.argv[1], sys.argv[2], sys.argv[3]
    results = os.path.join(base, "results")

    table = pd.read_csv(table_tsv, sep="\t", index_col=0)  # samples x features (QIIME2 view)
    tax = pd.read_csv(tax_tsv, sep="\t", index_col=0)
    tax = tax["Taxon"] if "Taxon" in tax.columns else tax.iloc[:, 0]

    if len(table.index.intersection(tax.index)) == 0 and len(table.columns.intersection(tax.index)) > 0:
        table = table.T  # make features x samples
    table = table.loc[table.index.intersection(tax.index)]
    tax = tax.loc[table.index]

    # keep only Bacteria domain
    keep = tax.apply(lambda s: str(s).startswith(("d__Bacteria", "k__Bacteria")))
    table = table.loc[keep]
    tax = tax.loc[keep]

    # remove mitochondria / chloroplast / inappropriate organelles
    bad = tax.apply(
        lambda s: bool(re.search(r"(Mitochondria|Chloroplast|Cyanobacteria/Chloroplast)", str(s), re.IGNORECASE))
    )
    table = table.loc[~bad]
    tax = tax.loc[~bad]

    # collapse to genus
    genus = tax.apply(genus_of)
    genus.name = "genus"
    table["__genus__"] = genus
    genus_count = table.groupby("__genus__").sum()
    genus_count.index.name = "genus"
    genus_count = genus_count.astype(int)

    genus_relabund = genus_count.div(genus_count.sum(axis=0), axis=1).replace(np.nan, 0.0)

    genus_count.to_csv(os.path.join(results, "genus_count.tsv"), sep="\t")
    genus_relabund.to_csv(os.path.join(results, "genus_relabund.tsv"), sep="\t", float_format="%.8g")

    kept_tax = pd.DataFrame({"Taxon": tax, "genus": genus})
    kept_tax.to_csv(os.path.join(results, "taxonomy_filtered.tsv"), sep="\t")

    print(f"features_kept={len(tax)} genera={len(genus_count)}")
    print(f"wrote genus_count.tsv and genus_relabund.tsv")


if __name__ == "__main__":
    main()

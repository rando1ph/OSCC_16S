#!/usr/bin/env python3
"""Helpers for converting QIIME2 artifacts to plain tables.

Requires the QIIME2-2026.7 conda python.
"""
import sys
import os

os.environ.setdefault(
    "R_HOME",
    os.path.expanduser("~/miniconda3/envs/rachis-qiime2-2026.7/lib/R"))

import qiime2


def table_to_tsv(qza, out_tsv):
    import pandas as pd
    t = qiime2.Artifact.load(qza)
    df = t.view(pd.DataFrame)
    df = df.fillna(0).astype(int)
    df.to_csv(out_tsv, sep="\t")
    print(f"table {qza} -> {out_tsv} shape={df.shape}")


def main():
    cmd = sys.argv[1]
    if cmd == "table-to-tsv":
        table_to_tsv(sys.argv[2], sys.argv[3])
    elif cmd == "taxa-to-tsv":
        import pandas as pd
        tax = qiime2.Artifact.load(sys.argv[2])
        df = tax.view(pd.DataFrame)
        df.to_csv(sys.argv[3], sep="\t")
        print(f"taxonomy -> {sys.argv[3]} rows={len(df)}")
    else:
        print("unknown command", cmd)
        sys.exit(1)


if __name__ == "__main__":
    main()

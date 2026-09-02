#!/usr/bin/env python3
"""Paired Tumor vs Matched_Normal genus-level statistics for one cohort.

Usage: genus_stats.py <cohort_base_dir>
Outputs (in <base>/results/):
  paired_genus_stats.tsv  per-genus paired statistics
  clr_matrix.tsv          sample x genus CLR (pseudocount applied)
  sample_status.tsv       per-sample inclusion / exclusion reason
"""
import sys
import os

import numpy as np
import pandas as pd
from scipy import stats

PSEUDOCOUNT = 0.5
LOW_DEPTH_THRESHOLD = 2000


def bh_fdr(pvals):
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    if n == 0:
        return np.array([])
    order = np.argsort(pvals)
    ranked = pvals[order]
    ranked_q = ranked * n / np.arange(1, n + 1)
    ranked_q = np.minimum.accumulate(ranked_q[::-1])[::-1]
    ranked_q = np.minimum(ranked_q, 1.0)
    qvals = np.empty(n)
    qvals[order] = ranked_q
    return qvals


def main():
    base = sys.argv[1]
    res = os.path.join(base, "results")
    meta = pd.read_csv(os.path.join(base, "metadata", "analysis_metadata.tsv"), sep="\t")
    meta = meta[["sample-id", "patient_id", "group"]]
    meta.columns = ["sample", "patient", "group"]

    count = pd.read_csv(os.path.join(res, "genus_count.tsv"), sep="\t", index_col=0)
    relabund = pd.read_csv(os.path.join(res, "genus_relabund.tsv"), sep="\t", index_col=0)

    # depth status / presence
    depth_path = os.path.join(res, "depth_status.tsv")
    depth = pd.read_csv(depth_path, sep="\t", index_col=0) if os.path.exists(depth_path) else None
    if depth is not None and "depth_status" not in depth.columns:
        depth["depth_status"] = np.where(
            depth["non_chimeric_reads"] >= LOW_DEPTH_THRESHOLD, "OK", "LOW_DEPTH")

    samples = count.columns
    status_rows = []
    reasons = {}
    for s in samples:
        if depth is not None and s in depth.index:
            dstat = depth.loc[s, "depth_status"]
            nreads = depth.loc[s, "non_chimeric_reads"]
        else:
            dstat = "absent"
            nreads = 0
        if dstat == "LOW_DEPTH":
            reasons[s] = f"LOW_DEPTH ({nreads} < {LOW_DEPTH_THRESHOLD})"
        elif dstat != "OK":
            reasons[s] = f"not_OK ({dstat})"
        else:
            reasons[s] = None
        status_rows.append([s, nreads, "included" if reasons[s] is None else "excluded"])

    # apply paired exclusion: if either member of a pair is excluded, exclude both
    pairs = {}
    for _, row in meta.iterrows():
        pairs.setdefault(row["patient"], {}).setdefault(row["group"], row["sample"])

    for patient, grp in pairs.items():
        t = grp.get("Tumor")
        n = grp.get("Matched_Normal")
        if t is None or n is None:
            continue
        if reasons.get(t) is not None or reasons.get(n) is not None:
            for s in (t, n):
                if reasons.get(s) is None:
                    reasons[s] = "pair_member_excluded"
                    for r in status_rows:
                        if r[0] == s:
                            r[2] = "excluded"

    # complete pairs among samples with genus data
    complete_pairs = []
    for patient, grp in pairs.items():
        t = grp.get("Tumor")
        n = grp.get("Matched_Normal")
        if t is None or n is None:
            continue
        if t not in samples or n not in samples:
            continue
        if reasons.get(t) is not None or reasons.get(n) is not None:
            continue
        complete_pairs.append((patient, t, n))

    # CLR transform
    X = count.T.copy()  # samples x genera
    Xc = X + PSEUDOCOUNT
    clr = np.log(Xc) - np.log(Xc).mean(axis=1).values[:, None]
    clr = pd.DataFrame(clr, index=X.index, columns=X.columns)

    # per-genus paired statistics
    rows = []
    for genus in count.index:
        effect = []
        for _patient, t, n in complete_pairs:
            effect.append(clr.loc[t, genus] - clr.loc[n, genus])
        effect = np.array(effect, dtype=float)
        n_pairs = len(effect)
        if n_pairs == 0:
            p = np.nan
            q = np.nan
            med_effect = np.nan
        else:
            med_effect = float(np.median(effect))
            if n_pairs >= 2 and np.std(effect) > 0:
                p = float(stats.wilcoxon(effect).pvalue)
            else:
                p = float(np.nan)

        med_tumor = float(relabund.loc[genus, [t for _, t, _ in complete_pairs]].median())
        med_normal = float(relabund.loc[genus, [n for _, _, n in complete_pairs]].median())

        rows.append([genus, n_pairs, med_tumor, med_normal, med_effect, p, np.nan])

    df = pd.DataFrame(rows, columns=[
        "genus", "complete_pairs", "median_tumor_relabund",
        "median_normal_relabund", "median_paired_clr_effect", "wilcoxon_p", "bh_fdr_q"])
    df["bh_fdr_q"] = bh_fdr(df["wilcoxon_p"].fillna(1.0))

    df.to_csv(os.path.join(res, "paired_genus_stats.tsv"), sep="\t", float_format="%.6g",
              index=False)

    clr.to_csv(os.path.join(res, "clr_matrix.tsv"), sep="\t", float_format="%.6g")

    status = pd.DataFrame(status_rows, columns=["sample", "non_chimeric_reads", "status"])
    status["reason"] = status["sample"].map(reasons).fillna("")
    status.to_csv(os.path.join(res, "sample_status.tsv"), sep="\t", index=False)

    ex = []
    for patient, grp in pairs.items():
        t = grp.get("Tumor")
        n = grp.get("Matched_Normal")
        for s in (t, n):
            if s is not None and reasons.get(s) is not None:
                ex.append([patient, s, reasons[s]])
    pd.DataFrame(ex, columns=["patient", "sample", "reason"]).to_csv(
        os.path.join(res, "pair_exclusions.tsv"), sep="\t", index=False)

    print(f"genera={len(df)} complete_pairs={len(complete_pairs)} "
          f"excluded_samples={len(ex)}")


if __name__ == "__main__":
    main()

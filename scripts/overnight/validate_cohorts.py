#!/usr/bin/env python3
"""Validate three cohort-level genus outputs against their underlying tables.

Recomputes genus relative abundance from genus counts and the per-genus paired
statistics from genus_count + analysis_metadata + depth_status, then compares
against the on-disk validated result files.  Reports PASS/FAIL per check.

Usage: validate_cohorts.py [cohort1 cohort2 ...]
Exit code 0 only if every check passes.
"""
import sys
import os

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
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


def recompute_stats(base):
    res = os.path.join(base, "results")
    meta = pd.read_csv(os.path.join(base, "metadata", "analysis_metadata.tsv"), sep="\t")
    meta = meta[["sample-id", "patient_id", "group"]]
    meta.columns = ["sample", "patient", "group"]

    count = pd.read_csv(os.path.join(res, "genus_count.tsv"), sep="\t", index_col=0)
    relabund = pd.read_csv(os.path.join(res, "genus_relabund.tsv"), sep="\t", index_col=0)
    assert list(count.index) == list(relabund.index), "genus index mismatch"

    depth = pd.read_csv(os.path.join(res, "depth_status.tsv"), sep="\t", index_col=0)
    depth.columns = [c if c != depth.columns[0] else "sample-id" for c in depth.columns]
    depth["depth_status"] = np.where(
        depth["non_chimeric_reads"] >= LOW_DEPTH_THRESHOLD, "OK", "LOW_DEPTH")

    samples = count.columns
    reasons = {}
    for s in samples:
        if s in depth.index:
            dstat = depth.loc[s, "depth_status"]
            nreads = depth.loc[s, "non_chimeric_reads"]
        else:
            dstat = "absent"
            nreads = 0
        if dstat == "LOW_DEPTH":
            reasons[s] = f"LOW_DEPTH ({nreads})"
        elif dstat != "OK":
            reasons[s] = f"not_OK ({dstat})"

    pairs = {}
    for _, row in meta.iterrows():
        pairs.setdefault(row["patient"], {}).setdefault(row["group"], row["sample"])

    for patient, grp in pairs.items():
        t, n = grp.get("Tumor"), grp.get("Matched_Normal")
        if t is None or n is None:
            continue
        if reasons.get(t) is not None or reasons.get(n) is not None:
            for s in (t, n):
                reasons.setdefault(s, "pair_member_excluded")

    complete_pairs = []
    for patient, grp in pairs.items():
        t, n = grp.get("Tumor"), grp.get("Matched_Normal")
        if t is None or n is None:
            continue
        if t not in samples or n not in samples:
            continue
        if reasons.get(t) is not None or reasons.get(n) is not None:
            continue
        complete_pairs.append((patient, t, n))

    X = count.T.copy()
    Xc = X + PSEUDOCOUNT
    clr = np.log(Xc) - np.log(Xc).mean(axis=1).values[:, None]
    clr = pd.DataFrame(clr, index=X.index, columns=X.columns)

    rows = []
    for genus in count.index:
        effect = [clr.loc[t, genus] - clr.loc[n, genus] for _, t, n in complete_pairs]
        effect = np.array(effect, dtype=float)
        n_pairs = len(effect)
        if n_pairs == 0:
            p, q, med_effect = np.nan, np.nan, np.nan
        else:
            med_effect = float(np.median(effect))
            if n_pairs >= 2 and np.std(effect) > 0:
                p = float(stats.wilcoxon(effect).pvalue)
            else:
                p = np.nan
        med_tumor = float(relabund.loc[genus, [t for _, t, _ in complete_pairs]].median())
        med_normal = float(relabund.loc[genus, [n for _, _, n in complete_pairs]].median())
        rows.append([genus, n_pairs, med_tumor, med_normal, med_effect, p, np.nan])

    df = pd.DataFrame(rows, columns=[
        "genus", "complete_pairs", "median_tumor_relabund",
        "median_normal_relabund", "median_paired_clr_effect", "wilcoxon_p", "bh_fdr_q"])
    df["bh_fdr_q"] = bh_fdr(df["wilcoxon_p"].fillna(1.0))
    return df, len(complete_pairs)


def validate(base):
    cohort = os.path.basename(base)
    res = os.path.join(base, "results")
    print(f"\n=== {cohort} ===")
    status = True

    def check(name, cond, detail=""):
        nonlocal status
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
        if not cond:
            status = False

    # file presence and non-empty
    files = ["paired_genus_stats.tsv", "genus_count.tsv", "genus_relabund.tsv",
             "depth_status.tsv", "sample_status.tsv"]
    for f in files:
        p = os.path.join(res, f)
        check(f"{f} exists & non-empty", os.path.exists(p) and os.path.getsize(p) > 0)

    stat = pd.read_csv(os.path.join(res, "paired_genus_stats.tsv"), sep="\t")
    count = pd.read_csv(os.path.join(res, "genus_count.tsv"), sep="\t", index_col=0)
    relabund = pd.read_csv(os.path.join(res, "genus_relabund.tsv"), sep="\t", index_col=0)
    meta = pd.read_csv(os.path.join(base, "metadata", "analysis_metadata.tsv"), sep="\t")

    exp_cols = ["genus", "complete_pairs", "median_tumor_relabund",
                "median_normal_relabund", "median_paired_clr_effect", "wilcoxon_p", "bh_fdr_q"]
    check("stats columns", list(stat.columns) == exp_cols)

    # recompute
    rec, n_complete = recompute_stats(base)
    check("complete_pairs max == recomputed", int(stat["complete_pairs"].max()) == n_complete,
          f"on-disk max={int(stat['complete_pairs'].max())} recomputed={n_complete}")
    check("stats genus set == count genus set",
          set(stat["genus"]) == set(count.index),
          f"stats n={len(stat)} count n={len(count)}")

    # compare recomputed vs on-disk values
    # NOTE: on-disk files were written with %.6g so ~6 significant figures.
    merged = stat.merge(rec, on="genus", suffixes=("_disk", "_rec"))
    rt = 1e-4
    check("median_tumor_relabund matches", np.allclose(
        merged["median_tumor_relabund_disk"], merged["median_tumor_relabund_rec"], rtol=rt, atol=1e-5))
    check("median_normal_relabund matches", np.allclose(
        merged["median_normal_relabund_disk"], merged["median_normal_relabund_rec"], rtol=rt, atol=1e-5))
    check("median_paired_clr_effect matches", np.allclose(
        merged["median_paired_clr_effect_disk"], merged["median_paired_clr_effect_rec"], rtol=rt, atol=1e-4))
    check("wilcoxon_p matches", np.allclose(
        merged["wilcoxon_p_disk"].fillna(-1), merged["wilcoxon_p_rec"].fillna(-1), rtol=rt, atol=1e-6))
    check("bh_fdr_q matches", np.allclose(
        merged["bh_fdr_q_disk"], merged["bh_fdr_q_rec"], rtol=rt, atol=1e-6))

    # relabund internal consistency
    cs = count.sum(axis=0)
    rs = relabund.sum(axis=0)
    check("genus_relabund columns sum to 1", bool(np.allclose(rs, 1.0, atol=1e-5)))
    check("genus index count == relabund", list(count.index) == list(relabund.index))
    check("count columns == relabund columns", list(count.columns) == list(relabund.columns))

    # sample counts vs metadata
    check("samples in count == metadata samples",
          set(count.columns) == set(meta["sample-id"]),
          f"count n={count.shape[1]} meta n={len(meta)}")

    # sample_status sanity
    ss = pd.read_csv(os.path.join(res, "sample_status.tsv"), sep="\t")
    check("sample_status non-empty & covers all samples",
          len(ss) == len(meta) and set(ss["sample"]) == set(meta["sample-id"]))
    if len(ss):
        n_inc = (ss["status"] == "included").sum()
        n_exc = (ss["status"] == "excluded").sum()
        check("sample_status included+excluded == total",
              n_inc + n_exc == len(ss), f"included={n_inc} excluded={n_exc}")

    # expected known values per cohort
    known = {"PRJNA666746": (480, 50), "PRJNA822685": (239, 56), "PRJNA813034": (283, 20)}
    if cohort in known:
        n_gen, n_pairs = known[cohort]
        check(f"expected {n_gen} genera", len(stat) == n_gen)
        check(f"expected {n_pairs} complete pairs", int(stat["complete_pairs"].max()) == n_pairs)

    return status


def main():
    args = sys.argv[1:] or ["PRJNA666746", "PRJNA822685", "PRJNA813034"]
    all_ok = True
    for c in args:
        ok = validate(os.path.join(ROOT, c))
        all_ok = all_ok and ok
    print("\nOVERALL:", "PASS" if all_ok else "FAIL")
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()

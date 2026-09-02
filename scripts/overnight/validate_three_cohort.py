#!/usr/bin/env python3
"""Final validation of three-cohort integration outputs in combined/three_cohort/.

Checks every deliverable for non-empty content, expected row counts, sensible
sample/pair/genera figures, coherent directions against the per-cohort stats
files, and consistency across tables.  Exit code 0 only if all checks pass.
"""
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TD = os.path.join(ROOT, "combined", "three_cohort")
COHORTS = ["PRJNA666746", "PRJNA822685", "PRJNA813034"]
EXPECT = {"PRJNA666746": dict(pairs=50, genera=480, sig=389),
          "PRJNA822685": dict(pairs=56, genera=239, sig=168),
          "PRJNA813034": dict(pairs=20, genera=283, sig=0)}

FILES = ["pipeline_status.tsv", "cohort_summary.tsv", "cohort_effects.tsv",
         "genus_harmonization.tsv", "harmonization_unmatched.tsv",
         "cross_cohort_consistency.tsv", "consistent_genera.tsv",
         "candidate_ranking.tsv", "preliminary_report.md",
         "effect_heatmap.png", "top_candidates.png",
         "direction_dotplot.png", "pipeline_retention.png"]


def main():
    ok = True

    def check(name, cond, detail=""):
        nonlocal ok
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
        if not cond:
            ok = False

    print("=== combined/three_cohort deliverable validation ===")
    for f in FILES:
        p = os.path.join(TD, f)
        check(f"{f} exists and non-empty", os.path.exists(p) and os.path.getsize(p) > 0)

    summ = pd.read_csv(os.path.join(TD, "cohort_summary.tsv"), sep="\t")
    check("cohort_summary has 3 cohort rows", len(summ) == 3)
    for c in COHORTS:
        r = summ[summ["cohort"] == c]
        exp = EXPECT[c]
        if len(r):
            r = r.iloc[0]
            check(f"{c} complete_pairs == {exp['pairs']}", int(r["complete_pairs"]) == exp["pairs"],
                  f"got {int(r['complete_pairs'])}")
            check(f"{c} n_genera == {exp['genera']}", int(r["n_genera"]) == exp["genera"])
            check(f"{c} n_genera_sig_q0.1 == {exp['sig']}", int(r["n_genera_sig_q0.1"]) == exp["sig"])

    eff = pd.read_csv(os.path.join(TD, "cohort_effects.tsv"), sep="\t")
    exp_rows = sum(EXPECT[c]["genera"] for c in COHORTS)
    check(f"cohort_effects row count == {exp_rows}", len(eff) == exp_rows,
          f"got {len(eff)}")
    check("cohort_effects covers 3 cohorts", set(eff["cohort"]) == set(COHORTS))

    # spot-check effect estimates against the source per-cohort stats files
    for c in COHORTS:
        src = pd.read_csv(os.path.join(ROOT, c, "results", "paired_genus_stats.tsv"),
                          sep="\t").set_index("genus")
        sub = eff[eff["cohort"] == c].set_index("genus")
        common = sub.index.intersection(src.index)
        okc = np.allclose(sub.loc[common, "median_paired_clr_effect"],
                          src.loc[common, "median_paired_clr_effect"], rtol=1e-6, atol=1e-12) and \
            np.allclose(sub.loc[common, "bh_fdr_q"], src.loc[common, "bh_fdr_q"],
                        rtol=1e-6, atol=1e-12)
        check(f"cohort_effects {c} matches source stats (n={len(common)})", bool(okc))

    harm = pd.read_csv(os.path.join(TD, "genus_harmonization.tsv"), sep="\t")
    check("genus_harmonization non-empty rows", len(harm) >= 500,
          f"{len(harm)} labels in union catalogue")

    unm = pd.read_csv(os.path.join(TD, "harmonization_unmatched.tsv"), sep="\t")
    n_un = len(unm)
    n_shared = len(harm[harm["status"] == "shared"])
    check("harmonization shared + unmatched == total",
          n_shared + n_un == len(harm), f"shared={n_shared} unmatched={n_un}")

    cons = pd.read_csv(os.path.join(TD, "cross_cohort_consistency.tsv"), sep="\t")
    check("consistency table non-empty (>=2 cohorts)", len(cons) >= 100, f"{len(cons)} shared")
    check("consistency n_cohorts in {2,3}", set(cons["n_cohorts"]).issubset({2, 3}))
    check("consistency has all-3 rows", (cons["n_cohorts"] == 3).sum() >= 90)
    check("consistency direction_consistent counts",
          cons["direction_consistent"].isin([0, 1]).all())
    check("consistency prevalence within [0,1]",
          cons["mean_prevalence"].between(0, 1).all())

    # coherence of cohort_effects string direction vs effect pivot
    for _, r in cons.head(10).iterrows():
        g = r["genus"]
        for c in COHORTS:
            row = eff[(eff["genus"] == g) & (eff["cohort"] == c)]
            token = [t for t in str(r["cohort_effects"]).split(";")
                     if t.startswith(c.replace("PRJNA", ""))]
            if len(row) and len(token):
                v = float(row["median_paired_clr_effect"].iloc[0])
                sign_txt = "+" if v > 0 else "-"
                check(f"effect-string sign for {g} {c}", sign_txt in token[0])

    csf = pd.read_csv(os.path.join(TD, "consistent_genera.tsv"), sep="\t")
    check("consistent_genera non-empty", len(csf) >= 100, f"{len(csf)} rows")
    check("consistent_genera all direction_consistent==1",
          (csf["direction_consistent"] == 1).all())
    check("consistent_genera all have q<0.1 in >=1 cohort", (csf["n_sig_q0.1"] >= 1).all())

    cand = pd.read_csv(os.path.join(TD, "candidate_ranking.tsv"), sep="\t")
    check("candidate_ranking rows == consistency rows", len(cand) == len(cons))
    check("rank is sequential 1..N", list(cand["rank"]) == list(range(1, len(cand) + 1)))
    check("score finite & >0", np.isfinite(cand["score"]).all() and (cand["score"] > 0).all())
    check("tier present", cand["tier"].notna().all())

    # candidate table directions agree with consistency table
    m = cand.merge(cons, on="genus", suffixes=("_cand", "_cons"))
    check("ranking direction == consistency direction",
          (m["direction_cand"] == m["direction_cons"]).all())

    report = open(os.path.join(TD, "preliminary_report.md")).read()
    for kw in ["Greengenes 13_8", "DADA2", "sensitivity", "candidate",
               "harmonization", "not causal"]:
        check(f"report mentions '{kw}'", kw in report)

    # non-empty / sanity counts for report metrics paragraphs
    print("\n=== quick contents ===")
    print("shared (>=2 cohorts):", len(cons),
          "| all-3:", int((cons["n_cohorts"] == 3).sum()),
          "| consistent:", int(cons["direction_consistent"].sum()),
          "| consistent & sig:", len(csf))
    print("pipeline_status rows:", len(pd.read_csv(os.path.join(TD, "pipeline_status.tsv"),
                                                   sep="\t")))

    print("\nOVERALL:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

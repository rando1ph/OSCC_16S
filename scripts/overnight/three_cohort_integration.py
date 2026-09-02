#!/usr/bin/env python3
"""Three-cohort genus-level integration, harmonization, ranking and reporting.

Reads the validated per-cohort results:
  <COHORT>/results/paired_genus_stats.tsv
  <COHORT>/results/genus_count.tsv
  <COHORT>/results/sample_status.tsv
  <COHORT>/metadata/analysis_metadata.tsv
and produces, under <ROOT>/combined/three_cohort/:

  pipeline_status.tsv          three-cohort pipeline state
  cohort_summary.tsv           per-cohort summary (incl. taxonomy route)
  cohort_effects.tsv           long-format genus x cohort effect estimates
  genus_harmonization.tsv      union catalogue of genus labels across cohorts
  harmonization_unmatched.tsv  database-specific / unmatched genus labels
  cross_cohort_consistency.tsv shared-genus consistency table
  consistent_genera.tsv        shared genera, consistent direction, q<0.1 in >=1
  candidate_ranking.tsv        defensible candidate ranking
  effect_heatmap.png           heatmap of median paired CLR effects
  top_candidates.png           bar plot of top ranked genera
  direction_dotplot.png        per-cohort effect dot plot for top genera
  pipeline_retention.png       mean reads/pairs per stage per cohort
  preliminary_report.md        updated methods + results + caveats

PRJNA666746 and PRJNA822685 are DADA2_ASV cohorts with SILVA 138.2 taxonomy;
PRJNA813034 is an OTU_VSEARCH_97 (Greengenes 13_8) sensitivity cohort.  Raw
counts are never pooled across cohorts.  Per-cohort effect estimates are kept.
"""
import os
import sys
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
COHORTS = ["PRJNA666746", "PRJNA822685", "PRJNA813034"]
OUT = os.path.join(ROOT, "combined", "three_cohort")
PSEUDOCOUNT = 0.5
Q_THRESH = 0.1
P_THRESH = 0.05

PIPELINE = {
    "PRJNA666746": dict(route="DADA2_ASV", taxonomy="SILVA_138.2",
                        design="primary", database="SILVA"),
    "PRJNA822685": dict(route="DADA2_ASV", taxonomy="SILVA_138.2",
                        design="primary", database="SILVA"),
    "PRJNA813034": dict(route="OTU_VSEARCH_97", taxonomy="Greengenes_13_8",
                        design="sensitivity", database="Greengenes"),
}

KNOWN = {"PRJNA666746": (480, 50), "PRJNA822685": (239, 56),
         "PRJNA813034": (283, 20)}


def genus_key(label):
    """Return a comparison key for a genus label (strip database prefix)."""
    if label in ("Unclassified",):
        return "Unclassified"
    return re.sub(r"^[dgk]__", "", label) if label.startswith(
        ("g__", "d__", "k__")) else label


def read_stats(c):
    p = os.path.join(ROOT, c, "results", "paired_genus_stats.tsv")
    return pd.read_csv(p, sep="\t")


def read_counts(c):
    p = os.path.join(ROOT, c, "results", "genus_count.tsv")
    return pd.read_csv(p, sep="\t", index_col=0)


def read_status(c):
    p = os.path.join(ROOT, c, "results", "sample_status.tsv")
    return pd.read_csv(p, sep="\t")


def read_meta(c):
    p = os.path.join(ROOT, c, "metadata", "analysis_metadata.tsv")
    m = pd.read_csv(p, sep="\t")
    return m[["sample-id", "patient_id", "group"]]


def load_retention(c):
    p = os.path.join(ROOT, c, "processed", "retention.tsv")
    if os.path.exists(p):
        return pd.read_csv(p, sep="\t")
    return None


def db_specificity(label):
    """Heuristic note for database-specific / informative genus labels."""
    s = genus_key(label)
    if s == "Unclassified":
        return "pseudo-label (taxonomy-unclassified reads)"
    if re.search(r"_group$|_UCG-|_CAG-|wastewater-sludge|Incertae_Sedis|"
                 r"Family_XIII|_R-7_group|_NK[0-9A-Z]+_group", s):
        return "SILVA-style lineage/informative label (no literal Greengenes counterpart)"
    if re.match(r"^[0-9]", s) or re.match(r"^[A-Z][0-9]{2,3}$", s):
        return "Greengenes-style environmental/candidate label (no literal SILVA counterpart)"
    if s.startswith("["):
        return "SILVA bracketed label"
    return "plain genus name"


def prevalence_by_group(c):
    """Detection prevalence of each genus among included tumor / normal samples."""
    cnt = read_counts(c)
    meta = read_meta(c)
    status = read_status(c)
    inc = set(status.loc[status["status"] == "included", "sample"])
    grp = dict(zip(meta["sample-id"], meta["group"]))
    tum = [s for s in cnt.columns if s in inc and grp.get(s) == "Tumor"]
    nor = [s for s in cnt.columns if s in inc and grp.get(s) == "Matched_Normal"]
    pos = cnt[list(cnt.columns)].gt(0)
    prev_tumor = pos[tum].mean(axis=1) if tum else pd.Series(0.0, index=cnt.index)
    prev_normal = pos[nor].mean(axis=1) if nor else pd.Series(0.0, index=cnt.index)
    prev_all = pos[list(cnt.columns)].mean(axis=1)
    out = pd.DataFrame({"prev_tumor": prev_tumor, "prev_normal": prev_normal,
                        "prev_all": prev_all})
    return out


def cohort_summary(c):
    stat = read_stats(c)
    cnt = read_counts(c)
    status = read_status(c)
    ret = load_retention(c)
    route = PIPELINE[c]
    n_pairs = int(stat["complete_pairs"].max())
    n_sig = int((stat["bh_fdr_q"] < Q_THRESH).sum())
    n_nom = int((stat["wilcoxon_p"] < P_THRESH).sum())
    n_exc = int((status["status"] == "excluded").sum())
    n_low = int(((status["status"] == "excluded") & (
        status["reason"].astype(str).str.contains("LOW_DEPTH"))).sum()) if len(status) else 0
    row = {"cohort": c, "analysis_route": route["route"],
           "taxonomy": route["taxonomy"], "design": route["design"],
           "database": route["database"],
           "samples_processed": len(status),
           "complete_pairs": n_pairs, "excluded_samples": n_exc,
           "low_depth_excluded": n_low, "n_genera": len(cnt),
           "n_genera_nominal_p0.05": n_nom, "n_genera_sig_q0.1": n_sig}
    if ret is not None and len(ret):
        row["mean_raw_pairs"] = round(float(ret["raw_pairs"].mean()), 1)
        row["mean_final_retained_pairs"] = round(float(ret["final_retained_pairs"].mean()), 1)
    # bacterial-mapped reads after QC
    row["mean_bacterial_reads_postqc"] = round(float(cnt.sum(axis=0).mean()), 1)
    # cohorts that mark very low bacterial content
    if c == "PRJNA813034":
        row["low_bacterial_reads_samples"] = int((cnt.sum(axis=0) < 2000).sum())
    return row


def main():
    os.makedirs(OUT, exist_ok=True)

    summaries = []
    effects = []
    prev_by_cohort = {}
    for c in COHORTS:
        summaries.append(cohort_summary(c))
        s = read_stats(c).copy()
        s["cohort"] = c
        s["analysis_route"] = PIPELINE[c]["route"]
        s["taxonomy"] = PIPELINE[c]["taxonomy"]
        s["design"] = PIPELINE[c]["design"]
        effects.append(s[["cohort", "analysis_route", "taxonomy", "design", "genus",
                          "complete_pairs", "median_tumor_relabund",
                          "median_normal_relabund", "median_paired_clr_effect",
                          "wilcoxon_p", "bh_fdr_q"]])
        prev_by_cohort[c] = prevalence_by_group(c)

    summ = pd.DataFrame(summaries)
    summ.to_csv(os.path.join(OUT, "cohort_summary.tsv"), sep="\t", index=False,
                float_format="%.6g")

    eff = pd.concat(effects, ignore_index=True)
    eff.to_csv(os.path.join(OUT, "cohort_effects.tsv"), sep="\t", index=False,
               float_format="%.6g")

    # ---- genus harmonization (union catalogue) -----------------------------
    labels_by = {c: set(read_counts(c).index) for c in COHORTS}
    union = sorted(set().union(*labels_by.values()), key=lambda x: x.lower())
    silva_both = labels_by["PRJNA666746"] & labels_by["PRJNA822685"]
    gg = labels_by["PRJNA813034"]
    harm_rows = []
    for lab in union:
        in666 = lab in labels_by["PRJNA666746"]
        in822 = lab in labels_by["PRJNA822685"]
        ingg = lab in labels_by["PRJNA813034"]
        n = int(in666) + int(in822) + int(ingg)
        if n >= 2:
            status = "shared"
        else:
            if in666 or in822:
                status = "SILVA-only" if not ingg else "shared"
            else:
                status = "Greengenes-only"
        note = db_specificity(lab)
        harm_rows.append({
            "genus": lab,
            "genus_name": genus_key(lab),
            "PRJNA666746_SILVA": int(in666),
            "PRJNA822685_SILVA": int(in822),
            "PRJNA813034_Greengenes": int(ingg),
            "n_cohorts": n,
            "in_both_SILVA": int(in666 and in822),
            "status": status,
            "note": note,
        })
    harm = pd.DataFrame(harm_rows)
    harm.to_csv(os.path.join(OUT, "genus_harmonization.tsv"), sep="\t", index=False,
                float_format="%.6g")

    # ---- unmatched / database-specific documentation -----------------------
    un_rows = []
    for _, r in harm.iterrows():
        if r["status"] == "shared":
            continue
        where = []
        if r["PRJNA666746_SILVA"]:
            where.append("PRJNA666746")
        if r["PRJNA822685_SILVA"]:
            where.append("PRJNA822685")
        if r["PRJNA813034_Greengenes"]:
            where.append("PRJNA813034")
        un_rows.append({"genus": r["genus"], "genus_name": r["genus_name"],
                        "status": r["status"], "detected_in": ";".join(where),
                        "note": r["note"]})
    pd.DataFrame(un_rows).to_csv(os.path.join(OUT, "harmonization_unmatched.tsv"),
                                 sep="\t", index=False)

    # ---- shared genera metrics (>=2 cohorts with effect) -------------------
    pivot_eff = eff.pivot(index="genus", columns="cohort",
                          values="median_paired_clr_effect")
    pivot_p = eff.pivot(index="genus", columns="cohort", values="wilcoxon_p")
    pivot_q = eff.pivot(index="genus", columns="cohort", values="bh_fdr_q")
    pivot_prev = pd.DataFrame({c: prev_by_cohort[c]["prev_all"]
                               for c in COHORTS}).reindex(pivot_eff.index)
    pivot_medt = eff.pivot(index="genus", columns="cohort",
                           values="median_tumor_relabund")

    consistency_rows = []
    for genus in union:
        if genus == "Unclassified":
            continue
        present = [c for c in COHORTS if genus in labels_by[c]]
        if len(present) < 2:
            continue
        effs = pivot_eff.loc[genus]
        valid = [c for c in present if not np.isnan(effs[c])]
        if len(valid) < 2:
            continue
        signs = np.array([np.sign(effs[c]) for c in valid])
        signs = signs[signs != 0]
        if len(signs) < 2:
            # cannot assess direction agreement
            continue
        n_same = int(max((signs == 1).sum(), (signs == -1).sum()))
        frac_agree = n_same / len(signs)
        direction = "up" if (signs == 1).sum() >= (signs == -1).sum() else "down"
        n_up = int((signs == 1).sum())
        n_down = int((signs == -1).sum())
        all_consistent = int(n_up == len(signs) or n_down == len(signs))
        qs = [pivot_q.loc[genus, c] for c in valid if not np.isnan(pivot_q.loc[genus, c])]
        ps = [pivot_p.loc[genus, c] for c in valid if not np.isnan(pivot_p.loc[genus, c])]
        n_sig = int(np.sum(np.array(qs) < Q_THRESH)) if qs else 0
        n_nom = int(np.sum(np.array(ps) < P_THRESH)) if ps else 0
        # sign-matched, significance-weighted evidence across cohorts
        evidence = 0.0
        for c in valid:
            e = float(effs[c])
            match = (e > 0 and direction == "up") or (e < 0 and direction == "down")
            if not match:
                continue
            q = pivot_q.loc[genus, c]
            p = pivot_p.loc[genus, c]
            if not np.isnan(q) and q < Q_THRESH:
                evidence += 1.0
            elif not np.isnan(p) and p < P_THRESH:
                evidence += 0.5
        min_q = float(np.min(qs)) if qs else np.nan
        abs_eff = [abs(float(effs[c])) for c in valid]
        mean_abs = float(np.mean(abs_eff))
        mean_eff = float(np.mean([float(effs[c]) for c in valid]))
        prev_vals = [pivot_prev.loc[genus, c] for c in valid
                     if not np.isnan(pivot_prev.loc[genus, c])]
        mean_prev = float(np.mean(prev_vals)) if prev_vals else 0.0
        in_666 = genus in labels_by["PRJNA666746"]
        in_822 = genus in labels_by["PRJNA822685"]
        in_gg = genus in labels_by["PRJNA813034"]
        sil_count = int(in_666) + int(in_822)
        # per-cohort strings
        per = []
        for c in COHORTS:
            if c in valid:
                m = "+" if float(effs[c]) > 0 else "-"
                mark = "SILVA" if PIPELINE[c]["database"] == "SILVA" else "GG"
                tag = "SILVA" if c in ("PRJNA666746", "PRJNA822685") else "GG"
                q = pivot_q.loc[genus, c]
                star = "*" if not np.isnan(q) and q < Q_THRESH else ""
                per.append(f"{c.replace('PRJNA','')}({tag}{star}{m}{abs(float(effs[c])):.3g})")
            else:
                per.append(f"{c.replace('PRJNA','')}(NA)")
        consistency_rows.append({
            "genus": genus,
            "n_cohorts": len(valid),
            "n_DADA2_SILVA": sil_count,
            "in_sensitivity_GG": int(in_gg),
            "n_up": n_up, "n_down": n_down,
            "direction": direction,
            "direction_consistent": all_consistent,
            "frac_agree": frac_agree,
            "n_sig_q0.1": n_sig, "n_nominal_p0.05": n_nom,
            "evidence_hits": round(evidence, 1),
            "min_q": min_q,
            "mean_abs_effect": mean_abs,
            "mean_median_effect": mean_eff,
            "mean_prevalence": round(mean_prev, 4),
            "cohort_effects": ";".join(per),
        })

    cons = pd.DataFrame(consistency_rows)
    # order: consistency, sig, cohort count, prevalence, magnitude
    cons = cons.sort_values(
        ["direction_consistent", "n_sig_q0.1", "n_cohorts", "mean_prevalence",
         "mean_abs_effect"],
        ascending=[False, False, False, False, False]).reset_index(drop=True)
    cons.to_csv(os.path.join(OUT, "cross_cohort_consistency.tsv"), sep="\t",
                index=False, float_format="%.6g")

    # ---- consistent genera (same direction in >=2, q<0.1 in >=1) -----------
    cons_sig = cons[(cons["direction_consistent"] == 1) & (cons["n_sig_q0.1"] >= 1)]
    cons_sig.to_csv(os.path.join(OUT, "consistent_genera.tsv"), sep="\t",
                    index=False, float_format="%.6g")

    # ---- candidate ranking -------------------------------------------------
    # Transparent composite score (documented in the report):
    #   score = frac_agree
    #           * (1 + evidence_hits)              # sign-matched significance
    #           * (1 + 0.75*min(mean_abs_effect,3))# effect magnitude
    #           * (0.2 + 0.8*mean_prevalence)      # usable abundance
    #           * (1 + 0.35*(n_cohorts - 2))       # cohort coverage
    # discordant direction strongly penalised.
    score = []
    for _, r in cons.iterrows():
        consistency_f = float(r["frac_agree"])
        evidence = float(r["evidence_hits"])
        magnitude = 1.0 + 0.75 * min(float(r["mean_abs_effect"]), 3.0)
        abundance = 0.2 + 0.8 * float(r["mean_prevalence"])
        coverage = 1.0 + 0.35 * (int(r["n_cohorts"]) - 2)
        s = consistency_f * (1.0 + evidence) * magnitude * abundance * coverage
        if r["direction_consistent"] != 1:
            s *= 0.4
        score.append(s)
    cons["score"] = score

    # tier assignment
    tier = []
    for _, r in cons.iterrows():
        gg = int(r["in_sensitivity_GG"])
        sil = int(r["n_DADA2_SILVA"])
        n = int(r["n_cohorts"])
        dc = int(r["direction_consistent"])
        if sil == 2 and gg and n == 3 and dc:
            tier.append("A_all3_consistent")
        elif sil == 2 and gg and n == 3:
            tier.append("B_all3_discordant")
        elif sil == 2 and not gg and dc:
            tier.append("C_DADA2_consistent_GGnotdetected")
        elif sil == 2 and not gg:
            tier.append("D_DADA2_discordant_GGnotdetected")
        elif sil == 1 and gg and dc:
            tier.append("E_SILVA+GG_consistent")
        elif dc and n == 2:
            tier.append("F_twoCohort_consistent")
        else:
            tier.append("Z_discordant_or_partial")
    cons["tier"] = tier
    cons = cons.sort_values(
        ["direction_consistent", "score"],
        ascending=[False, False]).reset_index(drop=True)
    cons["rank"] = np.arange(1, len(cons) + 1)
    cons.to_csv(os.path.join(OUT, "candidate_ranking.tsv"), sep="\t",
                index=False, float_format="%.6g")

    # ---- pipeline status ---------------------------------------------------
    status_rows = []
    now = pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S")
    for c in COHORTS:
        status_rows.append({"cohort": c, "stage": "genus_stats",
                            "status": "VALIDATED", "updated": now,
                            "note": f"{PIPELINE[c]['route']} / {PIPELINE[c]['taxonomy']}"})
    status_rows.append({"cohort": "all", "stage": "three_cohort_integration",
                        "status": "COMPLETED", "updated": now,
                        "note": "genus-level; counts not pooled"})
    status_rows.append({"cohort": "all", "stage": "ml_leave_one_cohort_out",
                        "status": "NOT_RUN", "updated": now,
                        "note": "optional; deferred"})
    pd.DataFrame(status_rows).to_csv(os.path.join(OUT, "pipeline_status.tsv"),
                                     sep="\t", index=False)

    # ---- plots -------------------------------------------------------------
    make_plots(eff, pivot_eff, cons, summ)
    write_report(summ, cons, cons_sig, harm, un_rows)

    # write updated report to canonical location as well
    with open(os.path.join(OUT, "preliminary_report.md")) as f:
        rep = f.read()
    with open(os.path.join(ROOT, "combined", "preliminary_report.md"), "w") as f:
        f.write(rep)

    print("three-cohort outputs written to", OUT)


def make_plots(eff, pivot_eff, cons, summ):
    top = cons.head(25)
    order = top["genus"].tolist()
    # ---- heatmap (top ranked x cohorts) ----
    hmap = pivot_eff.loc[pivot_eff.index.intersection(order)].T
    hmap = hmap[[g for g in order if g in hmap.columns]]
    if hmap.size:
        vmax = np.nanmax(np.abs(hmap.values[np.isfinite(hmap.values)]))
        fig, ax = plt.subplots(figsize=(max(6, 0.42 * len(hmap.columns) + 3),
                                        0.5 * len(hmap.index) + 2.5))
        im = ax.imshow(hmap.values, cmap="RdBu_r", aspect="auto", vmin=-vmax, vmax=vmax)
        ax.set_xticks(range(len(hmap.columns)))
        ax.set_xticklabels(hmap.columns, rotation=45, ha="right", fontsize=7)
        ax.set_yticks(range(len(hmap.index)))
        ax.set_yticklabels([y.replace("PRJNA", "") for y in hmap.index], fontsize=8)
        ax.set_title("Median paired CLR effect (Tumor - Matched_Normal)\ntop 25 ranked genera")
        plt.colorbar(im, ax=ax, label="median paired CLR effect")
        plt.tight_layout()
        plt.savefig(os.path.join(OUT, "effect_heatmap.png"), dpi=150)
        plt.close(fig)

    # ---- top candidates bar ----
    top20 = cons.head(20)
    fig, ax = plt.subplots(figsize=(9, 8))
    colors = ["#d62728" if d > 0 else "#1f77b4"
              for d in top20["mean_median_effect"]]
    ax.barh(top20["genus"][::-1], top20["mean_median_effect"][::-1], color=colors[::-1])
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("mean median paired CLR effect (Tumor - Matched_Normal)")
    ax.set_title("Top 20 ranked candidate genera (three-cohort)")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "top_candidates.png"), dpi=150)
    plt.close(fig)

    # ---- direction dot plot ----
    dd = cons.head(15)
    fig, ax = plt.subplots(figsize=(10, 8))
    y_pos = np.arange(len(dd))[::-1]
    for i, (_, r) in enumerate(dd.iterrows()):
        yi = y_pos[i]
        for c in COHORTS:
            col = c.replace("PRJNA", "")
            if r["genus"] in pivot_eff.index and c in pivot_eff.columns:
                v = pivot_eff.loc[r["genus"], c]
                if not np.isnan(v):
                    db = PIPELINE[c]["database"]
                    marker = "o" if db == "SILVA" else "D"
                    ax.scatter(v, yi, marker=marker, s=70,
                               color="#d62728" if v > 0 else "#1f77b4",
                               edgecolor="k", linewidth=0.4,
                               label=(c + (" (GG)" if db == "Greengenes" else "")) if i == 0 else None)
            else:
                ax.scatter(np.nan, yi, marker="x", color="grey")
    ax.axvline(0, color="k", lw=0.8)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(dd["genus"])
    ax.set_xlabel("median paired CLR effect (Tumor - Matched_Normal)")
    ax.set_title("Per-cohort effect directions for top 15 ranked genera\n"
                 "(circle = SILVA DADA2 cohort, diamond = Greengenes OTU sensitivity)")
    ax.legend(loc="lower right", fontsize=7)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "direction_dotplot.png"), dpi=150)
    plt.close(fig)

    # ---- pipeline retention ----
    stages = ["raw", "pretrim\n(cap)", "valid primer", "final retained\n(cap)",
              "bacteria reads\npost-QC"]
    fig, ax = plt.subplots(figsize=(9, 6))
    for i, c in enumerate(COHORTS):
        ret = load_retention(c)
        if ret is None:
            continue
        cnt = read_counts(c)
        raw = ret["raw_pairs"].mean() * 2
        pre = ret["pretrim_pairs"].mean() * 2
        val = ret["valid_primer_pairs"].mean() * 2
        fin = ret["final_retained_pairs"].mean() * 2
        qc = cnt.sum(axis=0).mean()
        means = [raw, pre, val, fin, qc]
        ax.plot(range(len(stages)), means, "o-", lw=2,
                color=plt.cm.Set1(i), label=c)
        for x, y in zip(range(len(stages)), means):
            ax.text(x, y, f"{int(y):,}", fontsize=7, ha="center", va="bottom")
    ax.set_xticks(range(len(stages)))
    ax.set_xticklabels(stages)
    ax.set_yscale("log")
    ax.set_ylabel("mean reads/sample (log scale)")
    ax.set_title("Mean reads retained across pipeline stages (three cohorts)")
    ax.legend(fontsize=7)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT, "pipeline_retention.png"), dpi=150)
    plt.close(fig)


def write_report(summ, cons, cons_sig, harm, un_rows):
    L = []
    L.append("# OSCC 16S three-cohort preliminary report (genus level)")
    L.append("")
    L.append(f"Generated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}")
    L.append(f"Cohorts integrated: {', '.join(COHORTS)}")
    L.append("")
    L.append("## Cohorts and pipelines")
    L.append("")
    L.append("- PRJNA666746: **DADA2 ASV cohort** (SILVA 138.2), 50 complete tumor/adjacent pairs.")
    L.append("- PRJNA822685: **DADA2 ASV cohort** (SILVA 138.2), 56 complete tumor/adjacent pairs "
             "after depth/pairing QC of 81.")
    L.append("- PRJNA813034: **OTU sensitivity cohort** (97% OTU workflow; Greengenes 13_8 taxonomy), "
             "20 complete tumor/adjacent pairs.")
    L.append("- Two DADA2 (SILVA) cohorts + one Greengenes 13_8 OTU sensitivity cohort. "
             "PRJNA813034 is treated as sensitivity validation only, not as an equivalent DADA2 "
             "replication.")
    L.append("")
    L.append("## Methods")
    L.append("")
    L.append("- Each cohort was processed independently. Raw/ASV/OTU counts were **never pooled** "
             "across cohorts; cohort-specific effect estimates are retained throughout.")
    L.append("- Genus-level CLR transform with a single consistent pseudocount of 0.5.")
    L.append("- Per-genus paired Wilcoxon signed-rank test on paired Tumor-Matched_Normal CLR "
             "differences; BH-FDR across genera within each cohort.")
    L.append("- Genus harmonization uses exact genus-name matching after stripping database "
             "rank prefixes. SILVA 138.2 and Greengenes 13_8 largely agree at genus rank, but "
             "naming differences exist (e.g., SILVA lineage-group labels such as "
             "`Burkholderia-Caballeronia-Paraburkholderia`, `Christensenellaceae_R-7_group`, "
             "`[Eubacterium]_yurii_group` have no literal Greengenes counterpart; Greengenes has "
             "`Burkholderia`, `g__1-68`, candidate/environmental labels not present in SILVA). "
             "These naming differences are **not** interpreted as biological absence.")
    L.append("- Cross-cohort conclusions prefer genera detectable in >=2 cohorts (shared labels); "
             "unmatched/database-specific labels are catalogued separately in "
             "`harmonization_unmatched.tsv`.")
    L.append("- Ranking score weights (documented in `candidate_ranking.tsv`): direction "
             "consistency across cohorts, statistical support (q<0.1 and p<0.05), effect "
             "magnitude, prevalence/usable abundance, and robustness across the two DADA2 "
             "cohorts plus the Greengenes sensitivity cohort. This is candidate "
             "prioritisation, not causal proof.")
    L.append("")
    L.append("## Cohort summary")
    L.append("")
    L.append(summ.to_markdown(index=False))
    L.append("")
    L.append("## Genus-label sharing / harmonization")
    L.append("")
    L.append(f"- Shared genus labels detected in >=2 cohorts: {len(cons)} (see "
             "`cross_cohort_consistency.tsv`).")
    L.append(f"- Genus labels found in all three cohorts: {(cons['n_cohorts'] == 3).sum()} "
             "biological genera (plus the 'Unclassified' pseudo-label).")
    L.append(f"- Database-specific (unmatched) labels: {len(un_rows)} (see "
             "`harmonization_unmatched.tsv`).")
    L.append("")
    L.append("## Consistent genera")
    L.append("")
    L.append(f"- Genera with the same effect direction in >=2 cohorts and q<0.1 in >=1 cohort: "
             f"{len(cons_sig)} (see `consistent_genera.tsv`).")
    L.append("")
    if len(cons_sig):
        L.append(cons_sig.head(20)[
            ["genus", "n_cohorts", "n_DADA2_SILVA", "in_sensitivity_GG", "direction",
             "n_sig_q0.1", "min_q", "mean_prevalence", "mean_abs_effect",
             "cohort_effects"]].to_markdown(index=False))
        L.append("")
    L.append("## Top 25 candidate genera (three-cohort ranking)")
    L.append("")
    L.append(cons.head(25)[
        ["rank", "tier", "genus", "n_cohorts", "n_DADA2_SILVA", "in_sensitivity_GG",
         "direction", "n_sig_q0.1", "n_nominal_p0.05", "min_q", "mean_prevalence",
         "mean_abs_effect", "mean_median_effect", "score"]].to_markdown(index=False))
    L.append("")
    L.append("Tiers: A = all 3 cohorts, consistent; B = all 3 cohorts, discordant; "
             "C = both DADA2 cohorts consistent, Greengenes label not detected; "
             "D = both DADA2 cohorts discordant, GG not detected; E = one SILVA + GG "
             "consistent; F = two cohorts consistent; Z = discordant/partial.")
    L.append("")
    L.append("## Outputs")
    L.append("")
    L.append("- `combined/three_cohort/pipeline_status.tsv` - pipeline state")
    L.append("- `combined/three_cohort/cohort_summary.tsv` - cohort summary")
    L.append("- `combined/three_cohort/cohort_effects.tsv` - long-format cohort effect estimates")
    L.append("- `combined/three_cohort/genus_harmonization.tsv` - union genus catalogue")
    L.append("- `combined/three_cohort/harmonization_unmatched.tsv` - database-specific labels")
    L.append("- `combined/three_cohort/cross_cohort_consistency.tsv` - shared-genus consistency")
    L.append("- `combined/three_cohort/consistent_genera.tsv` - consistent, q<0.1 in >=1")
    L.append("- `combined/three_cohort/candidate_ranking.tsv` - ranked candidates")
    L.append("- `combined/three_cohort/effect_heatmap.png`, `top_candidates.png`, "
             "`direction_dotplot.png`, `pipeline_retention.png`")
    L.append("")
    L.append("## Limitations and caveats")
    L.append("")
    L.append("- PRJNA813034 is an OTU sensitivity cohort using a different denoising/OTU "
             "approach (97% OTU vs DADA2 ASV) and a different taxonomy database (Greengenes "
             "13_8 vs SILVA 138.2); effect signs there are corroborating, not equivalent.")
    L.append("- Taxonomy-database heterogeneity: a genus label absent from one database/cohort "
             "may reflect naming conventions, resolution, or primer/platform differences rather "
             "than biological absence.")
    L.append("- PRJNA813034 has 20 pairs only; **no genus reaches q<0.1** in that cohort "
             f"({int((pd.read_csv(os.path.join(ROOT, 'PRJNA813034', 'results', 'paired_genus_stats.tsv'), sep='\\t')['wilcoxon_p'] < P_THRESH).sum())} "
             "genera are nominal p<0.05). Its role is direction/magnitude corroboration.")
    L.append("- PRJNA666746 yields a very large number of q<0.1 genera (most very rare taxa); "
             "the CLR pseudocount can make sparse low-abundance taxa appear testable. "
             "Prevalence and magnitude are therefore reported for every candidate.")
    L.append("- In PRJNA813034 two tumor samples (SRR18236732, SRR18236733, T4_P04 and T4_P03) "
             "had essentially all reads unassignable at kingdom level by Greengenes "
             "(~123 and ~253 bacterial-mapped reads of ~35k); they remain in the 20 pairs as "
             "validated but add noise.")
    L.append("- Cross-study, cross-primer, cross-platform comparisons at genus level are "
             "exploratory. Associations are correlational; this analysis provides candidate "
             "prioritisation, not causal proof.")
    L.append("")

    with open(os.path.join(OUT, "preliminary_report.md"), "w") as f:
        f.write("\n".join(L))


if __name__ == "__main__":
    main()

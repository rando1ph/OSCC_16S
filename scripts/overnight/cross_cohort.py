#!/usr/bin/env python3
"""Cross-cohort genus-level integration, ranking, report and plots.

Usage: cross_cohort.py <combined_dir> <cohort1> <cohort2> [<cohort3> ...]
Each cohort path is <ROOT>/<COHORT>.
"""
import sys
import os
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.expanduser("~/OSCC_16S")
COMBINED = os.path.join(ROOT, "combined")
PSEUDOCOUNT = 0.5


SENSITIVITY_ROUTE = "OTU_VSEARCH_97"
PRIMARY_ROUTE = "DADA2_ASV"


def analysis_route(cohort, root, sensitivity):
    if cohort in sensitivity:
        return SENSITIVITY_ROUTE
    marker = os.path.join(root, cohort, ".otu_done")
    if os.path.exists(marker):
        return SENSITIVITY_ROUTE
    return PRIMARY_ROUTE


def read_stats(cohort_dir):
    path = os.path.join(cohort_dir, "results", "paired_genus_stats.tsv")
    if not os.path.exists(path):
        return None
    return pd.read_csv(path, sep="\t")


def load_retention(cohort_dir):
    p = os.path.join(cohort_dir, "processed", "retention.tsv")
    if os.path.exists(p):
        return pd.read_csv(p, sep="\t")
    return None


def load_depth(cohort_dir):
    p = os.path.join(cohort_dir, "results", "dada2_retention.tsv")
    if os.path.exists(p):
        return pd.read_csv(p, sep="\t", index_col=0)
    return None


def cohort_summary(cohort, cd, route):
    ret = load_retention(cd)
    dep = load_depth(cd)
    stat = read_stats(cd)
    row = {"cohort": cohort, "analysis_route": route}
    if ret is not None and len(ret):
        row["samples_processed"] = len(ret)
        row["mean_raw_pairs"] = round(ret["raw_pairs"].mean(), 1)
        row["mean_pretrim_pairs"] = round(ret["pretrim_pairs"].mean(), 1)
        row["mean_valid_primer_pairs"] = round(ret["valid_primer_pairs"].mean(), 1)
        row["mean_final_retained_pairs"] = round(ret["final_retained_pairs"].mean(), 1)
        row["mean_valid_fraction"] = round(ret["valid_fraction"].mean(), 4)
    if dep is not None and len(dep):
        row["samples_dada2"] = len(dep)
        row["mean_input_reads"] = round(dep["input_reads"].mean(), 1)
        row["mean_non_chimeric_reads"] = round(dep["non_chimeric_reads"].mean(), 1)
        row["low_depth_samples"] = int((dep["non_chimeric_reads"] < 2000).sum())
    if stat is not None:
        row["complete_pairs"] = int(stat["complete_pairs"].max())
        row["n_genera"] = len(stat)
        row["n_genera_sig_q0.1"] = int((stat["bh_fdr_q"] < 0.1).sum())
    return row


def main():
    args = sys.argv[1:]
    global COMBINED
    root = os.environ.get("CROSS_ROOT", ROOT)
    sensitivity = set()
    if args and args[0] == "--root":
        root = args[1]
        args = args[2:]
    if args and args[0] == "--sensitivity":
        sensitivity = set(args[1].split(","))
        args = args[2:]
    cohorts = args  # cohort names
    COMBINED = os.path.join(root, "combined")
    os.makedirs(COMBINED, exist_ok=True)

    summaries = []
    effects_long = []
    for c in cohorts:
        cd = os.path.join(root, c)
        stat = read_stats(cd)
        if stat is None:
            continue
        route = analysis_route(c, root, sensitivity)
        summaries.append(cohort_summary(c, cd, route))
        tmp = stat.copy()
        tmp["cohort"] = c
        tmp["analysis_route"] = route
        effects_long.append(tmp)

    if not effects_long:
        print("no cohort stats available")
        sys.exit(1)

    # cohort_summary.tsv
    summ = pd.DataFrame(summaries)
    summ.to_csv(os.path.join(COMBINED, "cohort_summary.tsv"), sep="\t", index=False)

    # cohort_effects.tsv (genus x cohort long)
    eff = pd.concat(effects_long, ignore_index=True)
    eff = eff[["cohort", "analysis_route", "genus", "complete_pairs",
               "median_tumor_relabund",
               "median_normal_relabund", "median_paired_clr_effect",
               "wilcoxon_p", "bh_fdr_q"]]
    eff.to_csv(os.path.join(COMBINED, "cohort_effects.tsv"), sep="\t",
               float_format="%.6g", index=False)

    # pivot for heatmap
    pivot = eff.pivot(index="genus", columns="cohort", values="median_paired_clr_effect")
    pvals = eff.pivot(index="genus", columns="cohort", values="wilcoxon_p")
    qvals = eff.pivot(index="genus", columns="cohort", values="bh_fdr_q")

    # consistent / ranking
    rows = []
    for genus, s in eff.groupby("genus"):
        present = s[s["median_paired_clr_effect"].notna()]
        if len(present) == 0:
            continue
        signs = np.sign(present["median_paired_clr_effect"])
        signs = signs.replace(0, np.nan).dropna()
        n_cohorts = len(present)
        if len(signs) >= 2:
            dir_cons = int((signs == signs.iloc[0]).all())
            n_same_sign = int((signs == signs.iloc[0]).sum())
            direction = "up" if signs.iloc[0] > 0 else "down"
        else:
            dir_cons = 0
            n_same_sign = len(signs)
            direction = "up" if len(signs) == 1 and signs.iloc[0] > 0 else (
                "down" if len(signs) == 1 else "na")
        n_sig = int((present["bh_fdr_q"] < 0.1).sum())
        mean_abs = float(np.abs(present["median_paired_clr_effect"]).mean())
        mean_eff = float(present["median_paired_clr_effect"].mean())
        min_q = float(present["bh_fdr_q"].min())
        consistency = n_same_sign / n_cohorts
        support = 1.0 + n_sig
        magnitude = 1.0 + mean_abs
        score = consistency * support * magnitude
        rows.append({
            "genus": genus,
            "n_cohorts": n_cohorts,
            "n_same_sign": n_same_sign,
            "direction_consistent": dir_cons,
            "direction": direction,
            "n_sig_q0.1": n_sig,
            "min_q": min_q,
            "mean_abs_effect": mean_abs,
            "mean_median_effect": mean_eff,
            "score": score,
        })

    rank = pd.DataFrame(rows)
    rank = rank.sort_values(
        ["direction_consistent", "n_sig_q0.1", "n_cohorts", "score"],
        ascending=[False, False, False, False]).reset_index(drop=True)
    rank.to_csv(os.path.join(COMBINED, "candidate_ranking.tsv"), sep="\t",
                float_format="%.6g", index=False)

    # consistent genera
    cons = rank[(rank["n_cohorts"] >= 2) & (rank["direction_consistent"] == 1)
                & (rank["n_sig_q0.1"] >= 1)]
    # add cohort-level effect columns
    cons_out = []
    for _, r in cons.iterrows():
        g = r["genus"]
        sub = eff[eff["genus"] == g]
        cohorts_eff = ";".join(
            f"{co}{'[OTU]' if (sub.loc[i, 'analysis_route'] == SENSITIVITY_ROUTE) else ''}({ef:+.3g})"
            for i, co, ef in zip(
                sub.index, sub["cohort"], sub["median_paired_clr_effect"]))
        cons_out.append({
            "genus": g, "n_cohorts": r["n_cohorts"], "direction": r["direction"],
            "min_q": r["min_q"], "n_sig_q0.1": r["n_sig_q0.1"],
            "mean_abs_effect": r["mean_abs_effect"],
            "mean_median_effect": r["mean_median_effect"],
            "cohort_effects": cohorts_eff,
        })
    if cons_out:
        pd.DataFrame(cons_out).to_csv(
            os.path.join(COMBINED, "consistent_genera.tsv"), sep="\t",
            float_format="%.6g", index=False)
    else:
        pd.DataFrame(columns=["genus", "n_cohorts", "direction", "min_q",
                              "n_sig_q0.1", "mean_abs_effect",
                              "mean_median_effect", "cohort_effects"]).to_csv(
            os.path.join(COMBINED, "consistent_genera.tsv"), sep="\t", index=False)

    make_plots(eff, pivot, rank, cohorts, root)

    write_report(cohorts, summ, rank, cons_out)

    print("cross-cohort outputs written to", COMBINED)


def make_plots(eff, pivot, rank, cohorts, root):
    plt.style.use("default")

    # heatmap
    order = rank["genus"].head(25).tolist()
    hmap = pivot.loc[pivot.index.intersection(order)].T
    hmap = hmap[[g for g in order if g in hmap.columns]]
    fig, ax = plt.subplots(figsize=(max(4, 0.5 * len(hmap.columns) + 2),
                                    0.5 * len(hmap.index) + 2))
    im = ax.imshow(hmap.values, cmap="RdBu_r", aspect="auto", vmin=-np.nanmax(
        np.abs(hmap.values[np.isfinite(hmap.values)])) if np.isfinite(
            hmap.values).any() else -1,
        vmax=np.nanmax(np.abs(hmap.values[np.isfinite(hmap.values)])) if np.isfinite(
            hmap.values).any() else 1)
    ax.set_xticks(range(len(hmap.columns)))
    ax.set_xticklabels(hmap.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(hmap.index)))
    ax.set_yticklabels(hmap.index)
    ax.set_title("Median paired CLR Tumor - Matched_Normal effect")
    plt.colorbar(im, ax=ax, label="median paired CLR effect")
    plt.tight_layout()
    plt.savefig(os.path.join(COMBINED, "effect_heatmap.png"), dpi=150)
    plt.close(fig)

    # top candidates bar
    top = rank.head(20)
    fig, ax = plt.subplots(figsize=(10, 7))
    colors = ["#d62728" if d > 0 else "#1f77b4" for d in top["mean_median_effect"]]
    ax.barh(top["genus"], top["mean_median_effect"], color=colors)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("mean median paired CLR effect (Tumor - Matched_Normal)")
    ax.set_title("Top ranked candidate genera (cross-cohort)")
    plt.tight_layout()
    plt.savefig(os.path.join(COMBINED, "top_candidates.png"), dpi=150)
    plt.close(fig)

    # retention plot
    fig, ax = plt.subplots(figsize=(10, 6))
    stage_labels = ["raw", "pretrim\n(cap 150k)", "valid primer",
                    "final\n(cap 100k)", "dada2\nnon-chimeric"]
    all_means = []
    cohort_colors = plt.cm.Set1(np.linspace(0, 1, len(cohorts)))
    for i, c in enumerate(cohorts):
        cd = os.path.join(root, c)
        ret = load_retention(cd)
        dep = load_depth(cd)
        if ret is None:
            continue
        raw = ret["raw_pairs"].mean()
        pre = ret["pretrim_pairs"].mean()
        valid = ret["valid_primer_pairs"].mean()
        final = ret["final_retained_pairs"].mean()
        dada2 = dep["non_chimeric_reads"].mean() if dep is not None else np.nan
        means = [raw, pre, valid, final, dada2]
        all_means.append((c, means))
        ax.plot(range(len(stage_labels)), means, "o-", color=cohort_colors[i],
                label=c, lw=2)
        for x, y in zip(range(len(stage_labels)), means):
            if not np.isnan(y):
                ax.text(x, y, f"{int(y):,}", fontsize=7, ha="center", va="bottom")
    ax.set_xticks(range(len(stage_labels)))
    ax.set_xticklabels(stage_labels)
    ax.set_yscale("log")
    ax.set_ylabel("mean reads/pairs per sample (log scale)")
    ax.set_title("Pipeline retention across stages")
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(COMBINED, "pipeline_retention.png"), dpi=150)
    plt.close(fig)


def write_report(cohorts, summ, rank, cons_out):
    lines = []
    lines.append("# OSCC 16S preliminary multi-cohort report")
    lines.append("")
    lines.append(f"Generated: {pd.Timestamp.now()}")
    lines.append(f"Cohorts integrated: {', '.join(cohorts)}")
    lines.append("")
    lines.append("## Methods")
    lines.append("")
    lines.append("- Each cohort processed independently; SILVA 138.2 taxonomy; genus collapse.")
    lines.append("- Primary cohorts (PRJNA666746, PRJNA822685): DADA2 denoising (ASV route).")
    lines.append("- PRJNA813034: OTU-based sensitivity analysis (VSEARCH merge-pairs, "
                 "dereplication, de novo chimera filtering and 97% de novo OTU clustering, "
                 "following the source study's USEARCH 97% OTU methodology). DADA2 is not "
                 "applicable to this cohort because its discrete quality scores prevent "
                 "error-model estimation. PRJNA813034 results are treated as sensitivity "
                 "validation only.")
    lines.append(f"- Genera marked [OTU] in consistent_genera.tsv are from the sensitivity cohort.")
    lines.append(f"- Genus-level CLR transform applied with a single consistent pseudocount of {PSEUDOCOUNT}.")
    lines.append("- Per-genus paired Wilcoxon signed-rank test on Tumor-Matched_Normal paired CLR differences; BH-FDR across genera.")
    lines.append("- Samples with <2000 non-chimeric reads flagged LOW_DEPTH; if one member of a pair is excluded, the whole pair is removed.")
    lines.append("- Genera ranked primarily by cross-cohort effect-direction consistency, statistical support (q<0.1), and effect magnitude.")
    lines.append("- Association, not causation.")
    lines.append("")
    lines.append("## Cohort summary")
    lines.append("")
    lines.append(summ.to_markdown(index=False))
    lines.append("")
    if len(rank):
        lines.append("## Top 20 candidate genera")
        lines.append("")
        lines.append(rank.head(20).to_markdown(index=False))
    lines.append("")
    if cons_out:
        lines.append("## Consistent genera (same direction in >=2 cohorts, q<0.1 in >=1)")
        lines.append("")
        lines.append(pd.DataFrame(cons_out).to_markdown(index=False))
    lines.append("")
    lines.append("## Caveats")
    lines.append("")
    lines.append(f"- PRJNA813034 used an OTU (97% VSEARCH) route and is therefore a sensitivity "
                 "validation rather than an equivalent DADA2 replication; effect signs for this "
                 "cohort should be interpreted accordingly.")
    lines.append("- Cross-study comparisons at genus level are exploratory; primer/platform differences remain.")
    lines.append("- Raw and processed reads are excluded from Git; all tables/plots here are reproducible via scripts/overnight/.")
    lines.append("")

    with open(os.path.join(COMBINED, "preliminary_report.md"), "w") as f:
        f.write("\n".join(lines))


if __name__ == "__main__":
    main()

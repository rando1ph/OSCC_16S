#!/usr/bin/env python3
"""Leave-one-cohort-out regularized logistic regression.

Usage: ml_analysis.py <cohort1> <cohort2> <cohort3>
Requires all three cohorts to have genus CLR matrices and sample status.
Reports held-out ROC AUC and coefficient stability. Failure here must not
invalidate paired statistics.

Outputs: combined/ml_leave_one_cohort_out.tsv
"""
import sys
import os
import json
from collections import defaultdict

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score

ROOT = os.path.expanduser("~/OSCC_16S")
COMBINED = os.path.join(ROOT, "combined")


def load_clr(cohort, root):
    p = os.path.join(root, cohort, "results", "clr_matrix.tsv")
    if not os.path.exists(p):
        return None
    return pd.read_csv(p, sep="\t", index_col=0)


def load_status(cohort, root):
    p = os.path.join(root, cohort, "results", "sample_status.tsv")
    if not os.path.exists(p):
        return None
    return pd.read_csv(p, sep="\t")


def load_metadata(cohort, root):
    p = os.path.join(root, cohort, "metadata", "analysis_metadata.tsv")
    return pd.read_csv(p, sep="\t")[["sample-id", "group"]]


def main():
    args = sys.argv[1:]
    root = os.environ.get("CROSS_ROOT", ROOT)
    if args and args[0] == "--root":
        root = args[1]
        args = args[2:]
    cohorts = args
    COMBINED = os.path.join(root, "combined")
    clrs = {c: load_clr(c, root) for c in cohorts}
    stats = {c: load_status(c, root) for c in cohorts}
    metas = {c: load_metadata(c, root) for c in cohorts}

    if any(v is None for v in clrs.values()) or any(v is None for v in stats.values()):
        print("ML SKIPPED: missing clr_matrix or sample_status in one or more cohorts")
        sys.exit(0)

    rows = []
    coef_sign = defaultdict(list)

    for held in cohorts:
        train = [c for c in cohorts if c != held]
        Xtr_parts = []
        ytr = []
        for c in train:
            status = stats[c]
            keep = status[status["status"] == "included"]
            if len(keep) == 0:
                continue
            meta = metas[c].set_index("sample-id").loc[keep["sample"]]
            clr = clrs[c].loc[keep["sample"]]
            Xtr_parts.append(clr)
            ytr.extend((meta["group"] == "Tumor").astype(int).tolist())
        if len(Xtr_parts) < 2 or len(ytr) < 4:
            rows.append({"held_out": held, "roc_auc": np.nan,
                         "note": "insufficient training data"})
            continue
        Xtr = pd.concat(Xtr_parts)
        feature_cols = Xtr.columns.tolist()
        Xtr = Xtr.values
        ytr = np.asarray(ytr)

        scaler = StandardScaler().fit(Xtr)
        Xtr_s = scaler.transform(Xtr)
        clf = LogisticRegression(C=0.1, penalty="l2", solver="lbfgs",
                                 max_iter=2000)
        clf.fit(Xtr_s, ytr)

        # held-out evaluation (standardize with training-only scaler)
        status = stats[held]
        keep = status[status["status"] == "included"]
        if len(keep) == 0:
            rows.append({"held_out": held, "roc_auc": np.nan,
                         "note": "no included held-out samples"})
            continue
        meta = metas[held].set_index("sample-id").loc[keep["sample"]]
        Xte = clrs[held].loc[keep["sample"]].reindex(columns=feature_cols).fillna(0.0)
        yte = (meta["group"] == "Tumor").astype(int).values
        Xte_s = scaler.transform(Xte.values)

        if len(np.unique(yte)) < 2:
            rows.append({"held_out": held, "n_test": len(yte),
                         "roc_auc": np.nan, "note": "held-out set single-class"})
        else:
            auc = roc_auc_score(yte, clf.predict_proba(Xte_s)[:, 1])
            rows.append({
                "held_out": held,
                "train_cohorts": "+".join(train),
                "n_train": len(ytr),
                "n_test": len(yte),
                "n_features": len(feature_cols),
                "roc_auc": auc,
                "note": "",
            })

        for f, coef in zip(feature_cols, clf.coef_[0]):
            coef_sign[f].append(np.sign(coef))

    df = pd.DataFrame(rows)
    # coefficient stability across folds
    stab = []
    for f, signs in coef_sign.items():
        if len(signs) >= 2:
            same = sum(1 for s in signs if s == signs[0])
            stab.append({"feature": f, "n_folds": len(signs),
                         "fraction_same_sign": same / len(signs)})
    stab_df = pd.DataFrame(stab) if stab else pd.DataFrame(
        columns=["feature", "n_folds", "fraction_same_sign"])
    mean_stab = (stab_df["fraction_same_sign"].mean()
                 if len(stab_df) else np.nan)

    summary = {
        "n_folds": len(df),
        "mean_roc_auc": round(float(df["roc_auc"].mean()), 4) if df["roc_auc"].notna().any() else None,
        "median_roc_auc": round(float(df["roc_auc"].median()), 4) if df["roc_auc"].notna().any() else None,
        "mean_coef_sign_consistency": round(float(mean_stab), 4) if mean_stab == mean_stab else None,
        "n_features_with_coefs": len(stab_df),
        "pseudocount": 0.5,
        "model": "LogisticRegression(C=0.1, l2), StandardScaler fit on training only",
    }

    with open(os.path.join(COMBINED, "ml_leave_one_cohort_out.tsv"), "w") as f:
        df.to_csv(f, sep="\t", index=False, float_format="%.6g")
    with open(os.path.join(COMBINED, "ml_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print("ML results written; summary:", summary)


if __name__ == "__main__":
    main()

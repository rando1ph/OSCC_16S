import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

ranking = pd.read_csv(
    "combined/three_cohort/candidate_ranking.tsv",
    sep="\t"
)

# 前15名
df = ranking.sort_values("rank").head(15).copy()

# 去掉 g__ 前缀
df["genus_label"] = df["genus"].str.replace(r"^g__", "", regex=True)

# rank 1 放最上方
df = df.iloc[::-1].reset_index(drop=True)

y = np.arange(len(df))
h = 0.35

fig, ax = plt.subplots(figsize=(9, 8))

ax.barh(
    y - h/2,
    df["n_sig_q0.1"],
    height=h,
    label="Cohorts with q < 0.1"
)

ax.barh(
    y + h/2,
    df["n_nominal_p0.05"],
    height=h,
    label="Cohorts with p < 0.05"
)

ax.set_yticks(y)
ax.set_yticklabels(df["genus_label"])

ax.set_xlim(0, 3.2)
ax.set_xticks([0, 1, 2, 3])

ax.set_xlabel("Number of cohorts")
ax.set_ylabel("Genus")
ax.set_title("Statistical support across cohorts")

ax.legend()
ax.grid(axis="x", alpha=0.25)

plt.tight_layout()

plt.savefig(
    "combined/three_cohort/statistical_support.png",
    dpi=300,
    bbox_inches="tight"
)

plt.savefig(
    "combined/three_cohort/statistical_support.pdf",
    bbox_inches="tight"
)

print("Saved:")
print("combined/three_cohort/statistical_support.png")
print("combined/three_cohort/statistical_support.pdf")

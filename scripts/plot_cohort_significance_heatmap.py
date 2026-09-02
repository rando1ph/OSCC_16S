import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ranking = pd.read_csv(
    "combined/three_cohort/candidate_ranking.tsv",
    sep="\t"
)

effects = pd.read_csv(
    "combined/three_cohort/cohort_effects.tsv",
    sep="\t"
)

# Top 15
top = ranking.sort_values("rank").head(15).copy()
top_genera = top["genus"].tolist()

df = effects[effects["genus"].isin(top_genera)].copy()

# 显示名去掉 g__
df["genus_label"] = df["genus"].str.replace(r"^g__", "", regex=True)

# cohort顺序固定
cohorts = [
    "PRJNA666746",
    "PRJNA822685",
    "PRJNA813034"
]

# genus顺序按ranking，rank1在上面
genus_order = (
    top["genus"]
    .str.replace(r"^g__", "", regex=True)
    .tolist()
)

# 0 = 不显著
# 1 = p < 0.05
# 2 = q < 0.1
def support_level(row):
    if row["bh_fdr_q"] < 0.1:
        return 2
    elif row["wilcoxon_p"] < 0.05:
        return 1
    else:
        return 0

df["support"] = df.apply(support_level, axis=1)

matrix = (
    df.pivot(
        index="genus_label",
        columns="cohort",
        values="support"
    )
    .reindex(index=genus_order, columns=cohorts)
)

fig, ax = plt.subplots(figsize=(7, 9))

im = ax.imshow(
    matrix.values,
    aspect="auto",
    vmin=0,
    vmax=2
)

ax.set_xticks(np.arange(len(cohorts)))
ax.set_xticklabels(cohorts)

ax.set_yticks(np.arange(len(genus_order)))
ax.set_yticklabels(genus_order)

ax.set_title("Cohort-specific statistical support")

# 在格子里直接写结果，比猜颜色文明一点
for i in range(len(genus_order)):
    for j in range(len(cohorts)):
        value = matrix.iloc[i, j]

        if pd.isna(value):
            text = "NA"
        elif value == 2:
            text = "q<0.1"
        elif value == 1:
            text = "p<0.05"
        else:
            text = "NS"

        ax.text(
            j,
            i,
            text,
            ha="center",
            va="center",
            fontsize=8
        )

ax.set_xlabel("Cohort")
ax.set_ylabel("Genus")

plt.tight_layout()

plt.savefig(
    "combined/three_cohort/cohort_significance_heatmap.png",
    dpi=300,
    bbox_inches="tight"
)

plt.savefig(
    "combined/three_cohort/cohort_significance_heatmap.pdf",
    bbox_inches="tight"
)

print("Saved:")
print("combined/three_cohort/cohort_significance_heatmap.png")
print("combined/three_cohort/cohort_significance_heatmap.pdf")

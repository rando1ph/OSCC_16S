import pandas as pd
import matplotlib.pyplot as plt

# 读取候选排名和每个 cohort 的 effect
ranking = pd.read_csv(
    "combined/three_cohort/candidate_ranking.tsv",
    sep="\t"
)

effects = pd.read_csv(
    "combined/three_cohort/cohort_effects.tsv",
    sep="\t"
)

# 只取排名前15
top = ranking.sort_values("rank").head(15).copy()
top_genera = top["genus"].tolist()

df = effects[effects["genus"].isin(top_genera)].copy()

# 去掉 g__ 前缀，画图好看一点
df["genus_label"] = df["genus"].str.replace(r"^g__", "", regex=True)

# 按 candidate ranking 排序
order = (
    top["genus"]
    .str.replace(r"^g__", "", regex=True)
    .tolist()
)

# 反转，使 rank 1 在最上方
order = order[::-1]

cohorts = [
    "PRJNA666746",
    "PRJNA822685",
    "PRJNA813034"
]

fig, ax = plt.subplots(figsize=(9, 8))

for cohort in cohorts:
    sub = df[df["cohort"] == cohort]

    ax.scatter(
        sub["median_paired_clr_effect"],
        sub["genus_label"],
        label=cohort,
        s=55,
        alpha=0.8
    )

# 0线：左边 tumor-depleted，右边 tumor-enriched
ax.axvline(0, linewidth=1)

ax.set_yticks(range(len(order)))
ax.set_yticklabels(order)

ax.set_xlabel("Median paired CLR effect (Tumor − Matched Normal)")
ax.set_ylabel("Genus")
ax.set_title("Cross-cohort effects of top 15 candidate genera")

ax.legend(title="Cohort")
ax.grid(axis="x", alpha=0.25)

plt.tight_layout()

plt.savefig(
    "combined/three_cohort/cross_cohort_effect_dotplot.png",
    dpi=300,
    bbox_inches="tight"
)

plt.savefig(
    "combined/three_cohort/cross_cohort_effect_dotplot.pdf",
    bbox_inches="tight"
)

print("Saved:")
print("combined/three_cohort/cross_cohort_effect_dotplot.png")
print("combined/three_cohort/cross_cohort_effect_dotplot.pdf")

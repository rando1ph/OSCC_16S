import pandas as pd
import matplotlib.pyplot as plt

ranking = pd.read_csv(
    "combined/three_cohort/candidate_ranking.tsv",
    sep="\t"
)

# 只画前15名
df = ranking.sort_values("rank").head(15).copy()

df["genus_label"] = df["genus"].str.replace(r"^g__", "", regex=True)

fig, ax = plt.subplots(figsize=(9, 7))

ax.scatter(
    df["mean_prevalence"],
    df["mean_abs_effect"],
    s=90,
    alpha=0.8
)

for _, row in df.iterrows():
    ax.annotate(
        row["genus_label"],
        (row["mean_prevalence"], row["mean_abs_effect"]),
        xytext=(5, 4),
        textcoords="offset points",
        fontsize=8
    )

ax.set_xlabel("Mean prevalence across cohorts")
ax.set_ylabel("Mean absolute CLR effect")
ax.set_title("Prevalence and effect size of top candidate genera")

ax.grid(alpha=0.25)

plt.tight_layout()

plt.savefig(
    "combined/three_cohort/prevalence_effect.png",
    dpi=300,
    bbox_inches="tight"
)

plt.savefig(
    "combined/three_cohort/prevalence_effect.pdf",
    bbox_inches="tight"
)

print("Saved:")
print("combined/three_cohort/prevalence_effect.png")
print("combined/three_cohort/prevalence_effect.pdf")

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# Input / output
# ============================================================

input_file = "abundance_matrix.tsv"

results_file = "differential_results.tsv"
barplot_file = "abundance_barplot.png"
heatmap_file = "abundance_heatmap.png"


# ============================================================
# Read abundance matrix
# ============================================================

df = pd.read_csv(input_file, sep="\t")

samples = [c for c in df.columns if c != "Taxon"]

if len(samples) < 2:
    raise ValueError("Need at least two samples in the abundance matrix.")

# Convert to numeric, normalize to relative abundance per sample
for s in samples:
    df[s] = pd.to_numeric(df[s], errors="coerce").fillna(0)
    total = df[s].sum()
    if total > 0:
        df[s] = df[s] / total


# ============================================================
# Per-taxon summary stats across all samples
# ============================================================

pseudocount = 1e-9

df["mean_abundance"] = df[samples].mean(axis=1)
df["std_abundance"] = df[samples].std(axis=1)
df["max_abundance"] = df[samples].max(axis=1)
df["min_abundance"] = df[samples].min(axis=1)

df["fold_change_max_min"] = (
    (df["max_abundance"] + pseudocount) / (df["min_abundance"] + pseudocount)
)
df["log2_fold_change_max_min"] = np.log2(df["fold_change_max_min"])

df["dominant_sample"] = df[samples].idxmax(axis=1)


# ============================================================
# Save numerical results
# ============================================================

results_cols = (
    ["Taxon"]
    + samples
    + [
        "mean_abundance",
        "std_abundance",
        "fold_change_max_min",
        "log2_fold_change_max_min",
        "dominant_sample",
    ]
)

results = df[results_cols].copy()
results = results.sort_values("std_abundance", ascending=False)
results.to_csv(results_file, sep="\t", index=False)


# ============================================================
# Select top 15 most abundant taxa (by total across samples)
# ============================================================

df["total_abundance"] = df[samples].sum(axis=1)

plot_df = (
    df.sort_values("total_abundance", ascending=False)
    .head(15)
    .copy()
)

taxa = plot_df["Taxon"].astype(str).tolist()

x = np.arange(len(taxa))
n_samples = len(samples)
width = 0.8 / n_samples


# ============================================================
# BARPLOT — grouped bars, one group per sample
# ============================================================

fig, ax = plt.subplots(figsize=(14, 8))

for i, sample in enumerate(samples):
    offset = (i - (n_samples - 1) / 2) * width
    ax.bar(x + offset, plot_df[sample], width, label=sample)

ax.set_xlabel("Taxon")
ax.set_ylabel("Relative abundance")
ax.set_title("Relative Abundance by Sample (top 15 taxa)")

ax.set_xticks(x)
ax.set_xticklabels(taxa, rotation=45, ha="right")

ax.legend()

plt.tight_layout()
plt.savefig(barplot_file, dpi=300, bbox_inches="tight")
plt.close()


# ============================================================
# HEATMAP
# ============================================================

heatmap_df = plot_df.set_index("Taxon")[samples]

fig, ax = plt.subplots(figsize=(max(6, n_samples * 2), 10))

im = ax.imshow(heatmap_df.values, aspect="auto")

ax.set_xticks(np.arange(len(heatmap_df.columns)))
ax.set_xticklabels(heatmap_df.columns, rotation=45, ha="right")

ax.set_yticks(np.arange(len(heatmap_df.index)))
ax.set_yticklabels(heatmap_df.index)

ax.set_xlabel("Sample")
ax.set_ylabel("Taxon")
ax.set_title("Relative Abundance Heatmap (top 15 taxa)")

fig.colorbar(im, ax=ax, label="Relative abundance")

plt.tight_layout()
plt.savefig(heatmap_file, dpi=300, bbox_inches="tight")
plt.close()


# ============================================================
# Summary
# ============================================================

print()
print("==============================================")
print("DIFFERENTIAL ABUNDANCE")
print("==============================================")
print()
print(f"Samples : {', '.join(samples)}")
print(f"Taxa    : {len(df)}")
print()
print(f"Results : {results_file}")
print(f"Barplot : {barplot_file}")
print(f"Heatmap : {heatmap_file}")
print()
print("Differential abundance analysis complete.")
print("==============================================")

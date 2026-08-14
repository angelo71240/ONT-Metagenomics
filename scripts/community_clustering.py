import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.decomposition import PCA


# ------------------------------------------------------------
# Input
# ------------------------------------------------------------

input_file = "abundance_matrix.tsv"

df = pd.read_csv(input_file, sep="\t")

# First column = taxon names
df = df.set_index("Taxon")

# Samples are columns
samples = df.columns.tolist()

# Convert to numeric
abundance = df.apply(pd.to_numeric, errors="coerce").fillna(0)

# Remove taxa absent from all samples
abundance = abundance.loc[abundance.sum(axis=1) > 0]


# ------------------------------------------------------------
# Relative abundance normalization
# ------------------------------------------------------------

relative_abundance = abundance.div(
    abundance.sum(axis=0),
    axis=1
)

# Samples × taxa
sample_matrix = relative_abundance.T


# ------------------------------------------------------------
# Bray-Curtis distance
# ------------------------------------------------------------

bray_curtis = pdist(
    sample_matrix.values,
    metric="braycurtis"
)

bray_matrix = pd.DataFrame(
    squareform(bray_curtis),
    index=samples,
    columns=samples
)

bray_matrix.to_csv(
    "bray_curtis_distance.tsv",
    sep="\t"
)


# ------------------------------------------------------------
# PCoA
# ------------------------------------------------------------

distance_matrix = squareform(bray_curtis)

n = distance_matrix.shape[0]

centering = np.eye(n) - np.ones((n, n)) / n

B = -0.5 * centering @ (distance_matrix ** 2) @ centering

eigenvalues, eigenvectors = np.linalg.eigh(B)

order = np.argsort(eigenvalues)[::-1]

eigenvalues = eigenvalues[order]
eigenvectors = eigenvectors[:, order]

positive = eigenvalues > 0

eigenvalues_positive = eigenvalues[positive]
eigenvectors_positive = eigenvectors[:, positive]

if len(eigenvalues_positive) >= 2:

    coordinates = (
        eigenvectors_positive[:, :2]
        * np.sqrt(eigenvalues_positive[:2])
    )

    explained = (
        eigenvalues_positive[:2]
        / eigenvalues_positive.sum()
        * 100
    )

    pcoa_df = pd.DataFrame(
        coordinates,
        index=samples,
        columns=["PCoA1", "PCoA2"]
    )

    pcoa_df.to_csv(
        "pcoa_coordinates.tsv",
        sep="\t"
    )

    plt.figure(figsize=(8, 6))

    plt.scatter(
        pcoa_df["PCoA1"],
        pcoa_df["PCoA2"],
        s=100
    )

    for sample in samples:
        plt.text(
            pcoa_df.loc[sample, "PCoA1"],
            pcoa_df.loc[sample, "PCoA2"],
            sample
        )

    plt.xlabel(f"PCoA1 ({explained[0]:.2f}%)")
    plt.ylabel(f"PCoA2 ({explained[1]:.2f}%)")
    plt.title("PCoA Based on Bray-Curtis Distance")

    plt.tight_layout()
    plt.savefig(
        "pcoa_bray_curtis.png",
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()


# ------------------------------------------------------------
# Hierarchical clustering
# ------------------------------------------------------------

if len(samples) >= 2:

    linkage_matrix = linkage(
        bray_curtis,
        method="average"
    )

    plt.figure(figsize=(10, 6))

    dendrogram(
        linkage_matrix,
        labels=samples
    )

    plt.title("Hierarchical Clustering")
    plt.xlabel("Samples")
    plt.ylabel("Bray-Curtis Distance")

    plt.tight_layout()

    plt.savefig(
        "hierarchical_clustering.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()


# ------------------------------------------------------------
# PCA
# ------------------------------------------------------------

if len(samples) >= 2:

    pca = PCA(n_components=min(2, len(samples)))

    pca_coordinates = pca.fit_transform(sample_matrix.values)

    pca_df = pd.DataFrame(
        pca_coordinates,
        index=samples,
        columns=[
            f"PC{i+1}"
            for i in range(pca_coordinates.shape[1])
        ]
    )

    pca_df.to_csv(
        "pca_coordinates.tsv",
        sep="\t"
    )

    if pca_coordinates.shape[1] >= 2:

        plt.figure(figsize=(8, 6))

        plt.scatter(
            pca_df["PC1"],
            pca_df["PC2"],
            s=100
        )

        for sample in samples:
            plt.text(
                pca_df.loc[sample, "PC1"],
                pca_df.loc[sample, "PC2"],
                sample
            )

        explained_pca = pca.explained_variance_ratio_ * 100

        plt.xlabel(f"PC1 ({explained_pca[0]:.2f}%)")
        plt.ylabel(f"PC2 ({explained_pca[1]:.2f}%)")
        plt.title("PCA of Community Composition")

        plt.tight_layout()

        plt.savefig(
            "pca.png",
            dpi=300,
            bbox_inches="tight"
        )

        plt.close()


# ------------------------------------------------------------
# Heatmap
# ------------------------------------------------------------

# Keep top 20 taxa by total abundance
top_taxa = (
    relative_abundance.sum(axis=1)
    .sort_values(ascending=False)
    .head(20)
    .index
)

heatmap_df = relative_abundance.loc[top_taxa]

plt.figure(
    figsize=(
        max(8, len(samples) * 2),
        max(6, len(top_taxa) * 0.35)
    )
)

plt.imshow(
    heatmap_df.values,
    aspect="auto",
    interpolation="nearest"
)

plt.colorbar(
    label="Relative abundance"
)

plt.xticks(
    range(len(samples)),
    samples,
    rotation=45,
    ha="right"
)

plt.yticks(
    range(len(top_taxa)),
    top_taxa
)

plt.xlabel("Samples")
plt.ylabel("Taxa")
plt.title("Community Composition Heatmap")

plt.tight_layout()

plt.savefig(
    "community_heatmap.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()


# ------------------------------------------------------------
# Finish
# ------------------------------------------------------------

print("==============================================")
print("COMMUNITY CLUSTERING")
print("==============================================")
print(f"Samples : {len(samples)}")
print(f"Taxa    : {len(abundance)}")
print()
print("Generated:")
print("  bray_curtis_distance.tsv")
print("  pcoa_coordinates.tsv")
print("  pcoa_bray_curtis.png")
print("  hierarchical_clustering.png")
print("  pca_coordinates.tsv")
print("  pca.png")
print("  community_heatmap.png")
print("==============================================")

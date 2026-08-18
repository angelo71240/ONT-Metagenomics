# ONT Shotgun Metagenomics - Step 5 & 6: Differential Abundance and Community Clustering

A Nextflow (DSL2) pipeline that takes filtered Oxford Nanopore shotgun metagenomic
reads through taxonomic profiling (Kraken2 + Bracken), builds a combined abundance
matrix, and runs two downstream analyses:

- **Step 5 – Differential Abundance**: which taxa differ in relative abundance between samples
- **Step 6 – Community Clustering**: how similar/dissimilar the samples' microbial communities are (Bray-Curtis, PCoA, hierarchical clustering, PCA, heatmap)

This repo covers Steps 5-6 only. Steps 1-4 (assembly, gene prediction, binning,
read-based taxonomy) are handled in a separate pipeline.

## Workflow

```
Raw ONT FASTQ (per sample)
        │
        ▼
   NanoPlot (pre-filter QC)
        │
        ▼
   Filtlong (length/quality filtering)
        │
        ▼
   NanoPlot (post-filter QC)
        │
        ▼
   MultiQC (combined QC report)
        │
        ▼
   Kraken2 (taxonomic classification)
        │
        ▼
   Bracken (species-level abundance estimation)
        │
        ▼
   Abundance matrix (taxa × samples)
        │
        ├──▶ Step 5: Differential Abundance
        │      ├── differential_results.tsv
        │      ├── abundance_barplot.png
        │      └── abundance_heatmap.png
        │
        └──▶ Step 6: Community Clustering
               ├── bray_curtis_distance.tsv
               ├── pcoa_coordinates.tsv / pcoa_bray_curtis.png
               ├── hierarchical_clustering.png
               ├── pca_coordinates.tsv / pca.png
               └── community_heatmap.png
```

## Repo structure

```
.
├── main.nf                       # Pipeline entry point
├── nextflow.config                # Per-process conda environments
├── modules/
│   ├── nanoplot.nf
│   ├── filtlong.nf
│   ├── multiqc.nf
│   ├── kraken2.nf
│   ├── bracken.nf
│   ├── abundance_matrix.nf
│   ├── differential_abundance.nf
│   └── community_clustering.nf
├── scripts/
│   ├── make_abundance_matrix.py
│   ├── differential_abundance.py
│   └── community_clustering.py
├── data/                          # Input FASTQ (*.fastq.gz) — not tracked in git
└── results/                       # Pipeline output — not tracked in git
```

## Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 22.10, DSL2
- Conda/Miniconda, with the following environments (see [Setup](#setup)):
  - `nanoplot_env` – NanoPlot
  - `filtlong_env` – Filtlong
  - `kraken2_env` – Kraken2, Bracken
  - `community_env` – pandas, numpy, matplotlib, scipy, scikit-learn
- A Kraken2/Bracken reference database (path set in `params.kraken_db`)

## Setup

```bash
# NanoPlot
conda create -n nanoplot_env -c conda-forge -c bioconda nanoplot

# Filtlong
conda create -n filtlong_env -c conda-forge -c bioconda filtlong

# Kraken2 + Bracken (already used/tested for this project)
# conda env: kraken2_env

# Community clustering (Python stats/plotting)
conda create -n community_env python=3.11 pandas numpy matplotlib scipy scikit-learn -y
```

Each new shell needs Conda initialized once before `conda activate` works:

```bash
source /path/to/miniconda3/etc/profile.d/conda.sh
```

## Usage

1. Place input FASTQ files (`*.fastq.gz`) in `data/`
2. Set `params.kraken_db` in `main.nf` (or pass with `--kraken_db`) to your Kraken2/Bracken database path
3. Run the pipeline:

```bash
nextflow run main.nf -with-conda -resume
```

`-resume` lets Nextflow reuse already-completed steps (e.g. QC/Kraken2 output) when
adding a new sample instead of rerunning everything.

4. Outputs land in `results/`:

```
results/
├── nanoplot/
├── filtlong/
├── multiqc/
├── kraken2/
├── bracken/
├── abundance_matrix.tsv
├── differential/
│   ├── differential_results.tsv
│   ├── abundance_barplot.png
│   └── abundance_heatmap.png
└── community_clustering/
    ├── bray_curtis_distance.tsv
    ├── pcoa_coordinates.tsv
    ├── pcoa_bray_curtis.png
    ├── hierarchical_clustering.png
    ├── pca_coordinates.tsv
    ├── pca.png
    └── community_heatmap.png
```

## Method notes

- **Kraken2 + Bracken** were chosen over 16S-based tools because these are shotgun
  reads: every DNA fragment is sequenced, not a single marker gene, so read-based
  classification + abundance re-estimation is the appropriate approach.
- **Step 6 is called "Community Clustering," not "denoising."** Denoising (e.g.
  DADA2/Deblur) needs many reads piled up at the same genomic position to correct
  sequencing errors — that only happens with amplicon sequencing (16S), where PCR
  produces millions of copies of one gene. Shotgun reads come from random genomic
  locations and rarely overlap, so denoising isn't applicable; clustering samples by
  abundance profile (Bray-Curtis / PCoA / hierarchical clustering / PCA) is the shotgun
  equivalent.
- `fastp`/Trimmomatic were dropped from the QC stage in favor of NanoPlot + Filtlong,
  which are built for ONT long reads rather than Illumina short reads.

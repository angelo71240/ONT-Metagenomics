# ONT Shotgun Metagenomics - Step 5 & 6: Differential Abundance and Community Clustering

A Nextflow (DSL2) pipeline that takes filtered Oxford Nanopore shotgun metagenomic
reads through taxonomic profiling (Kraken2 + Bracken), builds a combined abundance
matrix, and runs two downstream analyses:

- **Step 5 - Differential Abundance**: which taxa differ in relative abundance between samples
- **Step 6 - Community Clustering**: how similar/dissimilar the samples' microbial communities are (Bray-Curtis, PCoA, hierarchical clustering, PCA, heatmap)

This repo covers Steps 5-6 only. Steps 1-4 (assembly, gene prediction, binning,
read-based taxonomy) are handled in a separate pipeline.

## ⚠️ About the test dataset

This pipeline was validated on **three ONT shotgun samples pulled from two unrelated
public studies** - they are **not** a biologically matched comparison group, they exist
purely to test that the Step 5-6 code runs correctly end-to-end on real ONT data before
being pointed at the actual project dataset.

| Sample | Study | Source environment | BioProject |
|---|---|---|---|
| `SRR13573807` | Kinema (fermented soybean, Nepal) | Fermented food microbiome | [PRJNA695113](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA695113) |
| `SRR13573808` | Kinema (fermented soybean, Bhutan) | Fermented food microbiome | [PRJNA695113](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA695113) |
| `ERR15076229` | Schönbuch Forest soil metagenome (Germany) | Temperate forest soil | [PRJEB89893](https://www.ebi.ac.uk/ena/browser/view/PRJEB89893) |

In other words: 2 of the 3 samples are the *same* environment (fermented soybean) from
different countries, and the 3rd is a *completely different* environment (forest soil).
Any "differences" the pipeline reports between these three are expected to be large and
are a sanity check of the code, not a real biological finding - a genuine differential
abundance / clustering analysis needs samples that share a real experimental question
(e.g. multiple replicates per condition).

Results from this run will be added to this README once the pipeline has been executed
on all three samples.

## Workflow

```
Raw ONT FASTQ (SRR13573807, SRR13573808, ERR15076229)
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
  - `nanoplot_env` - NanoPlot
  - `filtlong_env` - Filtlong
  - `kraken2_env` - Kraken2, Bracken
  - `community_env`- pandas, numpy, matplotlib, scipy, scikit-learn
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
  sequencing errors - that only happens with amplicon sequencing (16S), where PCR
  produces millions of copies of one gene. Shotgun reads come from random genomic
  locations and rarely overlap, so denoising isn't applicable; clustering samples by
  abundance profile (Bray-Curtis / PCoA / hierarchical clustering / PCA) is the shotgun
  equivalent.
- `fastp`/Trimmomatic were dropped from the QC stage in favor of NanoPlot + Filtlong,
  which are built for ONT long reads rather than Illumina short reads.

## Results

# Biocrust ONT Metagenomics - QC & Taxonomic Profiling Summary

This section summarizes read-level QC (NanoPlot, pre/post filtering) and taxonomic
profiling (Kraken2/Bracken) for three Oxford Nanopore samples: `SRR13573807`,
`SRR13573808`, and `ERR15076229`.

---

## 1. Read QC - NanoPlot (Pre- vs Post-Filtering)

NanoPlot was run before and after quality/length filtering (Step 5–6) to check how
filtering changed read yield, length, and quality distribution for each sample.

### 1.1 Pre-filtering (raw reads)

| Metric | SRR13573807 | SRR13573808 | ERR15076229 |
|---|---|---|---|
| Mean read length (bp) | 721.2 | 543.4 | 6,469.3 |
| Mean read quality (Q) | 9.1 | 8.9 | 12.4 |
| Median read length (bp) | 391.0 | 356.0 | 4,347.0 |
| Median read quality (Q) | 9.2 | 9.0 | 13.5 |
| Number of reads | 708,244 | 1,022,385 | 6,006,176 |
| N50 (bp) | 1,211 | 757 | 11,317 |
| Total bases | 510.8 Mb | 555.6 Mb | 38.9 Gb |
| Reads >Q10 | 27.1% (137.6 Mb) | 21.3% (111.0 Mb) | 87.0% (33.8 Gb) |

### 1.2 Post-filtering (cleaned reads)

| Metric | SRR13573807 | SRR13573808 | ERR15076229 |
|---|---|---|---|
| Mean read length (bp) | 1,392.7 | 1,062.9 | 8,468.5 |
| Mean read quality (Q) | 9.1 | 8.9 | 13.3 |
| Median read length (bp) | 941.0 | 816.0 | 6,509.0 |
| Median read quality (Q) | 9.3 | 9.0 | 13.9 |
| Number of reads | 283,971 | 349,071 | 4,129,408 |
| N50 (bp) | 1,701 | 1,142 | 11,936 |
| Total bases | 395.5 Mb | 371.0 Mb | 35.0 Gb |
| Reads >Q10 | 26.7% (106.9 Mb) | 19.5% (71.4 Mb) | 95.6% (32.5 Gb) |

**Interpretation**

- **Filtering removed mostly short, low-quality reads.** Read counts dropped sharply
  (SRR13573807: 708k → 284k, SRR13573808: 1.02M → 349k, ERR15076229: 6.0M → 4.1M),
  while mean/median read length and N50 *increased* in every sample - this is the
  expected signature of a length/quality filter removing short junk reads rather
  than trimming good reads down.
- **ERR15076229 is a fundamentally higher-quality, longer-read run** than the two
  SRR samples: mean quality ~12–13 vs ~9, and read lengths roughly an order of
  magnitude longer (N50 ~11–12 kb vs ~1–1.7 kb). Its total yield (Gb-scale) also
  dwarfs the other two (Mb-scale), so it will dominate any pooled analysis unless
  normalized.
- **SRR13573807 and SRR13573808 remain low-quality overall** (mean Q ~9, i.e. ~87.5%
  base-call accuracy) even after filtering, and almost nothing passes Q15+ in either
  sample. This is consistent with older/basic ONT chemistry or a rapid/low-accuracy
  kit rather than R10.4-class flow cells.
- **%>Q10 barely moves after filtering for the SRR samples** (27.1%→26.7%,
  21.3%→19.5%), because the filter here is doing more length-based than
  quality-based selection — worth flagging if your filtering step intended to
  enforce a strict quality cutoff.

---

## 2. Taxonomic Classification — Kraken2 / Bracken

Kraken2 assigns each read to a taxon based on k-mer matches against a reference
database; Bracken then re-estimates species-level abundance from those Kraken2
assignments, correcting for the fact that reads at higher taxonomic ranks get
"pushed down" to species proportionally.

### 2.1 Top classified taxa across all samples (Kraken2 report)

| Taxon | ERR15076229 | SRR13573807 | SRR13573808 |
|---|---|---|---|
| *Homo sapiens* | 0.03384 | 0.00165 | 0.00156 |
| *Bradyrhizobium erythrophlei* | 0.02098 | 0.00000 | 0.00000 |
| *Reyranella* sp. | 0.00783 | 0.00000 | 0.00000 |
| *Bradyrhizobium diazoefficiens* | 0.00756 | 0.00000 | 0.00000 |
| *Tunturiibacter gelidiferens* | 0.00746 | 0.00000 | 0.00000 |
| *Bradyrhizobium* sp. Ash2021 | 0.00699 | 0.00000 | 0.00000 |
| *Salmonella enterica* | 0.00662 | 0.00000 | 0.00000 |
| *Klebsiella pneumoniae* | 0.00654 | 0.00000 | 0.00030 |
| *Pseudomonas aeruginosa* | 0.00632 | 0.00000 | 0.00000 |

*(values = fraction of total reads assigned to that taxon)*

**Interpretation**

- **The three samples have almost no taxonomic overlap.** Every top taxon in
  ERR15076229 is essentially absent (0.000) in both SRR samples - this is a strong
  signal that ERR15076229 is a genuinely different community (soil/biocrust
  environmental taxa: *Bradyrhizobium*, *Reyranella*, *Tunturiibacter* - all
  soil/root-associated bacteria) compared to the SRR pair, which is dominated by
  *Bacillus*-group organisms (see §2.2).
- ***Homo sapiens* reads are present in all three samples** (0.3-3.4%), which is
  expected host/handler contamination in environmental ONT sequencing and should be
  filtered out before downstream diversity analysis if not already removed.
- **ERR15076229 also carries some reads classified as clinical/enteric organisms**
  (*Salmonella enterica*, *Klebsiella pneumoniae*, *Pseudomonas aeruginosa*) at low
  abundance — plausible in soil/biocrust (these genera have free-living soil
  relatives and Kraken2 species calls at low read counts can be noisy), but worth a
  sanity check against your database's false-positive rate at low abundance.

### 2.2 Top species-level abundance - Bracken (per sample)

> Note: values below are the `fraction_total_reads` column from each sample's
> Bracken output (the corrected relative-abundance estimate), which is the
> field that matters for comparison — the raw `kraken_assigned_reads`/`added_reads`
> columns in the pasted terminal output are misaligned from copy/paste and
> shouldn't be read directly from the text above.

**SRR13573807**

| Species | Relative abundance |
|---|---|
| *Bacillus subtilis* | 27.0% |
| *Bacillus licheniformis* | 3.3% |
| *Bacillus velezensis* | 1.4% |
| *Bacillus paralicheniformis* | 0.5% |
| *Bacillus spizizenii* | 0.5% |
| *Bacillus tequilensis* | 0.3% |
| *Bacillus amyloliquefaciens* | 0.3% |
| *Bacillus inaquosorum* | 0.2% |

**SRR13573808**

| Species | Relative abundance |
|---|---|
| *Bacillus anthracis* | 6.6% |
| *Bacillus paranthracis* | 7.0% |
| *Bacillus thuringiensis* | 5.3% |
| *Bacillus cereus* | 5.2% |
| *Bacillus mycoides* | 0.5% |
| *Bacillus wiedmannii* | 0.4% |
| *Bacillus tropicus* | 0.4% |

**ERR15076229**

| Species | Relative abundance |
|---|---|
| *Bradyrhizobium* sp. Ash2021 | 0.70% |
| *Bradyrhizobium* sp. NP1 | 0.24% |
| *Bradyrhizobium* sp. STM 3562 | 0.17% |
| *Bradyrhizobium* sp. McL0615 | 0.16% |
| *Bradyrhizobium* genosp. P | 0.15% |
| *Bradyrhizobium* sp. 170 | 0.15% |

**Interpretation**

- SRR13573807 and SRR13573808 are both **overwhelmingly *Bacillus*-dominated**, but
  by *different* members of the group: SRR13573807 is led by *B. subtilis* + *B.
  licheniformis*, while SRR13573808 is led by the closely related *B. cereus
  group* (*B. anthracis*/*paranthracis*/*thuringiensis*/*cereus* - these four are
  genomically near-identical and routinely hard for k-mer classifiers to separate,
  so treat the *species* call cautiously; the *B. cereus group* as a whole is the
  real signal).
- ERR15076229 shows **no single dominant species** - abundances are spread thinly
  across dozens of *Bradyrhizobium* strains/species (each <1%), which is typical of
  a genuinely diverse soil/biocrust community rather than a culture-dominated or
  low-complexity sample.
- This split (SRR pair = *Bacillus*-dominated, ERR = diverse soil community)
  is the same pattern that shows up in the ordination in §3.

### 2.3 Relative abundance bar chart & heatmap (top 15 taxa)

The bar chart and heatmap visualize relative abundance of the top 15 taxa (union
across all samples) side by side per sample.

*Bar Chart Plot *
<img width="1115" height="597" alt="image" src="https://github.com/user-attachments/assets/530221e6-070e-485a-9c22-ae7dbe592728" />

*Heat Map Plot*
<img width="750" height="620" alt="image" src="https://github.com/user-attachments/assets/0aa79dc0-3610-466a-9d83-13e2a5045f99" />

**Key patterns:**
- *Caldifermentibacillus hisashii* and *Bacillus subtilis* are the two largest bars
  overall, but almost entirely restricted to SRR13573807 (~42% and ~27%
  respectively) and, to a lesser extent, SRR13573808.
- *Bacillus paranthracis*, *B. anthracis*, *B. thuringiensis*, *B. cereus*, and
  *Proteus mirabilis* form a cluster of moderate-abundance taxa (~5–7% each) that
  are essentially unique to SRR13573808.
- *Escherichia coli* and *Homo sapiens* are the only taxa in this top-15 set with
  meaningful abundance in ERR15076229 (~6% and ~3%) - everything else in this
  particular taxon list is near-zero for that sample, which is why the heatmap's
  ERR15076229 column is almost entirely dark purple (low values) except those two
  rows.
- The heatmap makes the same point visually as the bar chart but is easier to scan
  for "which sample owns which taxon" - SRR13573807 and SRR13573808 share the
  yellow/green (high-abundance) cells for *B. subtilis*-group organisms, while
  ERR15076229 has no bright cells in this taxon list at all (its real dominant
  taxa - the *Bradyrhizobium* strains - aren't part of this particular top-15
  set, since it's a top-15 union chart, not per-sample top-15).

---

## 3. Beta Diversity — Bray-Curtis PCoA & Hierarchical Clustering

Both plots are built from the same **Bray-Curtis dissimilarity matrix** computed
on the full species-abundance table (Bracken output), so they should - and do -
tell a consistent story.

### 3.1 PCoA (Principal Coordinates Analysis)

*PCoA Plot*
<img width="700" height="499" alt="image" src="https://github.com/user-attachments/assets/7248d6e9-6551-48ce-9223-4c767d747a21" />

- **PCoA1 explains 78.1%** of the variance, **PCoA2 explains 21.9%** — together
  these two axes capture essentially all the compositional variation between the
  three samples (unsurprising with only 3 samples/2 possible axes).
- **ERR15076229 sits far to the left, alone**, separated from the other two samples
  almost entirely along PCoA1. This is the dominant axis of variation and reflects
  the *Bacillus*-dominated vs. *Bradyrhizobium/soil*-dominated community split
  described in §2.2.
- **SRR13573807 and SRR13573808 sit close together on PCoA1** (both near the
  right side) but are separated from each other along **PCoA2**, reflecting their
  different *Bacillus* sub-communities (*B. subtilis*-group vs. *B. cereus*-group).

### 3.2 Hierarchical clustering (on the same Bray-Curtis distances)

*Hierarchical Clustering Plot*
<img width="955" height="552" alt="image" src="https://github.com/user-attachments/assets/cfb1605c-6726-457c-a30b-539abae42035" />

- **ERR15076229 branches off first**, at a distance of ~0.98 (i.e. it is almost
  maximally dissimilar - Bray-Curtis ranges 0-1 - from either SRR sample).
- **SRR13573807 and SRR13573808 join at a much lower distance (~0.58)**, confirming
  they are compositionally more similar to *each other* than either is to
  ERR15076229, even though §2.2 shows they're led by different *Bacillus* species.

**Bottom line:** all three analyses (Bracken abundance, PCoA, dendrogram) agree -
ERR15076229 is compositionally distinct from the SRR13573807/SRR13573808 pair,
which likely reflects either a different sample type/site, a different sequencing
run/library prep, or a genuinely different microbial community (soil/biocrust
diversity vs. a *Bacillus*-dominated community, possibly enrichment or
lower-diversity substrate). 

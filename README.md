# Metagenomic Population Tracking with ANIr  
*A KDE-based framework for tracking sequence-discrete microbial populations across metagenomes*

This repository implements the full analysis workflow described in:

**Conrad, R. et al.**  
*An ANIr-based methodology to determine whether two sequence-discrete populations are identical or distinct in metagenomic datasets.*  
(in review)

Zenodo archive (code + data snapshot):  
https://doi.org/10.5281/zenodo.15231818

---

## Overview

This repository provides a complete, end-to-end computational workflow for **tracking microbial populations across metagenomic datasets using Average Nucleotide Identity of recruited reads (ANIr)**. The framework is designed to determine whether two metagenomic samples contain:

- The **same population** (shared strain-level lineage), or  
- **Distinct but closely related populations**, even when 16S identity or species-level resolution fails.

The pipeline is built around three core concepts:

1. **Recruitment-based population definition**  
   Metagenomic reads are recruited to reference genomes or MAGs using BLAST+ or MagicBlast. The distribution of percent identity values (`pident`) is used to separate target populations from non-target background recruitment.

2. **KDE-based threshold detection**  
   A **kernel density estimate (KDE)** is applied to the `pident` distribution, and **local minima** in the density curve are used as objective sequence-identity cutoffs separating true population members from non-members.

3. **Statistical population comparison**  
   Once target populations are isolated, the framework applies:
   - **ANIr and mANIr** (mean and median percent identity)
   - **Bootstrap resampling** to estimate confidence intervals
   - **Permutation testing** to assess whether two samples contain statistically indistinguishable populations
   - **Differential gene-content testing** using normalized per-gene coverage

Together, these tools allow you to answer:

> *Are two metagenomes recruiting the same microbial population, or two distinct ones?*

---

## Conceptual Workflow

At a high level, the population-tracking workflow follows these stages:

1. **Input data preparation**
   - Reference genomes or MAGs (FASTA)
   - Gene annotations (e.g., Prodigal CDS)
   - Metagenomic read sets (FASTQ)
   - Read recruitment tables (BLAST+ or MagicBlast)

2. **KDE-based population boundary detection**
   - KDE applied to `pident` distributions
   - Local minima used to define recruitment threshold
   - Target vs non-target reads separated objectively

3. **Population coverage & diversity metrics**
   - Breadth and depth of coverage
   - ANIr and mANIr calculated from thresholded reads

4. **Across-sample population matrices**
   - Presence/absence across samples
   - ANIr matrices across environments
   - Heatmaps and clustering

5. **Statistical population comparison**
   - Bootstrap testing of ANIr distributions
   - Permutation testing of ANIr differences
   - Identifies identical vs distinct populations

6. **Differential gene content**
   - Per-gene normalized coverage
   - Quantile-based detection of enriched/depleted genes
   - Functional differentiation between populations

---

## What This Repository Contains

This repository is organized into three main script families:

- **`scripts-gene-coverage/`**  
  Raw recruitment → coverage, ANIr, and per-gene depth calculations.

- **`scripts-populations/`**  
  KDE thresholding, population-level statistics, bootstrap/permutation testing, and differential gene analysis.  
  **This is the primary focus of the population-tracking workflow.**

- **`scripts-general/`**  
  Utility scripts, plotting helpers, and data-wrangling functions.

A legacy `old/` directory is preserved for transparency but is **not recommended for new analyses**.

A full breakdown of all major scripts is provided in later sections of this README.

---

## Software Requirements

### Core Environment
- Python ≥ 3.7.10  
- NumPy  
- Pandas  
- Matplotlib  
- Seaborn  
- SciPy  
- Statsmodels

### External Tools
- **BLAST+** (for `blastn` recruitment)
- **MagicBlast** (optional alternative recruiter)
- **Prodigal** or equivalent CDS caller (for gene-level analyses)

Exact versions used in the manuscript are documented in the Zenodo archive.

---

## Quickstart (High-Level)

This is the minimal outline of a complete run:

```bash
# 1. Recruit reads to MAGs or references
blastn -query reads.fq -db MAGs.fna -outfmt 6 > sampleA.MAG001.blast.tsv

# 2. Run KDE + recruitment plotting + ANIr estimation
python scripts-populations/06f_TabBlast_RecPlot_Mini_Auto_vX.py \
  --blast sampleA.MAG001.blast.tsv \
  --fasta MAG001.fna \
  --outdir results/MAG001/sampleA/

# 3. Build across-sample ANIr matrices

# 4. Run bootstrap and permutation tests between samples

# 5. Run differential gene-content analysis
```

A **fully reproducible, step-by-step tutorial for every script in `scripts-populations/` is provided in the next major section of this README.**

---

## How to Cite

If you use this repository in your work, please cite:

* The associated manuscript (once published)
* The Zenodo software archive:
  [https://doi.org/10.5281/zenodo.15231818](https://doi.org/10.5281/zenodo.15231818)

---

## Repository Structure

This repository is organized into three main script families plus supporting files:

```text
Metagenomic_Population_Tracking/
├── scripts-gene-coverage/     # Read recruitment → coverage, ANIr, RecPlot minis, gene-level depth
├── scripts-populations/       # Population thresholding, ANIr comparisons, bootstrap/permutation tests
├── scripts-general/           # Shared utilities, plotting helpers, data-wrangling functions
├── old/                       # Legacy scripts kept for transparency (not recommended for new use)
├── LICENSE                    # GPL-3.0
└── README.md                  # This document
````

In practice:

* Use **`scripts-gene-coverage/`** when you’re starting from **raw recruitment tables** and want coverage / ANIr / RecPlot minis / gene depth.
* Use **`scripts-populations/`** when you’re doing **population-level comparisons** (same vs different population across samples, differential gene content).
* Use **`scripts-general/`** for shared helper functions that some of the above scripts call.
* The **`old/`** directory contains earlier versions and exploratory scripts and is mostly for archival / provenance.

---

## Script Families

To make the workflow easier to follow, it helps to think in terms of functional “blocks” rather than individual scripts first:

1. **Recruitment & KDE / RecPlot minis**

   * Parse BLAST+/MagicBlast recruitment output.
   * Fit KDE to `%identity` (seqIDs) and identify local minima.
   * Apply a valley modifier and optional KDE bandwidth tweaks.
   * Compute breadth, depth, ANIr, and mANIr.
   * Generate **RecPlot minis** that overlay coverage and identity patterns.

2. **Population presence / absence and matrices**

   * Aggregate per-MAG×sample metrics into matrices (ANIr, breadth, depth).
   * Build presence/absence calls and sample × population tables.
   * Provide inputs for downstream clustering and metadata analyses.

3. **Population comparison tests (ANIr-based)**

   * Use **bootstrap** to generate ANIr distributions per sample.
   * Use **permutation tests** to evaluate significance of ANIr differences.
   * Return p-values and confidence intervals for “same vs different population”.

4. **Differential gene content**

   * Compute per-gene depth using positional information (`sstart`, `send`, gene coordinates).
   * Normalize gene coverage by genome-level abundance.
   * Identify genes in the extreme quantiles (e.g., bottom 2.5% / top 97.5%) of abundance differences between samples.

5. **Metadata integration & plotting**

   * Correlate ANIr and population presence/absence with environmental or experimental metadata.
   * Produce publication-ready figures (e.g., ANIr vs depth, basin, time).

---

## Key Scripts at a Glance

> 💡 Script names below correspond to those referenced in the manuscript.
> Depending on how you organize the repo, they may live under `scripts-gene-coverage/` or `scripts-populations/`.

| Script / family                             | Main purpose                                                   | Typical inputs                                                                             | Key outputs                                                                                                           |
| ------------------------------------------- | -------------------------------------------------------------- | ------------------------------------------------------------------------------------------ | --------------------------------------------------------------------------------------------------------------------- |
| `TabBlast_RecPlot_Mini_Auto_v2.py / v3.py`  | KDE-based thresholding, RecPlot minis, coverage & ANIr metrics | BLAST+ tabular recruitment file, reference genome/MAG FASTA                                | KDE plots, local-minimum threshold, RecPlot minis, per-MAG×sample table of breadth, depth, ANIr, mANIr, read counts   |
| `MagicBlast_RecPlotMini-Auto_subspecies.py` | Same as above, tuned for MagicBlast output                     | MagicBlast tabular recruitment file, reference genome/MAG FASTA                            | KDE + RecPlot minis + coverage / ANIr summary as above                                                                |
| `06h_PopX_Bootstrap_Test.py`                | Bootstrap test for ANIr distributions between two samples      | SeqID arrays (or recruitment tables) for sample A and B mapped to the same MAG / reference | Bootstrap ANIr distributions, confidence intervals, p-values for “same population?”                                   |
| `06h_PopX_Permutation_Test.py`              | Permutation test for ANIr differences between two samples      | Same as above (seqIDs for two samples for a single MAG)                                    | Null distribution of ANIr differences, observed difference, p-values                                                  |
| `08b_BlastPlus_CoverageMagic_Basic.py`      | Per-gene coverage and normalized gene abundance                | BLAST+ recruitment tables, gene (CDS) coordinates on contigs, genome abundance estimates   | Per-gene depth / coverage tables, gene-level normalized abundance per sample                                          |
| `08d_MagicBlast_Differential_Genes.py`      | Differential gene content between two samples                  | Per-gene coverage / normalized abundance tables for sample pairs                           | Gene list with difference metrics; genes flagged in bottom 2.5% / top 97.5% quantiles (significant abundance changes) |
| RecPlot aggregation scripts (e.g. “minis”)  | Arrange many RecPlot minis into composite figures              | Individual RecPlot outputs, sample metadata (e.g., basin, depth)                           | Multi-panel PDF/PNG grids of RecPlot minis arranged by geography / depth / other metadata                             |
| ANIr / metadata plotting scripts            | Relate ANIr and population metrics to environmental metadata   | ANIr matrices, presence/absence tables, sample metadata (depth, basin, time, etc.)         | Scatterplots, timeseries, and faceted plots linking population metrics to environment                                 |

You can think of the execution path like this:

1. **Start with** `TabBlast_RecPlot_Mini_Auto_*` / `MagicBlast_RecPlotMini-*` to define populations and compute coverage + ANIr.
2. **Aggregate** per-MAG×sample results into multi-sample matrices.
3. **Apply** `06h_PopX_*` scripts to test “same vs different” populations between sample pairs.
4. **Optionally run** `08b_*` and `08d_*` for differential gene content within those same sample pairs.
5. **Finish with** plotting/metadata scripts to put everything in ecological / experimental context.

---

## Software & Dependencies

These scripts were developed and tested under:

* **Python**: 3.7.10+
* **Core libraries**:

  * `numpy`
  * `pandas`
  * `matplotlib`
  * `seaborn`
  * `scipy`
  * `statsmodels`

**External tools:**

* **BLAST+** (`blastn`) – read recruitment against MAGs / reference genomes.
* **MagicBlast** – optional alternative recruiter.
* **Prodigal** (or equivalent) – gene prediction (CDS coordinates) for per-gene coverage / differential gene tests.
* Common command-line tools (`bash`, `gzip`, etc.) for file management and pre-processing.

A simple `conda` environment example:

```bash
conda create -n anir-populations python=3.8 \
  numpy pandas matplotlib seaborn scipy statsmodels
conda activate anir-populations

# Install BLAST+ and (optionally) MagicBlast via bioconda
conda install -c bioconda blast
conda install -c bioconda magicblast   # optional
conda install -c bioconda prodigal     # for gene prediction
```

---

## 3. Population-level analysis (`scripts-populations/`)

This folder contains the population-level analysis tools used in the manuscript, including:

1. **Population detection and RecPlot mini-figures** (ANIr, breadth, depth)

   * `06f_TabBlast_RecPlot_Mini_Auto_v4.py` 
   * `MagicBlast_RecPlotMini-Auto_subspecies.py` 

2. **Population ANIr significance tests**

   * `06h_PopX_Bootstrap_Test.py` 
   * `06h_PopX_Permutation_Test.py` 

3. **Gene-level coverage and differential gene content**

   * `08c_Low_Coverage_Genes.py` 
   * `08d_Differential_Genes.py` 

Below is a step-by-step guide to using these scripts, in the typical order you would run them.

---

### 3.1. Input expectations and directory layout

Most scripts assume:

* **One reference genome / rMAG FASTA** per run (can be multi-contig).
* **Tabular BLAST+ or MagicBlast read recruitment files**:

  * One file per **sample × reference genome**.
  * “Outfmt 6”-style with at least:

    * `% identity` in column 3 (`pident`)
    * subject start/stop in columns 9–10 (`sstart`, `send`) – 0-based vs 1-based is handled consistently in the scripts.
* **Consistent naming convention for sample IDs** encoded in the filenames.

For the GoM case, filenames look like:

```text
GoM_1200m_2016-01-01_MAG001.tab
GoM_80m_2016-01-01_MAG001.tab
...
```

Internally, the scripts parse:

* `Sample` = everything before the first `.` and before the first `-`
* `depth` = the third `_`-separated field (e.g. `1200m`)
* `loca` = the first two `_`-separated fields (e.g. `GoM_1200m`) 

If your sample naming scheme is different, you **must** edit the lines that set `Sample`, `depth`, and `loca` in the plotting scripts.

---

### 3.2. Step 1 – Detect target populations and build RecPlot minis

#### 3.2.1. `06f_TabBlast_RecPlot_Mini_Auto_v4.py`

**Purpose**

* Reads a **reference FASTA** and a folder of **filtered BLAST tabular files** (one per sample for that reference).
* Identifies a **population cutoff** in the seqID distribution via KDE and valley detection.
* Filters alignments above the population cutoff.
* Calculates:

  * **Breadth of coverage** (fraction of genome positions with ≥1 population read).
  * **Depth of coverage** (truncated mean coverage).
  * **ANIr** and **MNIr** (mean and median percent identity of the population reads).
* Produces:

  * Per-sample **RecPlot mini** PDFs (pident vs genome position).
  * A summary `*_stats.tsv` file.
  * A `*_popx_readmap_values.tsv` file used downstream by the bootstrap and permutation tests. 

**Key inputs**

* `-f / --input_fasta_file`
  Reference genome (or rMAG) FASTA used for read recruitment.
* `-b / --input_blast_dir`
  Directory containing tabular BLAST result files (one per sample).
* `-o / --out_dir`
  Output directory (created if it does not exist).
* `-k / --kde_bandwidth` *(optional, float, default ~0.3)*
  KDE bandwidth for smoothing seqID distributions.
* `-v / --valley_modifier` *(optional, float, default 3.0)*
  Value added to the KDE valley to define the population cutoff (in % identity).
* `-d / --draw_threshold` *(flag)*
  If present, draw the horizontal population cutoff line on RecPlots.
* `-t / --breadth_threshold` *(optional, float, default 10.0)*
  Minimum breadth (%) required to call a population “detected”. 

**Example command**

```bash
# Example: one rMAG, BLAST-based read recruitment files per sample
python 06f_TabBlast_RecPlot_Mini_Auto_v4.py \
    -f refs/MAG001.fasta \
    -b blasts/MAG001/ \
    -o results/MAG001_popx/ \
    -k 0.25 \
    -v 3.0 \
    -d \
    -t 10.0
```

**Outputs**

From the example above, you’ll get (under `results/MAG001_popx/`):

* `MAG001_stats.tsv`

  * Columns: `Sample`, `Detected` (0/1), `Depth`, `Breadth`, `ANIr`, `MNIr`, `PopCutoff`. 
* `MAG001_popx_readmap_values.tsv`

  * Contains per-sample **population-level metrics** and the raw seqID array for the population reads:
  * Columns (for each sample):
    `Sample`, `PopCutoff`, `ExactMatchRatio(EMR)`, `BreadthTotal`, `DepthTotal`,
    `Breadth100`, `Depth100`, `Breadth99`, `Depth99`,
    `pident above population threshold` (comma-separated list of seqIDs). 
  * This file is the **input** for `06h_PopX_Bootstrap_Test.py` and `06h_PopX_Permutation_Test.py`.
* Per-sample plots:

  * `SampleX_RecPlot_Mini.pdf`

    * KDE-based population cutoff (optional line).
    * SeqID vs genome position; ANIr, MNIr, depth, and breadth annotated in the plot. 
  * `SampleX_Pop_Cutoff.png`

    * Histogram and KDE of seqIDs with the chosen cutoff line.

---

#### 3.2.2. `MagicBlast_RecPlotMini-Auto_subspecies.py`

This script is the **MagicBlast-friendly** counterpart of `06f_TabBlast_RecPlot_Mini_Auto_v4.py`. It:

* Expects MagicBlast-style tabular recruitment instead of BLAST+.
* Performs the **same** KDE valley-based population cutoff search, breadth/depth/ANIr/MNIr calculations, and RecPlot mini generation. 

The CLI interface mirrors `06f_TabBlast_RecPlot_Mini_Auto_v4.py`:

* `-f / --input_fasta_file`
* `-b / --input_blast_dir`
* `-o / --out_dir`
* `-k / --kde_bandwidth`
* `-v / --valley_modifier`
* `-d / --draw_threshold`
* `-t / --breadth_threshold`

You can inspect `python MagicBlast_RecPlotMini-Auto_subspecies.py -h` for the exact help text, but usage is effectively identical: point it at a reference FASTA and a directory of MagicBlast tabular results, and it will produce the same types of outputs (including the `*_popx_readmap_values.tsv` file used by downstream tests). 

---

### 3.3. Step 2 – Bootstrap tests of ANIr (population stability)

#### `06h_PopX_Bootstrap_Test.py`

**Purpose**

Test the **stability of ANIr** for each detected population by:

* Bootstrapping subsamples of the seqID array (population reads) from `*_popx_readmap_values.tsv`.
* Building a **bootstrap distribution of ANIr** for each sample.
* Estimating how extreme the observed ANIr is relative to the bootstrap distribution. 

The underlying method (as described in the manuscript):

* For each sample, repeatedly sample a fraction of the seqID array (`bootsize`, default 2% of reads) with replacement.
* Compute the mean seqID (ANIr) for each bootstrap.
* From the bootstrap distribution, compute:

  * mean, standard deviation, min, max.
  * probability of observing the measured ANIr (p-value) under this distribution. 

**Required input**

* `*_popx_readmap_values.tsv` output from `06f_TabBlast_RecPlot_Mini_Auto_v4.py` or `MagicBlast_RecPlotMini-Auto_subspecies.py`. 

**Arguments**

* `-i / --input_popx_file`
  PopX file (e.g., `MAG001_popx_readmap_values.tsv`).
* `-t / --target_sample`
  Name of the **target sample** (must match a `Sample` in the input file).
* `-o / --output_prefix`
  Prefix for all output files (e.g., `MAG001_bootstrap`).
* `-b / --bootstrap_boot_size` *(optional, default 0.02)*
  Fraction of the sample’s reads per bootstrap. `0.02` ≈ 2% of the reads. 
* `-n / --bootstrap_number_boots` *(optional, default 10000)*
  Number of bootstrap replicates. 

**Example command**

```bash
python 06h_PopX_Bootstrap_Test.py \
    -i results/MAG001_popx/MAG001_popx_readmap_values.tsv \
    -t GoM_1200m_2016-01-01 \
    -o results/MAG001_popx/MAG001_bootstrap \
    -b 0.02 \
    -n 10000
```

**Outputs**

* `MAG001_bootstrap_<Sample>_PopX_Bootstrap.png` (per sample)

  * Histogram of bootstrapped ANIr values with fitted normal PDF.
  * Annotated with mean, 3×SD, measured ANIr, and p-value relative to the bootstrap distribution. 
* `MAG001_bootstrap_bootstrap_data.tsv`

  * Columns (per sample):
    `Sample`, `Boot_Mean`, `Boot_Stdev`, `Boot_Min`, `Boot_Max`, `Boot_pvalue`.
  * A placeholder row for the target sample is included with `"Target"` entries. 

---

### 3.4. Step 3 – Permutation tests of ANIr differences between samples

#### `06h_PopX_Permutation_Test.py`

**Purpose**

Test whether the **ANIr difference between two samples** (for the same reference genome) is larger than expected if both samples came from the same underlying population.

* Implements a **permutation test** on seqIDs (pidents) drawn from `*_popx_readmap_values.tsv`. 

Conceptually (also described in the manuscript):

* Choose a **target sample** A.
* For each other sample B:

  * Combine the seqID arrays from A and B.
  * Randomly shuffle (permute) the combined array.
  * Split into two new arrays of sizes matching A and B.
  * Compute ANIr for each array and record their difference.
  * Repeat many times (default 10,000) to obtain a null distribution of ANIr differences under the hypothesis A and B are the same population.
  * Compare the **observed** |ANIr(A) − ANIr(B)| to this null distribution and compute a p-value. 

**Required input**

* Same `*_popx_readmap_values.tsv` file used by the bootstrap script. 

**Arguments**

* `-i / --input_popx_file`
  PopX file (e.g., `MAG001_popx_readmap_values.tsv`).
* `-t / --target_sample`
  Target sample name (the reference population).
* `-o / --output_prefix`
  Prefix for all outputs (e.g., `MAG001_permutation`).
* `-p / --number_of_permutations` *(optional, default 10000)*
  Number of permutations to run. 

**Example command**

```bash
python 06h_PopX_Permutation_Test.py \
    -i results/MAG001_popx/MAG001_popx_readmap_values.tsv \
    -t GoM_1200m_2016-01-01 \
    -o results/MAG001_popx/MAG001_perm \
    -p 10000
```

**Outputs**

* `MAG001_perm_permutation_data.tsv`

  * Per sample (relative to target sample):
    `Sample`, `ANIr`, `MNIr`, `EMR`, `BreadthTotal`, `DepthTotal`, `Breadth100`, `Depth100`, `Breadth99`, `Depth99`, `sample_size`,
    `Measured_Diff`, `Perm_Mean_Diff`, `Perm_diff_Stdev`, `pvalue_Measured_Diff`. 
* `MAG001_perm_<Sample>_permutation_test.png` (one per non-target sample)

  * Histogram of simulated ANIr differences from permutations.
  * Fitted normal PDF.
  * Lines for:

    * 3×SD of simulated differences.
    * Observed ANIr difference.
  * Annotated p-value for the observed difference. 

---

### 3.5. Step 4 – Gene-level coverage: low-coverage genes within a single sample

#### `08c_Low_Coverage_Genes.py`

**Purpose**

Identify **significantly low coverage genes** within a sample, based on the distribution of normalized gene coverage values (`TAD`, typically from `*_gene_tad.tsv` produced by the coverage-magic script).

* Uses a **non-parametric approach**:

  * Computes the **quantile** (default 0.025) of the gene coverage distribution.
  * Flags all genes with coverage ≤ that quantile as **significantly low** (α = 0.025). 

**Inputs**

* `*_gene_tad.tsv` file (one per sample) from coverage magic, with at least:

  * `Gene` ID.
  * `TAD` (normalized coverage depth per gene). 
* Optional:

  * `genes.faa` FASTA containing gene sequences.
  * `*.annotations` file from MicrobeAnnotator (tabular, gene as first column). 

**Arguments**

* `-i / --input_gene_tad_file`
  Input gene coverage file (e.g., `MAG001_GoM_1200m_gene_tad.tsv`).
* `-o / --output_prefix`
  Prefix for output files.
* `-s / --significance_alpha` *(optional, default 0.025)*
  Quantile used as the low-coverage cutoff.
* `-a / --annotation_file` *(optional)*
  MicrobeAnnotator annotations file.
* `-f / --fasta_file` *(optional)*
  Gene FASTA file. 

**Example command**

```bash
python 08c_Low_Coverage_Genes.py \
    -i coverage_magic/MAG001_GoM_1200m_gene_tad.tsv \
    -o results/MAG001_lowcov/MAG001_GoM_1200m \
    -s 0.025 \
    -a annotations/MAG001.annotations \
    -f genes/MAG001_genes.faa
```

**Outputs**

* `MAG001_GoM_1200m_gene_coverage.png`

  * Histogram of gene coverage values.
  * Vertical line at the low-coverage cutoff q (e.g., 2.5th percentile), annotated with:

    * α (significance level),
    * cutoff value,
    * # of genes ≤ cutoff,
    * # of genes with coverage = 0. 
* `MAG001_GoM_1200m_sig_genes.tsv`

  * Significantly low coverage genes (≤ cutoff).
  * Columns: `Gene`, `TAD`, and optional annotation fields. 
* `MAG001_GoM_1200m.fasta` *(optional)*

  * FASTA of significantly low coverage genes if `--fasta_file` was provided. 

---

### 3.6. Step 5 – Differential gene coverage between samples

#### `08d_Differential_Genes.py`

**Purpose**

Identify genes whose **relative coverage differs significantly between samples**, for a given reference genome.

As described in the manuscript:

* For each gene and each pair of samples, compute the **difference in normalized gene coverage**.
* Use the distribution of differences across all genes to find:

  * genes in the **bottom 2.5%** (underenriched),
  * genes in the **top 2.5%** (overenriched),
    corresponding to a two-sided α = 0.05. 

**Inputs**

* A directory containing `*_gene_tad.tsv` files from coverage magic:

  * One file per sample for the same reference genome; they must share a consistent basename suffix (e.g., `*_gene_tad.tsv`). 
* Optional:

  * Annotations and gene FASTA, as above.

**Arguments**

* `-i / --input_directory`
  Directory with the `*_gene_tad.tsv` files for this reference. 
* `-t / --target_sample`
  The sample name to treat as the **reference** in pairwise comparisons (must match the prefix used in the gene TAD filenames).
* `-o / --output_directory`
  Output directory to store results.
* `-s / --significance_alpha` *(optional, default 0.025)*
  Lower quantile for significance; the script uses `α` and `1 − α` to define low and high cutoff percentiles (e.g., 0.025 and 0.975 for a 5% two-sided test). 
* `-a / --annotation_directory` *(optional)*
  Annotations file from MicrobeAnnotator.
* `-f / --fasta_directory` *(optional)*
  Gene FASTA file. 

**Example command**

```bash
python 08d_Differential_Genes.py \
    -i coverage_magic/MAG001/ \
    -t GoM_1200m_2016-01-01 \
    -o results/MAG001_diffgenes/ \
    -s 0.025 \
    -a annotations/MAG001.annotations \
    -f genes/MAG001_genes.faa
```

**What the script does (high level)**

1. **Parse all `*_gene_tad.tsv` files** in `input_directory`.

   * For each sample:

     * Reads gene IDs.
     * Reads **normalized TAD** values and raw TAD. 
   * Stores:

     * log-transformed normalized TADs (`lgnmtd`) for statistical comparisons,
     * normalized TADs (`nmtd`) and raw TADs (`td`) for reporting. 

2. **Select the target sample** and compute per-gene differences:

   * Difference = log(norm TAD_target) − log(norm TAD_sample). 
   * For each non-target sample, compute:

     * the low quantile (`lq`) at `α`,
     * the high quantile (`hq`) at `1 − α`,
     * counts of genes below `lq` and above `hq`. 

3. **Plot differential coverage distributions**:

   * For each sample vs target:

     * Histogram of per-gene differential log(TAD/AVG).
     * Vertical lines for `lq` and `hq`.
     * Text annotation summarizing α, cutoffs, and counts of significant genes. 

4. **Write out results**:

   * For each sample vs target, four main outputs:

     * `<target>-<sample>_diffcov_plot.png`
       Histogram of log(TAD/AVG) differences with quantile lines. 
     * `<target>-<sample>_diffcov_high.tsv`
       Genes with differential coverage ≥ `hq`. Columns:
       `Gene`, `log(TAD/AVG) Diff`, `Target_Norm_TAD`, `Test_Norm_TAD`, `Target_TAD`, `Test_TAD`, and optional annotations. 
     * `<target>-<sample>_diffcov_low.tsv`
       Genes with differential coverage ≤ `lq`, same columns as above. 
     * Optional FASTA files:
       `<target>-<sample>_diffcov_high.fa` and `<target>-<sample>_diffcov_low.fa` containing coding sequences for high and low differential genes if a gene FASTA was provided. 

---

### 3.7. Typical end-to-end usage pattern

For a single reference genome / rMAG:

1. **RecPlot + population detection**

   * Use BLAST+ or MagicBlast to recruit reads.
   * Run `06f_TabBlast_RecPlot_Mini_Auto_v4.py` or `MagicBlast_RecPlotMini-Auto_subspecies.py` on each reference genome to:

     * generate RecPlot minis per sample,
     * compute breadth/depth/ANIr/MNIr,
     * produce `*_popx_readmap_values.tsv`. 

2. **Population ANIr statistics**

   * Run `06h_PopX_Bootstrap_Test.py` to characterize **within-sample** ANIr uncertainty. 
   * Run `06h_PopX_Permutation_Test.py` to test **between-sample** ANIr differences for the target sample against all others. 

3. **Gene-level coverage and content differences**

   * Use coverage-magic (upstream script) to generate per-sample `*_gene_tad.tsv` files.
   * Optionally run `08c_Low_Coverage_Genes.py` on each sample to identify significantly undercovered genes. 
   * Run `08d_Differential_Genes.py` on the folder of gene TAD files to identify genes that differ significantly in relative coverage between the target sample and each other sample. 

---

### 3.8. Worked Example: One rMAG across multiple metagenomes

This worked example demonstrates a **complete population-tracking analysis** for a single reference genome (rMAG) across multiple metagenomic samples, starting from BLAST recruitment and ending with differential gene content.

We assume the following directory layout:

```text
project/
├── refs/
│   └── MAG001.fasta
├── blasts/
│   └── MAG001/
│       ├── GoM_1200m_2016-01-01_MAG001.tab
│       ├── GoM_80m_2016-01-01_MAG001.tab
│       └── GoM_10m_2016-01-01_MAG001.tab
├── coverage_magic/
│   └── MAG001/
│       ├── MAG001_GoM_1200m_gene_tad.tsv
│       ├── MAG001_GoM_80m_gene_tad.tsv
│       └── MAG001_GoM_10m_gene_tad.tsv
├── annotations/
│   └── MAG001.annotations
├── genes/
│   └── MAG001_genes.faa
└── results/
````

---

## Step 1 — Detect target populations and generate RecPlot minis

We first identify the target population in each metagenomic sample using KDE-based valley detection.

```bash
python scripts-populations/06f_TabBlast_RecPlot_Mini_Auto_v4.py \
  -f refs/MAG001.fasta \
  -b blasts/MAG001/ \
  -o results/MAG001_popx/ \
  -k 0.25 \
  -v 3.0 \
  -d \
  -t 10.0
```

### Outputs generated:

```text
results/MAG001_popx/
├── MAG001_stats.tsv
├── MAG001_popx_readmap_values.tsv
├── GoM_1200m_2016-01-01_RecPlot_Mini.pdf
├── GoM_80m_2016-01-01_RecPlot_Mini.pdf
├── GoM_10m_2016-01-01_RecPlot_Mini.pdf
├── GoM_1200m_2016-01-01_Pop_Cutoff.png
├── GoM_80m_2016-01-01_Pop_Cutoff.png
└── GoM_10m_2016-01-01_Pop_Cutoff.png
```

At this stage you now have:

* **Population cutoff** for each sample
* **ANIr, MNIr, depth, and breadth**
* **Publication-ready RecPlot minis**
* A master file used for all downstream population tests:

```text
MAG001_popx_readmap_values.tsv
```

---

## Step 2 — Bootstrap the target population ANIr distribution

We now quantify **within-sample ANIr uncertainty** for one focal population (e.g., the deep GoM sample):

```bash
python scripts-populations/06h_PopX_Bootstrap_Test.py \
  -i results/MAG001_popx/MAG001_popx_readmap_values.tsv \
  -t GoM_1200m_2016-01-01 \
  -o results/MAG001_popx/MAG001_bootstrap \
  -b 0.02 \
  -n 10000
```

### Outputs:

```text
results/MAG001_popx/
├── MAG001_bootstrap_bootstrap_data.tsv
└── MAG001_bootstrap_GoM_1200m_2016-01-01_PopX_Bootstrap.png
```

This gives you:

* The **bootstrapped ANIr distribution**
* The **mean, SD, and p-value** of the observed ANIr
* A **PDF+histogram figure** used to evaluate population stability

---

## Step 3 — Permutation tests between populations

We now test whether **other samples represent the same population** as the 1200 m target sample.

```bash
python scripts-populations/06h_PopX_Permutation_Test.py \
  -i results/MAG001_popx/MAG001_popx_readmap_values.tsv \
  -t GoM_1200m_2016-01-01 \
  -o results/MAG001_popx/MAG001_perm \
  -p 10000
```

### Outputs:

```text
results/MAG001_popx/
├── MAG001_perm_permutation_data.tsv
├── MAG001_perm_GoM_80m_2016-01-01_permutation_test.png
└── MAG001_perm_GoM_10m_2016-01-01_permutation_test.png
```

This tells you, for each comparison sample:

* Observed **ANIr difference**
* Mean and SD of the **null permutation distribution**
* The **p-value** testing if the populations are distinguishable

At this point you can formally say:

> “These two samples contain the same population”
> or
> “These samples contain statistically distinct populations.”

---

## Step 4 — Identify low-coverage genes within a single sample (optional)

If you want to identify genes that are significantly underrepresented in a single sample:

```bash
python scripts-populations/08c_Low_Coverage_Genes.py \
  -i coverage_magic/MAG001/MAG001_GoM_1200m_gene_tad.tsv \
  -o results/MAG001_lowcov/MAG001_GoM_1200m \
  -s 0.025 \
  -a annotations/MAG001.annotations \
  -f genes/MAG001_genes.faa
```

### Outputs:

```text
results/MAG001_lowcov/
├── MAG001_GoM_1200m_gene_coverage.png
├── MAG001_GoM_1200m_sig_genes.tsv
└── MAG001_GoM_1200m.fasta
```

These are genes in the **lowest 2.5% of normalized coverage**.

---

## Step 5 — Identify differential gene content between populations

Now we compare **gene content differences** between the 1200 m sample and all other samples:

```bash
python scripts-populations/08d_Differential_Genes.py \
  -i coverage_magic/MAG001/ \
  -t GoM_1200m_2016-01-01 \
  -o results/MAG001_diffgenes/ \
  -s 0.025 \
  -a annotations/MAG001.annotations \
  -f genes/MAG001_genes.faa
```

### Outputs (per comparison sample):

```text
results/MAG001_diffgenes/
├── GoM_1200m-GoM_80m_diffcov_plot.png
├── GoM_1200m-GoM_10m_diffcov_plot.png
├── GoM_1200m-GoM_80m_diffcov_high.tsv
├── GoM_1200m-GoM_80m_diffcov_low.tsv
├── GoM_1200m-GoM_10m_diffcov_high.tsv
├── GoM_1200m-GoM_10m_diffcov_low.tsv
├── GoM_1200m-GoM_80m_diffcov_high.fa
├── GoM_1200m-GoM_80m_diffcov_low.fa
├── GoM_1200m-GoM_10m_diffcov_high.fa
└── GoM_1200m-GoM_10m_diffcov_low.fa
```

Each comparison yields:

* A **distribution plot** of per-gene coverage differences
* A list of **over-enriched genes** (top 2.5%)
* A list of **under-enriched genes** (bottom 2.5%)
* Optional **FASTA outputs** for downstream functional analysis

---

## Summary of the full pipeline for one rMAG

| Stage                 | Script                                 | Output                               |
| --------------------- | -------------------------------------- | ------------------------------------ |
| Population detection  | `06f_TabBlast_RecPlot_Mini_Auto_v4.py` | ANIr, MNIr, depth, breadth, RecPlots |
| ANIr stability        | `06h_PopX_Bootstrap_Test.py`           | Bootstrap distributions              |
| Population comparison | `06h_PopX_Permutation_Test.py`         | Same vs different population tests   |
| Low-coverage genes    | `08c_Low_Coverage_Genes.py`            | Underrepresented genes               |
| Differential genes    | `08d_Differential_Genes.py`            | Enriched & depleted genes            |

This workflow reproduces the population-tracking, statistical validation, and functional differentiation results described in the manuscript.

---

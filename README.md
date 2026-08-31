# ANUBIS

**Ancient DNA Damage Profiler for MetaPhlAn SGBs**

[![PyPI version](https://badge.fury.io/py/anubis-adna.svg)](https://pypi.org/project/anubis-adna/)
[![Python ≥3.10](https://img.shields.io/badge/python-%E2%89%A53.10-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENCE)

ANUBIS integrates MetaPhlAn shotgun metagenomics profiles with ancient DNA (aDNA) damage assessment. It takes a MetaPhlAn relative-abundance table and the accompanying bz2-compressed SAM alignment, groups reads per species-genome-bin (SGB), and fits the characteristic C→T / G→A deamination model using [PyDamage](https://github.com/maxibor/pydamage). Benjamini–Hochberg FDR correction is applied across all taxa and, optionally, base-quality scores are rescaled or terminal bases are masked for downstream variant calling with StrainPhlAn.

---

## Features

- Per-SGB damage assessment - C→T deamination at 5′ end, G→A at 3′ end
- Geometric decay damage model with confidence intervals (PyDamage)
- Benjamini–Hochberg FDR correction across all analysed SGBs
- Predicted accuracy score (logistic regression GLM)
- Classical aDNA damage profile plots (PNG)
- Model-based base-quality rescaling for StrainPhlAn genotyping (`--rescale`)
- Flexible terminal base masking (`--mask-5prime N`, `--mask-3prime N`)
- Clade-specific mode: target one or more SGBs directly (`-c`/`--clade`)
- SAM header index cache for fast repeated runs
- Multi-threaded BAM sorting and PyDamage analysis

---

## Requirements

### System dependencies

These must be on your `PATH` before running ANUBIS:

| Tool | Version | Install |
|------|---------|---------|
| `samtools` | ≥1.15 | `conda install -c bioconda samtools` |
| `bzip2` / `bzcat` | any | usually pre-installed; `conda install bzip2` |

### Python dependencies

Installed automatically by pip:
`numpy`, `pandas`, `pysam`, `tqdm`, `matplotlib`, `statsmodels`, `pydamage ≥1.0`

---

## Installation

### From PyPI (recommended)

```bash
pip install anubis-adna
```

### Conda + pip (recommended for bioinformatics environments)

```bash
conda create -n anubis -c bioconda -c conda-forge \
    python=3.11 samtools bzip2 pysam numpy pandas tqdm matplotlib statsmodels
conda activate anubis
pip install anubis-adna
```

### From GitHub (latest development version)

```bash
pip install git+https://github.com/Mummy-Studies/ANUBIS.git
```

After installation the `anubis` command is available on your `PATH`.

---

## Quick start

```bash
anubis \
    -m sample.metaphlan.tsv \
    -s sample_metaphlan.sam.bz2 \
    -o results/
```

---

## Usage

### Basic — analyse the most abundant SGBs

```bash
anubis \
    -m sample.metaphlan.tsv \
    -s sample_metaphlan.sam.bz2 \
    -o results/ \
    --top-n 20 \
    --plot
```

### Clade-specific analysis

Target specific SGBs regardless of their abundance ranking. All of the following
formats are accepted:

```bash
# bare SGB ID
anubis -m ... -s ... -o results/ -c SGB13165

# MetaPhlAn clade suffix
anubis -m ... -s ... -o results/ -c t__SGB13165

# multiple SGBs (space- or comma-separated)
anubis -m ... -s ... -o results/ -c SGB13165 SGB6653
anubis -m ... -s ... -o results/ -c SGB13165,SGB6653
```

### Rescale base qualities for StrainPhlAn genotyping

ANUBIS can modify the BAM files produced per SGB in two complementary ways.
Both options write sorted, indexed BAMs to `<outdir>/rescaled/` which can be
passed directly to StrainPhlAn's `sample2markers.py`.

```bash
# Model-based quality rescaling (only applied to SGBs with significant damage)
anubis -m ... -s ... -o results/ --rescale

# Terminal masking only (independent of damage significance)
anubis -m ... -s ... -o results/ --mask-5prime 5 --mask-3prime 5

# Both combined (rescale first, then mask terminals)
anubis -m ... -s ... -o results/ \
    --rescale \
    --mask-5prime 5 --mask-3prime 5
```

---

## Options

### Main arguments

| Flag | Description |
|------|-------------|
| `-m`, `--metaphlan` | MetaPhlAn output TSV (required) |
| `-s`, `--sam` | MetaPhlAn SAM file, bz2-compressed (required) |
| `-o`, `--outdir` | Output directory (default: `metaphlan2damage_out`) |

### Filters

| Flag | Default | Description |
|------|---------|-------------|
| `-c`, `--clade` | — | Restrict to one or more specific SGBs (bypasses `--min-abundance` and `--top-n`) |
| `--min-reads N` | — | Discard SGBs with fewer than N aligned reads (applied after damage analysis) |
| `--min-abundance F` | — | Minimum relative abundance % (e.g. `0.01`) |
| `--top-n N` | — | Analyse only the N most abundant SGBs |

### Damage analysis

| Flag | Default | Description |
|------|---------|-------------|
| `--wlen N` | 30 | Damage modelling window length (bp) |
| `--threads N` | 4 | CPU threads for BAM sorting and PyDamage |
| `--min_mapq_val N` | 30 | Minimum mapping quality (MAPQ) for a read to be kept |

### Output

| Flag | Default | Description |
|------|---------|-------------|
| `--plot` | off | Save per-SGB damage profile plots (PNG) to `<outdir>/plots/` |
| `--rescale` | off | Rescale base qualities using the fitted damage model |
| `--rescale-threshold F` | 0.5 | Minimum `predicted_accuracy` required for rescaling |
| `--rescale-alpha F` | 0.05 | Maximum q-value required for rescaling |
| `--mask-5prime N` | 0 | Zero base quality for first N bases of every mapped read |
| `--mask-3prime N` | 0 | Zero base quality for last N bases of every mapped read |
| `--keep-bams` | off | Keep intermediate per-SGB BAMs in the output directory |

---

## Output files

| File | Description |
|------|-------------|
| `anubis_results.tsv` | Per-SGB damage metrics table (see column guide below) |
| `plots/<SGB>.png` | C→T / G→A damage profile with model fits (`--plot`) |
| `rescaled/<SGB>.bam` + `.bai` | Processed BAMs for StrainPhlAn (`--rescale` or `--mask-*`) |
| `<sam>.anubis_idx.json.gz` | SAM header index cache (auto-created beside the SAM file) |

### `anubis_results.tsv` column guide

| Column | Description |
|--------|-------------|
| `sgb_id` | SGB identifier (e.g. `SGB13165`) |
| `species` | Species name from MetaPhlAn lineage |
| `relative_abundance` | MetaPhlAn relative abundance (%) |
| `nb_reads_aligned` | Reads mapped to this SGB's marker genes |
| `coverage` | Mean read coverage across marker genes |
| `reflen` | Total marker gene reference length (bp) |
| `damage_model_pmax` | Maximum damage rate at position 1 (C→T at 5′ end) |
| `damage_model_pmax_stdev` | Standard deviation of the maximum damage rate |
| `damage_model_p` | Geometric decay rate of the damage signal |
| `damage_model_pmin` | Background (asymptotic) deamination rate |
| `damage_model_pmin_stdev` | Standard deviation of the background deamination rate |
| `null_model_p0` | Null model substitution rate |
| `pvalue` | Likelihood-ratio test p-value (damage vs. null model) |
| `qvalue` | BH-adjusted p-value (FDR across all SGBs in this run) |
| `predicted_accuracy` | GLM-predicted reliability of the damage call (0–1) |
| `RMSE` | Root mean square error of the damage model fit |
| `CtoT-0` … `CtoT-29` | Observed C→T frequency at each position from the 5′ end |
| `GtoA-0` … `GtoA-29` | Observed G→A frequency at each position from the 3′ end |

---

## Examples

The `examples/` directory contains:

- **`merge_anubis_tables.py`** — concatenate `anubis_results.tsv` files from multiple
  samples into a single long-format table for cross-sample analyses
- **`plot_multiple_damage.rmd`** — R Markdown workflow for visualising damage profiles
  across many samples

---

## Citation

If you use ANUBIS in your research, please cite:

> Bello K, Segata N, Maixner F, Sarhan M.
> **ANUBIS: Ancient DNA damage profiler for MetaPhlAn SGBs.**
> GitHub: https://github.com/Mummy-Studies/ANUBIS (2025)

PyDamage should also be cited:

> Borry M, Hübner A, Rohrlach AB, Warinner C.
> **PyDamage: automated ancient damage identification and estimation for contigs in ancient DNA de novo assembly.**
> *PeerJ* 2021, 9:e11845. https://doi.org/10.7717/peerj.11845

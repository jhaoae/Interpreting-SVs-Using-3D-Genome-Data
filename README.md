# Interpreting Structural Variations Using 3D Genome Data

This repository contains three post-processing scripts used to associate structural
variants (SVs) with promoter contacts and local insulation-score changes. It does not
include raw sequencing/contact data or the upstream SV-calling and contact-matrix
generation pipelines.

## Installation

Python 3.10 or newer is recommended.

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
mkdir -p results
```

## 1. High-contact promoter bins

`high_contact_promoter_bin.py` finds promoter-overlapping bins among the strongest
balanced contacts of each non-translocation SV bin.

```bash
python3 high_contact_promoter_bin.py \
  --vcf sample.vcf \
  --cool sample.cool \
  --promoter promoter.bed \
  --out results/high_contact \
  --percentile 99.99 \
  --dist 100000
```

Inputs:

- VCF: standard tab-delimited VCF; `SVTYPE` and `END` are read from INFO.
- COOL: a fixed-bin cooler matrix with balancing weights.
- promoter BED: four columns (`chrom`, `start`, `end`, `gene`), without a header.

Outputs are `<out>_detail.tsv` and `<out>_summary.tsv`.

## 2. Insulation-score changes

`insulation_change.py` compares tumor and baseline insulation scores in windows around
SVs and performs a paired Wilcoxon test.

```bash
python3 insulation_change.py \
  tumor.bedgraph normal.bedgraph sv_windows.bed \
  --out-prefix results/SV_insulation \
  --quantile 0.05
```

The two bedGraph files contain `chrom`, `start`, `end`, and `score`. The SV BED file
contains `chrom`, `window_start`, `window_end`, `sv_type`, `sv_start`, and `sv_end`,
without a header.

Outputs are `<prefix>_change_detail_score.tsv`, `<prefix>_change.tsv`, and
`<prefix>_significant_change.tsv`.

## 3. Background permutation test

`permutation.py` compares each selected SV window with same-chromosome windows sampled
from a BED file of regions without SVs.

```bash
python3 permutation.py \
  results/SV_insulation_significant_change.tsv \
  no_sv_regions.bed tumor.bedgraph normal.bedgraph \
  --n-bg 500 \
  --seed 1 \
  --out results/SV_insulation_permutation_test.tsv
```

The first input is the significant-change table from the previous step. The no-SV BED
contains `chrom`, `start`, and `end`, without a header. Setting `--seed` makes the
background sampling reproducible.

## Scope and limitations

- The high-contact script reports contacts above a user-selected percentile; it does not
  implement an unaffected-locus null model.
- Background windows in `permutation.py` are matched by chromosome and window length,
  but not by GC content, mappability, or other genomic covariates.
- The scripts do not apply multiple-testing correction. Threshold sensitivity and
  additional context-matched null models should be evaluated before making inferential
  claims.

# Step C: Structural Variant Annotation

This step takes clustered structural variant (SV) candidates and annotates them using the GENCODE v43 reference genome.

## Overview

The goal is to classify each SV based on its relationship to known genes and exons. The pipeline identifies events such as gene fusions, exonic deletions, and distinguishes them from normal splicing. It uses a pre-computed IntervalTree for high-performance annotation lookups.

## Requirements

pandas v2.3.2+
intervaltree v3.1.0+
python3

## Workflow

The procecess is divided into two main stages: a one-time pre-computation and the main annotation script.

**1. Pre-computation (Run this once)**

First, prepare the annotation database.

```bash
# Define gene regions for candidate search and initialize transcript tree for later queries
python init_bed_and_tree.py --gtf GTF --outdir OUTDIR_FOR_BED_AND_TREE
```

This performs two key actions:

1.  **Builds a Gene BED File**: It creates `gencode.v43.annotation.genes.bed`, a file containing the genomic coordinates for every gene, which can be used for broad, gene-level SV searches.
2.  **Generates the Interval Tree Cache**: It produces `tx_tree_cache.pkl`. This file is a pickled object of interval trees for interval-based queries (e.g. chr1:130000-150000 should return a list of all overlapping exons across verified GENCODE transcripts.)

N.B.: These files are now placed under `resources/` for easy reference.

**2. Run Annotation**

Once the cache is built, run the main annotation script.

```bash
usage: annotate.py [-h] --candidates CANDIDATES --cache CACHE -o OUTDIR

Annotate SV candidates using prebuilt transcript trees.

options:
  -h, --help                    show this help message and exit
  --candidates CANDIDATES       Path to SV candidates
  --cache CACHE                 Transcript tree file path
  -o OUTDIR, --outdir OUTDIR    Output directory for annotated VCF
```

Output VCF is named `sv_candidates.annotated.vcf` under the output directory of choice.

## File Descriptions

- `scripts/`: Contains all executable scripts for this step.
- `resources/`: Contains input data and the pre-computed annotation tree (`tx_tree_cache.pkl`).
- `results/`: Contains the final annotated output files.

## Annotation Output

The main output of this step is `results/annotated_regions.csv`. This file contains the original SV candidate information from Step B, plus two additional columns: `Annotation` and `Gene(s)`.

The `Annotation` column can have one of the following values:

- **`Gene_Fusion(Fusion)`**: The SV's coordinates overlap with two or more distinct genes, suggesting a large-scale rearrangement.
  - 1. Large deletions spanning >1 genes
  - 2. INS of known gene regions (might need de novo assembly from RNAseq data? SR should be pretty noisy for this ...)
- **`Exonic_Deletion`**: The SV occurs within a single gene and completely contains at least one of its exons.
- **`Canonical_Splicing`**: The SV's coordinates precisely match the boundaries of a known intron (the region between two consecutive exons). This represents a standard splicing event.
- **`Intronic/Novel_Splicing`**: An SV that occurs within a single gene but does not fit the criteria for `Exonic_Deletion` or `Canonical_Splicing`. This indicates a potential genomic deletion within an intron or a non-standard splicing event.
- **`Intergenic`**: The SV does not overlap with any known gene.
  - Probably not something we want to prioritize?

# Step C: Structural Variant Annotation

This step takes clustered structural variant (SV) candidates and annotates them using the GENCODE v43 reference genome.

## Overview

The goal is to classify each SV based on its relationship to known genes and exons. The pipeline identifies events such as gene fusions, exonic deletions, and distinguishes them from normal splicing. It uses a pre-computed IntervalTree for high-performance annotation lookups.

## Workflow

The process is divided into two main stages: a one-time pre-computation and the main annotation script.

**1. Pre-computation (Run this once)**

First, prepare the annotation database.

```bash
# Filter the raw GTF file, generate resources/cleaned_annotations.gtf
bash scripts/01_preprocess.sh

# Build the fast lookup tree, using resources/cleaned_annotations.gtf
python scripts/02_build_gene_table_and_trees.py
```
This will generate the `tx_tree_cache.pkl` file in the `resources/` directory (fast, reusable lookup database - the IntervalTree).

**2. Run Annotation**

Once the cache is built, run the main annotation script.

```bash
# Usage:
python scripts/annotate_svs.py \
    --input <path_to_sv_candidates.csv> \
    --tree resources/tx_tree_cache.pkl \
    --output results/annotated_regions.csv
```

## File Descriptions

-   `scripts/`: Contains all executable scripts for this step.
-   `resources/`: Contains input data and the pre-computed annotation tree (`tx_tree_cache.pkl`).
-   `results/`: Contains the final annotated output files.
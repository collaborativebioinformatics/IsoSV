# Step C: Structural Variant Annotation

This step takes clustered structural variant (SV) candidates and annotates them using the GENCODE v43 reference genome.

## Overview

The goal is to classify each SV based on its relationship to known genes and exons. The pipeline identifies events such as gene fusions, exonic deletions, and distinguishes them from normal splicing. It uses a pre-computed IntervalTree for high-performance annotation lookups.

## Workflow

The process is divided into two main stages: a one-time pre-computation and the main annotation script.

**1. Pre-computation (Run this once)**

First, prepare the annotation database.

```bash
# 1. Filter the raw GTF to create a clean, high-confidence annotation set.
# This generates resources/cleaned_annotations.gtf
bash scripts/01_preprocess.sh

# 2. Process the cleaned GTF to build fast lookup databases.
python scripts/02_build_gene_table_and_trees.py
```
This performs two key actions:
1.  **Builds a Gene BED File**: It creates `resources/gencode.v43.annotation.genes.bed`, a file containing the genomic coordinates for every gene, which can be used for broad, gene-level SV searches.
2.  **Generates the Interval Tree Cache**: It produces `resources/tx_tree_cache.pkl`. This file is a pre-compiled lookup database (an IntervalTree) that maps every transcript and its individual exons to their precise genomic locations, enabling high-speed annotation queries.

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

## Annotation Output

The main output of this step is `results/annotated_regions.csv`. This file contains the original SV candidate information from Step B, plus two additional columns: `Annotation` and `Gene(s)`.

The `Annotation` column can have one of the following values:

-   **`Gene_Fusion`**: The SV's coordinates overlap with two or more distinct genes, suggesting a large-scale rearrangement.
-   **`Exonic_Deletion`**: The SV occurs within a single gene and completely contains at least one of its exons.
-   **`Canonical_Splicing`**: The SV's coordinates precisely match the boundaries of a known intron (the region between two consecutive exons). This represents a standard splicing event.
-   **`Intronic/Novel_Splicing`**: An SV that occurs within a single gene but does not fit the criteria for `Exonic_Deletion` or `Canonical_Splicing`. This indicates a potential genomic deletion within an intron or a non-standard splicing event.
-   **`Intergenic`**: The SV does not overlap with any known gene.
# Step C: Structural Variant Annotation

This step takes clustered structural variant (SV) candidates and annotates them using the GENCODE v43 reference genome.

## Overview

The goal is to classify each SV based on its relationship to known genes and exons. The pipeline identifies events such as gene fusions, exonic deletions, and distinguishes them from normal splicing. It uses a pre-computed IntervalTree for high-performance annotation lookups.

## Requirements

pandas v2.3.2+
intervaltree v3.1.0+
python3

## Workflow

The process is divided into two main stages: a one-time pre-computation and the main annotation script.

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
# Usage:
python scripts/annotate.py \
    --candidates <path_to_sv_candidates.tsv> \
    --cache resources/tx_tree_cache.pkl \
    -o <output_directory>

# Example with mock data from Step B:
python scripts/annotate.py \
    --candidates ../step_b_IsoClustering/test_data/results/GIAB002_chr22_region_LongReadSVs_merged_SVs.tsv \
    --cache resources/tx_tree_cache.pkl \
    -o results/
```
The script will generate an annotated VCF file named `GIAB002_chr22_region_LongReadSV.annotated.vcf` inside the specified output directory.

## File Descriptions

- `scripts/`: Contains all executable scripts for this step.
- `resources/`: Contains input data and the pre-computed annotation tree (`tx_tree_cache.pkl`).
- `results/`: Contains the final annotated output files.

## Annotation Output

The main output is an annotated VCF file. The `INFO` column of the VCF contains our functional annotations.

-   **`SVTYPE`**: This standard field retains the original structural variant type (`DEL` or `INS`).
-   **`ANN_TYPE`**: This is a custom field that provides the detailed functional classification.

The possible values for `ANN_TYPE` depend on the original `SVTYPE`:

#### For Deletions (`SVTYPE=DEL`):
-   **`Gene_Fusion`**: The deletion's coordinates overlap with two or more distinct genes.
-   **`Exonic_Deletion`**: The deletion occurs within a single gene and completely contains at least one of its exons.
-   **`Canonical_Splicing`**: The deletion's coordinates precisely match the boundaries of a known intron.
-   **`Intronic/Novel_Splicing`**: The deletion occurs within a single gene but does not fit the criteria above.
-   **`Intergenic`**: The deletion does not overlap with any known gene.

#### For Insertions (`SVTYPE=INS`):
-   **`Insertion_in_Fusion_Region`**: The insertion point falls in a complex region with overlapping annotations for two or more genes.
-   **`Intragenic_Insertion`**: The insertion point falls within the boundaries of a single gene.
-   **`Intergenic`**: The insertion point does not overlap with any known gene.

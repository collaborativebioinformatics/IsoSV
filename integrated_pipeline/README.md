# IsoSV Integrated Pipeline

This directory contains the integrated pipeline for IsoSV, combining parsing, clustering, and annotation into a single step.

## Requirements

- Python 3.9+
- `pysam`
- `intervaltree`
- `pandas`

Install the requirements using conda:
```bash
conda create -n isosv_pipeline python=3.9 pysam intervaltree pandas -c bioconda -c conda-forge
conda activate isosv_pipeline
```

## Usage

```bash
python integrated_pipeline/scripts/run_pipeline.py <input.bam> -o <output_prefix> [options]
```

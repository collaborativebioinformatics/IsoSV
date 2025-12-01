# IsoSV Integrated Pipeline

This directory contains the integrated pipeline for IsoSV, combining parsing, clustering, and annotation into a single step.

## Requirements

- Python 3.9+
- `pysam`
- `intervaltree`

Install the requirements using conda:
```bash
conda create -n isosv python=3.9 pysam intervaltree -c bioconda -c conda-forge
conda activate isosv
```

## Usage

```bash
chmod +x ./IsoSV
usage: IsoSV [-h] -o OUT_PREFIX --gtf GTF [--reference REFERENCE] [--min-mapq MIN_MAPQ]
             [--min-ins MIN_INS] [--min-del MIN_DEL] [--min-clip MIN_CLIP] [--bp-window BP_WINDOW]
             [--save-parser-output] [--max-readnames MAX_READNAMES]
             in_bam

Integrated pipeline for RNA-seq based SV detection, clustering, and annotation.

optional arguments:
  -h, --help                show this help message and exit

Positional arguments:
  in_bam                    Input BAM/CRAM

Required arguments:
  -o OUT_PREFIX, --out-prefix OUT_PREFIX
                            Output prefix for TSVs
  --gtf GTF                 Path to plain GTF file for building gene annotation BED and cache
  --reference REFERENCE     Reference FASTA (required for CRAM)

Options:
  --min-mapq MIN_MAPQ
  --min-ins MIN_INS         Minimum insertion length to report
  --min-del MIN_DEL         Minimum deletion length to report
  --min-clip MIN_CLIP       Minimum soft/hard clip length to report
  --bp-window BP_WINDOW     Breakpoint clustering window (bp)
```

# isoSV_BAMparser

- [Google Doc - Daily Schedule](https://docs.google.com/document/d/1EMqbb5DUvDwu5YHkBu7oTPW8S4y2RMiDstHkjeXXnfE/edit?usp=sharing)
- [Google Doc - Team 9](https://docs.google.com/document/d/1i5qklL01o8b1E8FYtd3IBWgXcedRDMfNrrp3aHU8AwE/edit?tab=t.0)
- [SV Hackathon 2025](https://fritzsedlazeck.github.io/blog/2025/hackathon-2025/)

**Long Read**
- [GIAB002_chr22_region.bam](https://drive.google.com/drive/folders/1y48dxKJYkRXxDcEt6kTNbaEs6qCvUD8I?usp=drive_link)

**Short Read**
- [chr22.bam](https://drive.google.com/drive/folders/1udSRBNAhaS4xwd-hQ0pE4PgLi8kSQbXB)

---

## Aim : Detect SV-sized indels and candidate breakpoints in short-read RNA/DNA BAMs using:

- CIGAR `I`/ `D` (ins/del)
- Soft/hard clips at read ends (`S` / `H`)
- Split reads via `SA` tag (supplementary alignment)

---

## Short Read:
* `sr_isoSV_parser.py` — short-read focused parser
## Long Read:
* `lr_isoSV_parser.py`  - 

### Output:
  - `<prefix>.indels.tsv`
  - `<prefix>.softclips.tsv`
  - `<prefix>.splits.tsv`

---

## Requirements

* Python 3.8+
* `pysam` (install via pip or conda)
* `samtools` (recommended for indexing / manual inspection — optional for running the parser)

Minimal `pip` install:

```bash
python -m pip install --upgrade pip
pip install pysam
```

Or with conda:

```bash
conda create -n isoparser python=3.9 pysam -c bioconda -c conda-forge
conda activate isoparser
```

---

## How to download example BAM files

Below are a few example commands to download common small example BAMs you can use for testing. Replace URLs if you have other data sources.

### (1) UCSC `bamExample.bam` (small Illumina example)

```bash
wget http://genome.ucsc.edu/goldenPath/help/examples/bamExample.bam
wget http://genome.ucsc.edu/goldenPath/help/examples/bamExample.bam.bai
```

Notes:

* Older example files may contain auxiliary-field oddities; the parser is defensive and will still operate even if `samtools view` sometimes prints errors.
* If `samtools` complains about corruption when streaming, you can optionally create a cleaned BAM using a helper (not required for `iso_parser.py`).

### (Haven't tested (2) and (3) below)
##### (2)  Gencode v43 annotation (if you want to annotate clusters)

```bash
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz
```

Store that file somewhere accessible and pass `--gtf /path/to/gencode.v43.annotation.gtf.gz` to the parser if you want gene/exon annotation for clusters.

##### (3) Example short/long-read BAMs from public archives

* NCBI SRA: use `prefetch` + `fasterq-dump` + alignment, or download pre-aligned BAMs if available from project pages.
* ENA/1000 Genomes: many projects host BAM/CRAM; use `wget`/`curl` with direct links or use the ENA browser.

---

## How to run (examples)

Basic:

```bash
python sr_isoSV_parser.py input.bam -o results
```

**With gene annotation (if gene BED uses column 7 for gene symbol, for example):**
```bash
python sr_isoSV_parser.py input.bam -o results \
  --gene-bed gencode.v43.annotation.genes.bed --gene-name-col 7
```

**Quick debug (recommended SR-friendly small thresholds):**
```bash
python sr_isoSV_parser.py input.bam -o results \
  --min-ins 1 --min-del 1 --min-clip 1 --bp-window 25
```

**Full run example:**
```bash
python sr_isoSV_parser.py chr22.bam -o chr22_ShortReadSVs \
  --min-ins 20 --min-del 20 --min-clip 20 --bp-window 50 \
  --gene-bed gencode.v43.annotation.genes.bed --gene-name-col 7
```

Short reads need different defaults than long reads — suggestions:

- Small indel discovery (sensitive): `--min-ins 1 --min-del 1` (great for tiny indels, expect many calls)
- Balanced / typical SR indel detection: `--min-ins 5 --min-del 5`
- Conservative (reduce false positives): `--min-ins 20 --min-del 20`
- Soft-clip: `--min-clip` often ~ `8–25` for SR depending on read length and goal
- Clustering window: `--bp-window` ~ `8–25` (tighter than LR script).
- Internal S (soft-clip): for short reads you usually care about end-clips; the script already reports left/right end clips and ignores internal S as breakpoint signal by default.

## Outputs and column definitions


When the parser finishes it writes a per-read candidates TSV and several helper files. The per-read TSV has the following header columns:




### `<prefix>.indels.tsv`
```
chrom  pos  end  svtype  length  support  median_mapq  example_reads  example_ins_seq  genes
```

Column meanings:
- `chrom` — reference (e.g. `chr17` or `17`)
- `pos` — 1-based start (INS: base before inserted seq; DEL: first deleted base)
- `end` — inclusive end (for DEL: pos + length - 1; for INS: = pos)
- `svtype` — `INS` or `DEL`
  - `INS` — insertion present in the read relative to the reference (CIGAR `I`)
  - `DEL` — deletion present in the read relative to the reference (CIGAR `D`)
- `length` — integer length
- `support` — cluster support (distinct reads)
- `median_mapq` — median MAPQ of cluster
- `example_reads` — comma-separated supporting read names (max limited via --max-readnames)
- `example_ins_seq` — an example inserted sequence if available
- `genes` — comma-separated overlapping gene names (empty if none)

### `<prefix>.shortclips.tsv`
```
chrom  pos  side  support  median_mapq  median_cliplen  example_reads  genes
```

Column meanings:
- `side` — `L` or `R` (left/right clip)
- `pos` — 1-based coordinate for the clip (left = reference_start+1, right = reference_end)
- `median_cliplen` — median clip length for the cluster

### `<prefix>.splits.tsv`
```
chrom1  pos1  chrom2  pos2  support  median_mapq  notes  example_reads  genes_left  genes_right
```

Column meanings:
- `chrom1,pos1` — primary side breakpoint (1-based)
- `chrom2,pos2` — supplementary alignment side (1-based)
- `genes_left` / `genes_right` — overlapping gene names on each side (empty string if none)
  - exact overlap only (no nearest-gene fallback by default). The gene loader tries to match `chr` vs non-`chr` names and will pull gene names from the BED column you specify.

### ORIENTATION
- strand orientation of read and mate in `READ/MATE` form using `F` (forward) and `R` (reverse)`. See the orientation mapping below.

Below is the shorthand → TSV mapping we use:
  - `+`    = forward strand → `F`
  - `-`    = reverse strand → `R`
  - `+/-`  = read on forward, mate on reverse → `F/R`
  - `-/+`  = read on reverse, mate on forward → `R/F`
  - `+/+`  = both forward → `F/F`
  - `-/-`  = both reverse → `R/R`
- **SEQ** — read sequence (query sequence) if available, or `.` when missing.

### Quick check: indel size distribution from CIGAR
Use this to see how many insertions (`I`) and deletions (`D`) of each size exist in a BAM (quick way to pick thresholds):

Perl one-liner (fast & compact):
```bash
samtools view bamExample.bam \
  | perl -ne 'chomp; @f=split(/\t/); $c=$f[5]; while($c=~/(\d+)([ID])/g){ print "$2\t$1\n"; }' \
  | sort | uniq -c | sort -nr | head
```

---

### > To test: `python sr_isoSV_parser.py chr22.bam -o chr22_ShortReadSVs \--min-ins 20 --min-del 20 --min-clip 20 --bp-window 50 --gene-bed gencode.v43.annotation.genes.bed --gene-name-col 7`

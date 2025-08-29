# isoSV_BAMparser

- [Google Doc - Daily Schedule](https://docs.google.com/document/d/1EMqbb5DUvDwu5YHkBu7oTPW8S4y2RMiDstHkjeXXnfE/edit?usp=sharing)
- [Google Doc - Team 9](https://docs.google.com/document/d/1i5qklL01o8b1E8FYtd3IBWgXcedRDMfNrrp3aHU8AwE/edit?tab=t.0)
- [SV Hackathon 2025](https://fritzsedlazeck.github.io/blog/2025/hackathon-2025/)

- [GIAB002_chr22_region.bam](https://drive.google.com/drive/folders/1y48dxKJYkRXxDcEt6kTNbaEs6qCvUD8I?usp=drive_link)

---

* `iso_parser.py` — the single main script. It performs robust parsing, per-read candidate extraction, optional clustering, and optional GTF annotation (if you supply a GTF).

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


**Minimal (write outputs next to `results.tsv`):**


```bash
python iso_parser.py bamExample.bam results.tsv
```


**Put all outputs into `results/` directory:**


```bash
python iso_parser.py bamExample.bam results.tsv --out-dir results/
```


**Use low thresholds for quick debugging + cluster with support >=2:**


```bash
python iso_parser.py bamExample.bam demo_candidates.tsv --out-dir results/ \
--insert 1 --delete 1 --softclip 1 --sample-parse 2000 \
--cluster-window 50 --min-support 2 --salvage
```

**Report discordant paired-end events (mate-unmapped, not-proper, large insert):**
```bash
python iso_parser.py input.bam results.tsv --report-discordant --out-dir results/
```

## Outputs and column definitions


When the parser finishes it writes a per-read candidates TSV and several helper files. The per-read TSV has the following header columns:


```
chrom	pos_start	pos_end	read_name	sv_type	sv_len	note	read_len	cigar	MAPQ	TLEN	MATE_UNMAPPED	IS_PROPER_PAIR	ORIENTATION	SEQ
```


Column meanings:
- **chrom** — reference sequence name (e.g. `1`, `chr1`, `21`, `X`) where the read aligned.
- **pos_start** — 0-based reference coordinate where the event begins (for insertions this is the base *before* the inserted sequence; for deletions it is the first deleted base; for soft-clips it's the position on the reference adjacent to the clipped sequence).
  - For `INS`: the base immediately before the inserted sequence (insertion occurs between bases), so pos_end == pos_start.
  - For `DEL`: first deleted base (0-based).
  - For `SOFTCLIP`: reference position adjacent to clipped bases (left/right/internal note further clarifies).
  - For `LARGE_N`: start of the skipped region (N) on the reference.

- **pos_end** — 0-based inclusive end coordinate of the event on the reference. For insertions this is equal to `pos_start` (the insertion occurs between bases). For events that span bases (deletions, LARGE_N), this is `pos_start + sv_len - 1`.
- **read_name** — query/read identifier (QNAME) from the BAM (useful to inspect supporting reads manually in IGV or samtools).
- **sv_type** — event type detected in the read. Typical values:
  - `INS` — insertion present in the read relative to the reference (CIGAR `I`)
  - `DEL` — deletion present in the read relative to the reference (CIGAR `D`)
  - `SOFTCLIP` — soft-clipped bases at the read ends (CIGAR `S`) — possible breakpoint
  - `LARGE_N` — very large `N` (skipped region) in the CIGAR (often large introns, reported only if `--large-N` threshold is set)
  - `SPLIT` — split / supplementary alignment indicated via `SA` tag
  - `MATE_UNMAPPED` — paired read mate is unmapped
  - `DISC_PAIR` — discordant pair (not proper)
  - `LARGE_INSERT` — abnormally large template length (TLEN)

- **sv_len** — integer length of the event reported (number of inserted or deleted bases, soft-clip length, or skipped `N` length). For `SPLIT` this will be `0` and the `note` explains the cause.
- **note** — short human-readable note describing why the read was reported (e.g. `CIGAR_I`, `CIGAR_D`, `CIGAR_S_left`, `SA_tag`, or `SALVAGED_CIGAR_I`). Salvaged detections (from reads that caused parse errors) are prefixed with `SALVAGED_`.
- **read_len** — length of the read sequence (number of bases in the read). This helps to disambiguate short vs long reads and to assess reliability of soft-clip signals.
- **cigar** — the exact CIGAR string for the alignment as present in the BAM (e.g. `12M1I23M`). This is written verbatim so you can reproduce the detection logic later.
- **MAPQ** — mapping quality of the read (integer). A missing/unknown value is given as `.`.
- **TLEN** — template length (TLEN/TLEN field from BAM). For paired-end reads this is the signed insert length; stored as absolute for reporting where relevant; missing/zero often shown as `.`.
- **MATE_UNMAPPED** — `1` if mate is unmapped, `0` otherwise.
- **IS_PROPER_PAIR** — `1` if `is_proper_pair` is true (flag 0x2), `0` otherwise.
- **ORIENTATION** — strand orientation of read and mate in `READ/MATE` form using `F` (forward) and `R` (reverse)`. See the orientation mapping below.
  Below is the shorthand → TSV mapping we use:
  - +    = forward strand → F
  - -    = reverse strand → R
  - +/-  = read on forward, mate on reverse → F/R
  - -/+  = read on reverse, mate on forward → R/F
  - +/+  = both forward → F/F
  - -/-  = both reverse → R/R
- **SEQ** — read sequence (query sequence) if available, or `.` when missing.








## Other files


- `candidates.tsv.errors.log` — logged parse and salvage errors (if any). Useful for debugging malformed BAMs.
- `candidates.tsv.clusters.tsv` and `candidates.tsv.clusters.txt` — output of per-read clustering (identical contents, TSV with cluster coordinates, type, support and list of supporting read names). Generated only when `--cluster-window` is used.
- `candidates.tsv.summary.txt` — short JSON + human-readable summary of counts (total reads examined, total candidates, per-type counts, skipped reads).

<!-- ### Cluster and annotate with Gencode v43 (example):
``` bash 
python iso_parser.py input.bam candidates.tsv --cluster-window 50 --min-support 2 --gtf gencode.v43.annotation.gtf.gz
``` -->

### Quick check: indel size distribution from CIGAR

Use this to see how many insertions (`I`) and deletions (`D`) of each size exist in a BAM (quick way to pick thresholds):

Perl one-liner (fast & compact):
```bash
samtools view bamExample.bam \
  | perl -ne 'chomp; @f=split(/\t/); $c=$f[5]; while($c=~/(\d+)([ID])/g){ print "$2\t$1\n"; }' \
  | sort | uniq -c | sort -nr | head
```

---

> To test: `python iso_parser.py bamExample.bam demo_candidates.tsv --insert 1 --delete 1 --softclip 1 --sample-parse 2000`

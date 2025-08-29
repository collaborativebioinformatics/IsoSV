# IsoSV
A pipeline for detecting structural variants from RNA-Seq Data

<img width="250" height="250" alt="image" src="https://github.com/user-attachments/assets/0a7a755e-688c-418d-8018-4077c9115364" />

## Team
## Team
- [Fritz Sedlazeck](https://github.com/fritzsedlazeck)
- [Rupesh Kesharwani](https://github.com/unique379r)
- [Van Truong](https://github.com/van-truong)
- [Minhang Xu](https://github.com/MinhangXu)
- [Memoona Rasheed](https://github.com/MemoonaRasheed)
- [Aisha Yousaf](https://github.com/AishaYousaf)
- [Bigy Ambat](https://github.com/bigyambat)
- [Louis She](https://github.com/snakesch)
- [Pu Kao (Paul)](https://github.com/isthatgopro)
- [Bharati Jadhav](https://github.com/bharatij)
- [Saolendra](https://github.com/sailepradh)
- [Siyu Wang]
- [Farhang Jaryani]
- [Kirtan Dave]
- [Christopher M. Grochowski](https://github.com/cgrochowski)
- [Yousuf Bahit]


## Overview

Detecting structural variants (SVs) from RNA-seq data presents unique challenges. Unlike DNA sequencing, RNA reads span spliced transcripts, resulting in complex CIGAR patterns that include skipped regions (N), soft-clips, insertions, deletions, and split alignments. Standard DNA-based SV callers often misinterpret these signals, leading to missed or misclassified events. To overcome this, we developed a pipeline i.e. IsoSV that scans RNA-seq BAM files to identify candidate SVs by parsing CIGAR operations, split-read (SA) tags, and, optionally, exon annotations. The method distinguishes expected introns from potential novel splice junctions or structural rearrangements, reporting each event in both TSV and VCF formats to enable comprehensive detection and downstream analysis of SVs from RNA-seq data. At the end we are validating these SV finding with known DNA variants. It involves three steps: [IsoParser](https://github.com/collaborativebioinformatics/IsoSV/tree/main/step_a_IsoParser), [IsoClustering](https://github.com/collaborativebioinformatics/IsoSV/tree/main/step_b_IsoClustering) and [IsoAnnotator](https://github.com/collaborativebioinformatics/IsoSV/tree/main/step_c_IsoAnnotator).

## Gene/Transcript Fusion vs RNA SV

RNA structural variants are any transcript-level rearrangements observed in RNA-seq reads, whereas transcript fusions are specific chimeric transcripts joining exons from two separate genes, often reflecting underlying DNA rearrangements.

## Installation

### 🚀 **Getting Started**

IsoSV is implemented as a three-stage workflow:

Step 1: IsoParser - Identifies candidate events by scanning BAM alignments for signatures of SVs encoded in CIGAR strings and supplementary alignment tags. It flags long insertions, deletions, skipped regions, and soft clips while recording read-level support.

[Running IsoParser](https://github.com/collaborativebioinformatics/IsoSV/blob/main/step_a_IsoParser/README.md)

Step 2: IsoClustering - Processes candidate events to merge overlapping signals across reads. It uses an interval-tree data structure, querying and consolidation of nearby events. This clustering step is used to distinguish true biological SVs from noise in alignment.

[Running IsoClustering](https://github.com/collaborativebioinformatics/IsoSV/blob/main/step_b_IsoClustering/README.md)

Step 3: IsoAnnotator - overlays SVs with transcript annotations, classifying them as exon deletions, gene fusions, canonical splice events, or intronic rearrangements. This  annotation ensures that SV calls are interpretable and aligned with known genes.

[Running IsoAnnotator](https://github.com/collaborativebioinformatics/IsoSV/blob/main/step_c_IsoAnnotator/readme.md)




    ```

## Workflow

The IsoSV workflow for structural variant analysis was designed to identify and evaluate candidate variants from RNA sequencing data. The process begins with parsing BAM files from RNASeq data (Long Read and Short Read), where input BAM alignments are filtered based on mapping quality (XXXX), and candidate variants larger than 30-50 bp are extracted. These candidates are exported into text or VCF files. Next, a data structure for genomic intervals is used to cluster similar candidate entries, generating consolidated variant calls in text or VCF format. To enrich for biological relevance, candidate regions are prioritized using external known annotation resources such as BED or GFF files are incorporated to annotate candidates, producing an updated VCF with annotated structural variants. For visualization, the resulting VCF files are inspected in IGV. Finally, benchmark datasets (HG002) are used for evaluation as a truthset using bedtools and RNASeq coverage profiles to validate expression of the candidate events. To benchmark performance, we constructed a GIAB truth set of large indels by combining two sources: (i) structural variant calls annotated with SVTYPE and SVLEN, and (ii) indel calls derived from reference/alternative allele length differences (≥30 bp). Each variant was represented as a ±10 bp window in BED format. Overlap between Clair3-RNA candidate large indels and GIAB truth sets was assessed using bedtools intersect. Validation was defined as any overlap between candidate and GIAB indels within this positional tolerance.

<img width="1026" height="799" alt="IsoSV" src="https://github.com/user-attachments/assets/38bc9eb9-8306-4a79-b0bf-d975375897a9" />

## Presentation

https://docs.google.com/presentation/d/1-pTwId0y6V8OCrv-FYEJuwpxVq7wN8G9hY5ynpy-XS0/edit?usp=sharing

## Examples

<img width="1466" height="870" alt="Screenshot 2025-08-28 at 3 17 39 PM" src="https://github.com/user-attachments/assets/49a444be-b40e-4d44-9d59-4b623191f727" />







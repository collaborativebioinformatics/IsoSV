# IsoSV

<img width="250" height="250" alt="image" src="https://github.com/user-attachments/assets/0a7a755e-688c-418d-8018-4077c9115364" />

## Team
Fritz Sedlazeck
Christopher Grochowski 
Rupesh Kesharwani
Siyu Wang
Farhang Jaryani
Van Truong
Minhang Xu
Memoona Rasheed
Kirtan Dave
Aisha Yousaf
Bigy Ambat
Louis SHE
Pu Kao (Paul)
Yousuf Bahit
Bharati Jadhav 

## Overview

Detecting structural variants (SVs) from RNA-seq data presents unique challenges. Unlike DNA sequencing, RNA reads span spliced transcripts, resulting in complex CIGAR patterns that include skipped regions (N), soft-clips, insertions, deletions, and split alignments. Standard DNA-based SV callers often misinterpret these signals, leading to missed or misclassified events. To overcome this, we developed a pipeline i.e. IsoSV that scans RNA-seq BAM files to identify candidate SVs by parsing CIGAR operations, split-read (SA) tags, and, optionally, exon annotations. The method distinguishes expected introns from potential novel splice junctions or structural rearrangements, reporting each event in both TSV and VCF formats to enable comprehensive detection and downstream analysis of SVs from RNA-seq data. At the end we are validating these SV finding with known DNA variants. 

## Gene/Transcript Fusion vs RNA SV
RNA structural variants are any transcript-level rearrangements observed in RNA-seq reads, whereas transcript fusions are specific chimeric transcripts joining exons from two separate genes, often reflecting underlying DNA rearrangements.

## Presentation

https://docs.google.com/presentation/d/1-pTwId0y6V8OCrv-FYEJuwpxVq7wN8G9hY5ynpy-XS0/edit?usp=sharing


## Installation
### 🚀 **Getting Started**

*(This section will be updated soon)*

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/collaborativebioinformatics/IsoSV.git
    ```
2.  **Set up the environment:**
    ```bash
    # Command to be added
    ```
3.  **Run the pipeline:**
    ```bash
    # Command to be added
    ```




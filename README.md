# IsoSV

<img width="250" height="250" alt="image" src="https://github.com/user-attachments/assets/0a7a755e-688c-418d-8018-4077c9115364" />

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
- [Christopher Grochowski]
- [Yousuf Bahit]

## Overview

Detecting structural variants (SVs) from RNA-seq data presents unique challenges. Unlike DNA sequencing, RNA reads span spliced transcripts, resulting in complex CIGAR patterns that include skipped regions (N), soft-clips, insertions, deletions, and split alignments. Standard DNA-based SV callers often misinterpret these signals, leading to missed or misclassified events. To overcome this, we developed a pipeline i.e. IsoSV that scans RNA-seq BAM files to identify candidate SVs by parsing CIGAR operations, split-read (SA) tags, and, optionally, exon annotations. The method distinguishes expected introns from potential novel splice junctions or structural rearrangements, reporting each event in both TSV and VCF formats to enable comprehensive detection and downstream analysis of SVs from RNA-seq data. At the end we are validating these SV finding with known DNA variants. 

## Gene/Transcript Fusion vs RNA SV

RNA structural variants are any transcript-level rearrangements observed in RNA-seq reads, whereas transcript fusions are specific chimeric transcripts joining exons from two separate genes, often reflecting underlying DNA rearrangements.

## Workflow

<img width="1026" height="799" alt="IsoSV" src="https://github.com/user-attachments/assets/38bc9eb9-8306-4a79-b0bf-d975375897a9" />

## Presentation

https://docs.google.com/presentation/d/1-pTwId0y6V8OCrv-FYEJuwpxVq7wN8G9hY5ynpy-XS0/edit?usp=sharing


## Installation

### 🚀 **Getting Started**

*(This section will be updated soon)*

1.  **Clone the repository:**
   
    ```bash
    git clone https://github.com/collaborativebioinformatics/IsoSV.git
    
    ```
3.  **Set up the environment:**
    ```bash
    # Command to be added
    ```
4.  **Run the pipeline:**
    ```bash
    # Command to be added
    ```




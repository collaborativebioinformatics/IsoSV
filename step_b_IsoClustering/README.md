# IsoSV: RNA-seq Structural Variant Discovery

## Overview

IsoSV is a Python-based pipeline for discovering structural variants (SVs) from RNA-seq alignments. This implementation focuses on **Sections 3 & 4** of the project: **Candidate Structuring (Intervals)** and **Clustering & Deduplication**.

## What You've Built

### 1. Interval Tree Infrastructure (`code/isosv/struct/intervals.py`)

- **`Candidate`** dataclass: Represents per-read SV candidates with genomic coordinates, SV type, and metadata
- **`ClusterSV`** dataclass: Represents clustered/consensus SVs with support counts and confidence intervals
- **`CandidateStructurer`** class: Manages candidates using interval trees for efficient spatial queries
- **`intervaltree` package integration**: Robust interval tree implementation for fast overlap queries

### 2. Testing & Validation (`code/main.py`)

- **Multi-format support**: Handles per-read candidates and clustered summary formats
- **Auto-detection**: Automatically detects input TSV format based on column headers
- **Interval tree testing**: Loads real SV data and tests tree operations
- **Advanced clustering**: Uses highest support read as baseline with exact coordinate positioning

## Installation

### Prerequisites
- Python 3.7+
- pip package manager

### Install Dependencies
```bash
pip install intervaltree
```

**Note**: The `pandas` dependency has been removed in favor of robust manual TSV parsing.

### Verify Installation
```bash
python -c "import intervaltree; print('intervaltree installed successfully')"
```

## Key Features

### Interval Tree Operations
- **Fast spatial queries** by chromosome and SV type
- **Efficient overlap detection** for clustering candidates
- **Point queries** for specific genomic positions
- **Region queries** for genomic intervals

### Data Handling
- **Flexible input formats**: 
  - Per-read candidate format (detailed read-level data)
  - Clustered summary format (consensus calls)
  - Auto-detection of input format
- **Automatic coordinate normalization**: Ensures consistent chromosome naming
- **INS expansion**: Converts point insertions to intervals using length information
- **Robust parsing**: Handles comments, missing values, and inconsistent column counts

### Clustering Capabilities
- **Configurable window size** for overlap detection and read grouping
- **Minimum support threshold** for filtering low-confidence clusters
- **Baseline selection**: Automatically selects highest support read as cluster reference
- **Exact positioning**: Cluster boundaries match baseline read coordinates precisely
- **Automatic directory creation** for output files
- **TSV output** in standardized format with clustering metadata

## File Structure

```
step_b_IsoClustering/
├── code/
│   ├── isosv/
│   │   └── struct/
│   │       └── intervals.py          # Core interval tree implementation
│   └── main.py                       # Main testing and clustering tool
├── test_data/
│   └── chr21_SVs.txt                 # Example SV data
└── README.md
```

## Data Formats

### Input: Per-Read Candidate Format
```
chrom	pos_start	pos_end	read_name	sv_type	sv_len	note	read_len	cigar	MAPQ	TLEN	MATE_UNMAPPED	IS_PROPER_PAIR	ORIENTATION	SEQ
chr22	17062902	17062902	read1	INS	322	note1	150	50M	60	300	0	1	F/F	ATCG...
```

### Input: Clustered Summary Format
```
chrom	pos	end	svtype	length	support	median_mapq	example_reads	example_ins_seq	genes
chr22	17062902	17062902	INS	322	3	60	read1,read2,read3	ATCG...	gene1
```

### Output: Merged Clusters Format
```
clusterid	chrom	start	end	type	support	median_sv_len	read_positions	baseline_read	baseline_support	window_expansion
cluster_001	chr22	17062902	17062902	INS	6	322	chr22:17062900:17062902,chr22:17062901:17062903,chr22:17062902:17062904	read1	3	0
```

**Note**: The `support` column shows the **sum of support values** from all reads in the cluster, not just the count of reads. For example, if 2 reads with support values 3 and 2 are merged, the cluster support will be 5.

**Advanced Clustering Strategy**: 
- **Baseline Selection**: Automatically identifies the read with the highest individual support value
- **Position Calculation**: 
  - `start = baseline_start` (exact coordinate)
  - `end = baseline_end` (exact coordinate)
- **Precision**: Cluster boundaries match the exact genomic coordinates of the highest confidence read
- **Window Usage**: Window size determines read grouping but does not expand cluster boundaries
- **Quality-Driven**: Highest confidence read determines cluster position, ensuring accuracy

## Usage

### Command-Line Interface

The main tool is `main.py` which provides a comprehensive interface for testing and clustering:

```bash
python code/main.py --input <input_file> [options]
```

### Required Arguments
- `--input, -i`: Input TSV file with SV candidates

### Optional Arguments
- `--format`: TSV format to use
  - `auto` (default): Auto-detect format based on headers
  - `per-read`: Force per-read candidate format
  - `clustered`: Force clustered summary format
- `--window, -w`: Window size for overlap testing (default: 100)
- `--min-support, -s`: Minimum support value required for clustering (default: 2)
- `--output, -o`: Output TSV file for merged intervals (default: print to console)

### Examples

#### Basic Usage
```bash
# Auto-detect format and test clustering
python code/main.py --input test_data/chr21_SVs.txt

# With custom window size
python code/main.py --input my_svs.tsv --window 200

# Force specific format
python code/main.py --input my_svs.tsv --format per-read
```

#### Clustering with Output
```bash
# Save merged intervals to file
python code/main.py --input my_svs.tsv --output merged_clusters.tsv

# High-quality clusters (min support = 5)
python code/main.py --input my_svs.tsv --min-support 5 --output high_quality.tsv

# Full custom configuration
python code/main.py --input my_svs.tsv --window 500 --min-support 3 --output filtered.tsv
```

#### Short Form
```bash
python code/main.py -i my_svs.tsv -w 200 -s 3 -o clusters.tsv
```

### Programmatic Usage

```python
from code.isosv.struct.intervals import CandidateStructurer, Candidate

# Initialize
structurer = CandidateStructurer()

# Add candidates
candidate = Candidate(chrom="chr21", pos=1000, end=1050, svtype="DEL", 
                     svlen=50, read="read1", strand="+", mapq=60, 
                     cigar="50M", is_primary=True)
structurer.add_candidate(candidate)

# Query
results = structurer.query_region("chr21", "DEL", 1000, 1100)
```

## Output

### Console Output
- **Format detection**: Shows detected input format
- **Parsing summary**: Reports total lines, parsed, and skipped
- **Clustering results**: Shows original vs. merged counts
- **Sample data**: Displays first 10 merged intervals
- **File status**: Confirms output file creation

### File Output
- **TSV format**: Standard tab-separated values
- **Header row**: Column names for easy parsing
- **All clusters**: Complete dataset (not just first 10)
- **Directory creation**: Automatically creates output directories

## Error Handling

- **Missing dependencies**: Clear error messages for missing packages
- **File parsing**: Robust handling of malformed TSV files
- **Directory creation**: Automatic creation of output directories
- **Format validation**: Checks for required columns in input files

## Performance

- **Interval trees**: O(log n) query performance for spatial searches
- **Memory efficient**: Scales well with chromosome-sized datasets
- **Fast clustering**: Efficient overlap detection for merging
- **Baseline selection**: O(n) complexity for finding highest support read in each cluster

## Clustering Methodology

### **Two-Phase Approach**

#### **Phase 1: Read Grouping**
- **Window-based grouping**: Reads within `window_size` bp are grouped together
- **Chromosome + SV type isolation**: Separate groups for each chromosome and SV type combination
- **Spatial proximity**: Ensures related structural variants are grouped appropriately

#### **Phase 2: Baseline Selection & Positioning**
- **Support evaluation**: Each read's individual support value is examined
- **Baseline identification**: Read with highest support becomes the cluster reference point
- **Coordinate assignment**: Cluster start/end positions match baseline read coordinates exactly
- **No artificial expansion**: Window size does not affect final cluster boundaries

### **Example Workflow**
```
Input: 3 reads within 50bp window
- read1: pos=17062902, support=3 (highest)
- read2: pos=17062901, support=2
- read3: pos=17062902, support=1

Process:
1. Group reads (window=50bp)
2. Select baseline: read1 (support=3)
3. Set cluster coordinates: start=17062902, end=17062902
4. Calculate total support: 3+2+1 = 6
5. Output: cluster with exact baseline coordinates
```

### **Quality Assurance**
- **Minimum support threshold**: Filters out low-confidence clusters
- **Baseline validation**: Ensures each cluster has a high-confidence reference read
- **Reproducible positioning**: Same input always produces same cluster coordinates
- **Support aggregation**: Total cluster support reflects cumulative evidence strength

## What This Enables

1. **Precision SV Clustering**: Exact coordinate positioning based on highest confidence reads
2. **Quality-Driven Merging**: Automatic selection of best reference positions for clusters
3. **Spatial Queries**: Quick lookup of SVs in specific genomic regions
4. **Support for Large Datasets**: Interval trees scale well with chromosome-sized data
5. **Advanced Quality Filtering**: Configurable support thresholds with baseline validation
6. **Standardized Output**: TSV format with comprehensive clustering metadata
7. **Reproducible Results**: Consistent clustering strategy regardless of input order

## Next Steps

This implementation provides the foundation for:
- **Advanced clustering algorithms** (density-based, hierarchical) with baseline positioning
- **Consensus calculation** (weighted positions, confidence intervals) based on support values
- **Quality filtering** (MAPQ filtering, strand bias detection) with baseline validation
- **VCF output** generation for standard bioinformatics pipelines
- **Cluster refinement** using additional quality metrics beyond support values

The advanced clustering methodology ensures that each cluster is positioned at the most confident genomic location, while the interval tree infrastructure efficiently handles the spatial organization needed for clustering thousands of SV candidates across the genome.

# IsoSV: RNA-seq Structural Variant Discovery

## Overview

IsoSV is a Python-based pipeline for discovering structural variants (SVs) from RNA-seq alignments. This implementation focuses on **Sections 3 & 4** of the project: **Candidate Structuring (Intervals)** and **Clustering & Deduplication**.

## What You've Built

### 1. Interval Tree Infrastructure (`code/isosv/struct/intervals.py`)

- **`Candidate`** dataclass: Represents per-read SV candidates with genomic coordinates, SV type, and metadata
- **`ClusterSV`** dataclass: Represents clustered/consensus SVs with support counts and confidence intervals
- **`CandidateStructurer`** class: Manages candidates using interval trees for efficient spatial queries
- **`intervaltree` package integration**: Robust interval tree implementation for fast overlap queries

### 2. Testing & Validation (`code/test_intervals.py`)

- **Multi-format support**: Handles per-read candidates and clustered summary formats
- **Auto-detection**: Automatically detects input TSV format based on column headers
- **Interval tree testing**: Loads real SV data and tests tree operations
- **Clustering validation**: Tests interval merging with configurable parameters

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
- **Configurable window size** for overlap detection
- **Minimum support threshold** for filtering low-confidence clusters
- **Automatic directory creation** for output files
- **TSV output** in standardized format

## File Structure

```
step_b_IsoClustering/
├── code/
│   ├── isosv/
│   │   └── struct/
│   │       └── intervals.py          # Core interval tree implementation
│   └── test_intervals.py             # Main testing and clustering tool
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
clusterid	chrom	start	end	type	support	median_sv_len	read_positions
cluster_001	chr22	17062902	17062902	INS	3	322	chr22:17062900:17062902,chr22:17062901:17062903,chr22:17062902:17062904
```

## Usage

### Command-Line Interface

The main tool is `test_intervals.py` which provides a comprehensive interface for testing and clustering:

```bash
python code/test_intervals.py --input <input_file> [options]
```

### Required Arguments
- `--input, -i`: Input TSV file with SV candidates

### Optional Arguments
- `--format`: TSV format to use
  - `auto` (default): Auto-detect format based on headers
  - `per-read`: Force per-read candidate format
  - `clustered`: Force clustered summary format
- `--window, -w`: Window size for overlap testing (default: 100)
- `--min-support, -s`: Minimum support (reads) required for clustering (default: 2)
- `--output, -o`: Output TSV file for merged intervals (default: print to console)

### Examples

#### Basic Usage
```bash
# Auto-detect format and test clustering
python code/test_intervals.py --input test_data/chr21_SVs.txt

# With custom window size
python code/test_intervals.py --input my_svs.tsv --window 200

# Force specific format
python code/test_intervals.py --input my_svs.tsv --format per-read
```

#### Clustering with Output
```bash
# Save merged intervals to file
python code/test_intervals.py --input my_svs.tsv --output merged_clusters.tsv

# High-quality clusters (min support = 5)
python code/test_intervals.py --input my_svs.tsv --min-support 5 --output high_quality.tsv

# Full custom configuration
python code/test_intervals.py --input my_svs.tsv --window 500 --min-support 3 --output filtered.tsv
```

#### Short Form
```bash
python code/test_intervals.py -i my_svs.tsv -w 200 -s 3 -o clusters.tsv
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

## What This Enables

1. **Efficient SV Clustering**: Fast detection of nearby candidates for merging
2. **Spatial Queries**: Quick lookup of SVs in specific genomic regions
3. **Support for Large Datasets**: Interval trees scale well with chromosome-sized data
4. **Quality Filtering**: Configurable support thresholds for confidence
5. **Standardized Output**: TSV format ready for downstream analysis

## Next Steps

This implementation provides the foundation for:
- **Advanced clustering algorithms** (density-based, hierarchical)
- **Consensus calculation** (weighted positions, confidence intervals)
- **Quality filtering** (MAPQ filtering, strand bias detection)
- **VCF output** generation for standard bioinformatics pipelines

The interval tree infrastructure efficiently handles the spatial organization needed for clustering thousands of SV candidates across the genome.

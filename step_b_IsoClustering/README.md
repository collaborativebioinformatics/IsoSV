### 1. Interval Tree Infrastructure (`code/isosv/struct/intervals.py`)

- **`Candidate`** dataclass: Represents per-read SV candidates with genomic coordinates, SV type, and metadata
- **`ClusterSV`** dataclass: Represents clustered/consensus SVs with support counts and confidence intervals
- **`CandidateStructurer`** class: Manages candidates using interval trees for efficient spatial queries
- **`intervaltree` package integration**: Robust interval tree implementation for fast overlap queries

### 2. Testing & Validation (`code/test_intervals.py`)

- **Data converter**: Transforms `test_data/chr21_SVs.txt` → `test_data/chr21_SVs_converted.tsv`
- **Interval tree testing**: Loads real chr21 SV data and tests tree operations
- **Query validation**: Tests region queries, point queries, and overlap detection

## Key Features

### Interval Tree Operations
- **Fast spatial queries** by chromosome and SV type
- **Efficient overlap detection** for clustering candidates
- **Point queries** for specific genomic positions
- **Region queries** for genomic intervals

### Data Handling
- **Flexible input formats**: Handles both raw SV calls and clustered intervals
- **Automatic coordinate normalization**: Ensures consistent chromosome naming
- **INS expansion**: Converts point insertions to intervals using length information

### Query Capabilities
```python
# Region query
candidates = structurer.query_region("chr21", "DEL", 1000000, 2000000)

# Point query  
candidates = structurer.query_point("chr21", "INS", 5000000)

# Overlap detection
nearby = structurer.get_overlapping_candidates(candidate, window=100)
```

## File Structure
chrom start end type support median_sv_len reads genes
21 5227407 5227477 DEL 1 70
21 5276291 5276307 INS 1 2560


## Usage

### 1. Convert Data Format
```bash
python code/test_intervals.py --convert
```
Converts `test_data/chr21_SVs.txt` to clustered format.

### 2. Test Interval Trees
```bash
python code/test_intervals.py
```
Loads converted data and tests interval tree functionality.

### 3. Use in Your Code
```python
from isosv.struct.intervals import CandidateStructurer, Candidate

# Initialize
structurer = CandidateStructurer()

# Add candidates
candidate = Candidate(chrom="chr21", pos=1000, end=1050, ...)
structurer.add_candidate(candidate)

# Query
results = structurer.query_region("chr21", "DEL", 1000, 1100)
```

## Dependencies

```bash
pip install intervaltree pandas
```

## What This Enables

1. **Efficient SV Clustering**: Fast detection of nearby candidates for merging
2. **Spatial Queries**: Quick lookup of SVs in specific genomic regions
3. **Support for Large Datasets**: Interval trees scale well with chromosome-sized data
4. **Foundation for Clustering**: Ready for implementing the clustering algorithms in Section 4

## Next Steps

This implementation provides the foundation for:
- **Clustering algorithms** (single-linkage, density-based)
- **Consensus calculation** (median positions, confidence intervals)
- **Quality filtering** (support thresholds, MAPQ filtering)
- **VCF output** generation

The interval tree infrastructure efficiently handles the spatial organization needed for clustering thousands of SV candidates across the genome.

"""
Created: 08.28.2025 
Authors: Bigy Ambat, Memoona Rasheed
Maintainer: Bigy Ambat

Test script for the interval tree functionality using the intervaltree package.
Supports the new per-read TSV format with detailed SV information.
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), 'isosv'))

from isosv.struct.intervals import Candidate, CandidateStructurer
import argparse
import pandas as pd


def _add_chr_prefix(val: str) -> str:
    """Ensure chromosome names have 'chr' prefix for consistency."""
    s = str(val).strip()
    return s if s.startswith("chr") else f"chr{s}"




def _parse_strand(orientation: str) -> str:
    """Parse orientation field to determine read strand."""
    if pd.isna(orientation) or orientation == '.':
        return "+"
    
    orientation = str(orientation).strip()
    
    # Handle various orientation formats
    if orientation.startswith('F'):
        return "+"
    elif orientation.startswith('R'):
        return "-"
    elif orientation in ['+', '-']:
        return orientation
    else:
        return "+"  # Default to forward


def load_per_read_candidates(tsv_path: str):
    """Load candidates from the new per-read TSV format.
    
    Expected columns:
    chrom, pos_start, pos_end, read_name, sv_type, sv_len, note, read_len, 
    cigar, MAPQ, TLEN, MATE_UNMAPPED, IS_PROPER_PAIR, ORIENTATION, SEQ
    
    Handles comment lines starting with # and skips them.
    """
    try:
        # Read the file and skip comment lines
        with open(tsv_path, 'r') as f:
            lines = f.readlines()
        
        # Filter out comment lines and empty lines
        data_lines = [line.strip() for line in lines if line.strip() and not line.startswith('#')]
        
        if not data_lines:
            print("No data lines found after filtering comments")
            return []
        
        # Parse the header (first non-comment line)
        header = data_lines[0].split('\t')
        print(f"Found header with {len(header)} columns: {header}")
        
        # Check for required columns
        required_cols = ['chrom', 'pos_start', 'pos_end', 'read_name', 'sv_type']
        missing_cols = [col for col in required_cols if col not in header]
        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")
        
        # Create column index mapping
        col_indices = {col: header.index(col) for col in header}
        
        candidates = []
        for i, line in enumerate(data_lines[1:], 1):  # Skip header, start from line 1
            try:
                parts = line.split('\t')
                if len(parts) != len(header):
                    print(f"Warning: Line {i+1} has {len(parts)} parts, expected {len(header)}. Skipping.")
                    continue
                
                # Parse chromosome
                chrom = _add_chr_prefix(parts[col_indices['chrom']])
                
                # Parse positions
                pos_start = int(parts[col_indices['pos_start']])
                pos_end = int(parts[col_indices['pos_end']])
                
                # Parse SV type and length
                sv_type = str(parts[col_indices['sv_type']]).upper()
                sv_len = int(abs(int(parts[col_indices['sv_len']]))) if parts[col_indices['sv_len']] != '.' else 0
                
                # Handle special cases for different SV types
                if sv_type == 'INS':
                    # For insertions, pos_start == pos_end, expand by sv_len
                    if pos_start == pos_end and sv_len > 0:
                        pos_end = pos_start + sv_len
                    elif pos_start == pos_end:
                        pos_end = pos_start + 1  # Minimum span
                
                # Parse read information
                read_name = str(parts[col_indices['read_name']]) if parts[col_indices['read_name']] != '.' else f"read_{i}"
                strand = _parse_strand(parts[col_indices.get('ORIENTATION', 'F/F')])
                
                # Parse MAPQ (handle missing values)
                try:
                    mapq_str = parts[col_indices.get('MAPQ', '60')]
                    mapq = int(mapq_str) if mapq_str != '.' else 30
                except (ValueError, TypeError):
                    mapq = 30
                
                # Parse CIGAR
                cigar = str(parts[col_indices.get('cigar', '')]) if parts[col_indices.get('cigar', '')] != '.' else ""
                
                # Determine if primary alignment
                is_primary = True  # Default assumption
                
                # Parse additional flags
                flags = {
                    'note': str(parts[col_indices.get('note', '')]) if parts[col_indices.get('note', '')] != '.' else "",
                    'read_len': int(parts[col_indices.get('read_len', '0')]) if parts[col_indices.get('read_len', '0')] != '.' else 0,
                    'tlen': int(parts[col_indices.get('TLEN', '0')]) if parts[col_indices.get('TLEN', '0')] != '.' else 0,
                    'mate_unmapped': bool(int(parts[col_indices.get('MATE_UNMAPPED', '0')])) if parts[col_indices.get('MATE_UNMAPPED', '0')] != '.' else False,
                    'is_proper_pair': bool(int(parts[col_indices.get('IS_PROPER_PAIR', '0')])) if parts[col_indices.get('IS_PROPER_PAIR', '0')] != '.' else False,
                    'orientation': str(parts[col_indices.get('ORIENTATION', 'F/F')]) if parts[col_indices.get('ORIENTATION', 'F/F')] != '.' else "F/F",
                    'seq': str(parts[col_indices.get('SEQ', '')]) if parts[col_indices.get('SEQ', '')] != '.' else ""
                }
                
                candidate = Candidate(
                    chrom=chrom,
                    pos=pos_start,
                    end=pos_end,
                    svtype=sv_type,
                    svlen=sv_len,
                    read=read_name,
                    strand=strand,
                    mapq=mapq,
                    prim=is_primary,
                    cigar=cigar,
                    flags=flags
                )
                candidates.append(candidate)
                
            except Exception as e:
                print(f"Warning: Skipping line {i+1} due to error: {e}")
                continue
        
        print(f"Successfully parsed {len(candidates)} candidates")
        return candidates
        
    except Exception as e:
        print(f"Error loading TSV file: {e}")
        return []


def load_clustered_summary_candidates(tsv_path: str):
    """Load candidates from the clustered summary TSV format.
    
    Expected columns:
    chrom, pos, end, svtype, length, support, median_mapq, example_reads, example_ins_seq, genes
    
    This format appears to be pre-clustered with summary statistics.
    """
    try:
        # Read the file and skip comment lines
        with open(tsv_path, 'r') as f:
            lines = f.readlines()
        
        # Filter out comment lines and empty lines
        data_lines = [line.strip() for line in lines if line.strip() and not line.startswith('#')]
        
        if not data_lines:
            print("No data lines found after filtering comments")
            return []
        
        # Parse the header (first non-comment line)
        header = data_lines[0].split('\t')
        print(f"Found clustered summary header with {len(header)} columns: {header}")
        
        # Check for required columns
        required_cols = ['chrom', 'pos', 'end', 'svtype', 'length']
        missing_cols = [col for col in required_cols if col not in header]
        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")
        
        # Create column index mapping
        col_indices = {col: header.index(col) for col in header}
        
        candidates = []
        for i, line in enumerate(data_lines[1:], 1):  # Skip header, start from line 1
            try:
                parts = line.split('\t')
                if len(parts) != len(header):
                    print(f"Warning: Line {i+1} has {len(parts)} parts, expected {len(header)}. Skipping.")
                    continue
                
                # Parse chromosome
                chrom = _add_chr_prefix(parts[col_indices['chrom']])
                
                # Parse positions
                pos_start = int(parts[col_indices['pos']])
                pos_end = int(parts[col_indices['end']])
                
                # Parse SV type and length
                sv_type = str(parts[col_indices['svtype']]).upper()
                sv_len = int(abs(int(parts[col_indices['length']]))) if parts[col_indices['length']] != '.' else 0
                
                # Handle special cases for different SV types
                if sv_type == 'INS':
                    # For insertions, pos_start == pos_end, expand by sv_len
                    if pos_start == pos_end and sv_len > 0:
                        pos_end = pos_start + sv_len
                    elif pos_start == pos_end:
                        pos_end = pos_start + 1  # Minimum span
                
                # Parse read information
                read_name = str(parts[col_indices.get('example_reads', '')]) if parts[col_indices.get('example_reads', '')] != '.' else f"read_{i}"
                strand = "+"  # Default strand for clustered format
                
                # Parse MAPQ (handle missing values)
                try:
                    mapq_str = parts[col_indices.get('median_mapq', '60')]
                    mapq = int(mapq_str) if mapq_str != '.' else 30
                except (ValueError, TypeError):
                    mapq = 30
                
                # Parse CIGAR (not available in this format, use empty string)
                cigar = ""
                
                # Determine if primary alignment
                is_primary = True  # Default assumption
                
                # Parse additional flags
                support = int(parts[col_indices.get('support', '1')]) if parts[col_indices.get('support', '1')] != '.' else 1
                genes = str(parts[col_indices.get('genes', '')]) if parts[col_indices.get('genes', '')] != '.' else ""
                ins_seq = str(parts[col_indices.get('example_ins_seq', '')]) if parts[col_indices.get('example_ins_seq', '')] != '.' else ""
                
                flags = {
                    'support': support,
                    'genes': genes,
                    'ins_seq': ins_seq,
                    'format': 'clustered_summary',
                    'original_pos': pos_start,
                    'original_end': pos_end
                }
                
                candidate = Candidate(
                    chrom=chrom,
                    pos=pos_start,
                    end=pos_end,
                    svtype=sv_type,
                    svlen=sv_len,
                    read=read_name,
                    strand=strand,
                    mapq=mapq,
                    prim=is_primary,
                    cigar=cigar,
                    flags=flags
                )
                candidates.append(candidate)
                
            except Exception as e:
                print(f"Warning: Skipping line {i+1} due to error: {e}")
                continue
        
        print(f"Successfully parsed {len(candidates)} clustered summary candidates")
        return candidates
        
    except Exception as e:
        print(f"Error loading clustered summary TSV file: {e}")
        return []


def auto_detect_format_and_load(tsv_path: str):
    """Automatically detect the TSV format and load candidates accordingly.
    
    Tries to determine if it's per-read format or clustered summary format
    based on the column headers.
    """
    try:
        # Read the first few lines to detect format
        with open(tsv_path, 'r') as f:
            lines = f.readlines()
        
        # Filter out comment lines and empty lines
        data_lines = [line.strip() for line in lines if line.strip() and not line.startswith('#')]
        
        if not data_lines:
            print("No data lines found after filtering comments")
            return []
        
        # Parse the header
        header = data_lines[0].split('\t')
        print(f"Detected header with columns: {header}")
        
        # Check for format-specific columns
        per_read_indicators = ['pos_start', 'pos_end', 'read_name', 'sv_type']
        clustered_indicators = ['pos', 'end', 'svtype', 'length', 'support']
        
        per_read_score = sum(1 for col in per_read_indicators if col in header)
        clustered_score = sum(1 for col in clustered_indicators if col in header)
        
        print(f"Format detection scores - Per-read: {per_read_score}, Clustered: {clustered_score}")
        
        if per_read_score >= clustered_score:
            print("Detected per-read format, using load_per_read_candidates")
            return load_per_read_candidates(tsv_path)
        else:
            print("Detected clustered summary format, using load_clustered_summary_candidates")
            return load_clustered_summary_candidates(tsv_path)
            
    except Exception as e:
        print(f"Error in format detection: {e}")
        return []


def test_interval_tree_with_candidates(candidates):
    """Test the interval tree functionality with provided candidates."""
    if not candidates:
        print("No candidates to test!")
        return
    
    print(f"Testing with {len(candidates)} candidates")
    
    # Test structurer
    print("\nTesting candidate structurer with intervaltree...")
    structurer = CandidateStructurer()
    
    for candidate in candidates:
        structurer.add_candidate(candidate)
    
    print(f"Added candidates to structurer")
    print(f"Total candidates: {len(structurer.get_all_candidates())}")
    
    # Test tree statistics
    print("\nTree statistics:")
    stats = structurer.get_tree_stats()
    for key, stat in stats.items():
        print(f"  {key}: {stat['size']} intervals")
    
    # Test interval queries (driven by the data in the TSV)
    print("\n--- Testing interval queries (data-driven) ---")

    dels = [c for c in candidates if c.svtype == "DEL"]
    ins = [c for c in candidates if c.svtype == "INS"]
    softclips = [c for c in candidates if c.svtype == "SOFTCLIP"]
    splits = [c for c in candidates if c.svtype == "SPLIT"]

    # Region query: around the first DEL in the TSV
    if dels:
        d0 = dels[0]
        w = max(50, min(500, d0.svlen if d0.svlen and d0.svlen > 0 else 100))
        rstart, rend = max(0, d0.pos - w), d0.end + w
        hits = structurer.query_region(d0.chrom, d0.svtype, rstart, rend)
        print(f"Region query around first DEL {d0.chrom}:{d0.pos}-{d0.end} (±{w}): {len(hits)} hits")
        for c in hits[:10]:
            print(f"  {c.chrom}:{c.pos}-{c.end} {c.svtype}")
        if len(hits) > 10:
            print(f"  ... {len(hits)-10} more")
    else:
        print("No DEL entries in the TSV for region query.")

    # Point query: at the first INS position
    if ins:
        i0 = ins[0]
        p_hits = structurer.query_point(i0.chrom, i0.svtype, i0.pos)
        print(f"\nPoint query at first INS {i0.chrom}:{i0.pos}: {len(p_hits)} hits")
        for c in p_hits[:10]:
            print(f"  {c.chrom}:{c.pos}-{c.end} {c.svtype}")
    else:
        print("\nNo INS entries in the TSV for point query.")

    # Overlap query: use a mid DEL (or any candidate if not enough)
    target = dels[len(dels)//2] if len(dels) >= 3 else (dels[0] if dels else (ins[0] if ins else None))
    if target:
        ov_w = 100
        ostart, oend = max(0, target.pos - ov_w), target.end + ov_w
        ov_hits = structurer.query_overlap(target.chrom, target.svtype, ostart, oend)
        print(f"\nOverlap query around {target.chrom}:{target.pos}-{target.end} {target.svtype} (±{ov_w}): {len(ov_hits)} hits")
        for c in ov_hits[:10]:
            print(f"  {c.chrom}:{c.pos}-{c.end} {c.svtype}")
        if len(ov_hits) > 10:
            print(f"  ... {len(ov_hits)-10} more")
    else:
        print("\nNo suitable candidate found for overlap query.")
    
    # Test overlapping candidates with window
    print("\n--- Testing overlapping candidates with window ---")
    test_candidate = candidates[0]
    overlapping = structurer.get_overlapping_candidates(test_candidate, window=100)
    print(f"Candidates overlapping with {test_candidate.chrom}:{test_candidate.pos}-{test_candidate.end} "
          f"{test_candidate.svtype} (window=100): {len(overlapping)}")
    for c in overlapping[:5]:
        print(f"  {c.chrom}:{c.pos}-{c.end} {c.svtype} (read={c.read})")
    if len(overlapping) > 5:
        print(f"  ... {len(overlapping)-5} more")
    
    # Test tree visualization
    print("\n--- Tree visualization ---")
    if dels:
        print(structurer.get_tree_visualization(dels[0].chrom, "DEL"))
    elif ins:
        print(structurer.get_tree_visualization(ins[0].chrom, "INS"))
    elif softclips:
        print(structurer.get_tree_visualization(softclips[0].chrom, "SOFTCLIP"))
    elif splits:
        print(structurer.get_tree_visualization(splits[0].chrom, "SPLIT"))
    else:
        print("No candidates to visualize.")
    
    print("\nTest completed!")








if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="IsoSV Interval Tree Tester - Tests interval trees with multiple TSV formats",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Test with auto-detected format (recommended)
  python test_intervals.py --input my_svs.tsv
  
  # Test with custom window size
  python test_intervals.py --input my_svs.tsv --window 200
  
  # Force per-read format
  python test_intervals.py --input my_svs.tsv --format per-read
  
  # Force clustered format
  python test_intervals.py --input my_svs.tsv --format clustered
        """
    )
    
    parser.add_argument("--input", "-i", type=str, required=True,
                       help="Input TSV file with SV candidates")
    parser.add_argument("--format", choices=['auto', 'per-read', 'clustered'], default='auto',
                       help="TSV format to use (default: auto-detect)")
    parser.add_argument("--window", "-w", type=int, default=100,
                       help="Window size for overlap testing (default: 100)")
    
    args = parser.parse_args()

    print(f"Loading candidates from: {args.input}")
    
    if args.format == 'auto':
        candidates = auto_detect_format_and_load(args.input)
    elif args.format == 'per-read':
        candidates = load_per_read_candidates(args.input)
    elif args.format == 'clustered':
        candidates = load_clustered_summary_candidates(args.input)
    
    if candidates:
        test_interval_tree_with_candidates(candidates)
    else:
        print("Failed to load candidates. Exiting.")
        sys.exit(1)

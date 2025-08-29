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
    """
    try:
        df = pd.read_csv(tsv_path, sep='\t')
        print(f"Loaded TSV with {len(df)} rows and columns: {list(df.columns)}")
        
        # Check for required columns
        required_cols = ['chrom', 'pos_start', 'pos_end', 'read_name', 'sv_type']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")
        
        candidates = []
        for i, row in df.iterrows():
            try:
                # Parse chromosome
                chrom = _add_chr_prefix(row['chrom'])
                
                # Parse positions
                pos_start = int(row['pos_start'])
                pos_end = int(row['pos_end'])
                
                # Parse SV type and length
                sv_type = str(row['sv_type']).upper()
                sv_len = int(abs(row['sv_len'])) if pd.notna(row['sv_len']) else 0
                
                # Handle special cases for different SV types
                if sv_type == 'INS':
                    # For insertions, pos_start == pos_end, expand by sv_len
                    if pos_start == pos_end and sv_len > 0:
                        pos_end = pos_start + sv_len
                    elif pos_start == pos_end:
                        pos_end = pos_start + 1  # Minimum span
                
                # Parse read information
                read_name = str(row['read_name']) if pd.notna(row['read_name']) else f"read_{i+1}"
                strand = _parse_strand(row.get('ORIENTATION', '+'))
                
                # Parse MAPQ (handle missing values)
                try:
                    mapq = int(row['MAPQ']) if pd.notna(row['MAPQ']) and str(row['MAPQ']) != '.' else 30
                except (ValueError, TypeError):
                    mapq = 30
                
                # Parse CIGAR
                cigar = str(row['cigar']) if pd.notna(row['cigar']) else ""
                
                # Determine if primary alignment
                is_primary = True  # Default assumption
                
                # Parse additional flags
                flags = {
                    'note': str(row['note']) if pd.notna(row['note']) else "",
                    'read_len': int(row['read_len']) if pd.notna(row['read_len']) else 0,
                    'tlen': int(row['TLEN']) if pd.notna(row['TLEN']) and str(row['TLEN']) != '.' else 0,
                    'mate_unmapped': bool(row['MATE_UNMAPPED']) if pd.notna(row['MATE_UNMAPPED']) else False,
                    'is_proper_pair': bool(row['IS_PROPER_PAIR']) if pd.notna(row['IS_PROPER_PAIR']) else False,
                    'orientation': str(row['ORIENTATION']) if pd.notna(row['ORIENTATION']) else "+",
                    'seq': str(row['SEQ']) if pd.notna(row['SEQ']) else ""
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
                    cigar=cigar,
                    prim=is_primary,
                    flags=flags
                )
                candidates.append(candidate)
                
            except Exception as e:
                print(f"Warning: Skipping row {i+1} due to error: {e}")
                continue
        
        print(f"Successfully parsed {len(candidates)} candidates")
        return candidates
        
    except Exception as e:
        print(f"Error loading TSV file: {e}")
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
        description="IsoSV Interval Tree Tester - Tests interval trees with per-read SV candidates",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Test with per-read TSV format
  python test_intervals.py --input my_svs.tsv
  
  # Test with custom window size
  python test_intervals.py --input my_svs.tsv --window 200
        """
    )
    
    parser.add_argument("--input", "-i", type=str, required=True,
                       help="Input TSV file with per-read SV candidates")
    parser.add_argument("--window", "-w", type=int, default=100,
                       help="Window size for overlap testing (default: 100)")
    
    args = parser.parse_args()

    print(f"Loading per-read candidates from: {args.input}")
    candidates = load_per_read_candidates(args.input)
    if candidates:
        test_interval_tree_with_candidates(candidates)
    else:
        print("Failed to load candidates. Exiting.")
        sys.exit(1)

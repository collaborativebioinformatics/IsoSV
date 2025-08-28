"""
Test script for the interval tree functionality using the intervaltree package.
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), 'isosv'))

from isosv.struct.intervals import Candidate, CandidateStructurer


def create_test_candidates():
    """Create test candidates for demonstration."""
    candidates = []
    
    # Create some test candidates with overlapping regions
    test_data = [
        # Chrom, pos, end, svtype, svlen, read, strand, mapq, cigar, prim
        ("chr21", 1000, 1050, "DEL", 50, "read1", "+", 30, "50M50D50M", True),
        ("chr21", 1005, 1055, "DEL", 50, "read2", "+", 25, "50M50D50M", True),
        ("chr21", 1010, 1060, "DEL", 50, "read3", "-", 35, "50M50D50M", True),
        ("chr21", 2000, 2100, "INS", 100, "read4", "+", 40, "50M100I50M", True),
        ("chr21", 2005, 2105, "INS", 100, "read5", "+", 30, "50M100I50M", True),
        ("chr21", 3000, 3050, "DEL", 50, "read6", "+", 30, "50M50D50M", True),
        ("chr21", 3000, 3100, "DEL", 100, "read7", "-", 35, "50M100D50M", True),
        ("chr22", 1000, 1050, "DEL", 50, "read8", "+", 30, "50M50D50M", True),
    ]
    
    for i, (chrom, pos, end, svtype, svlen, read, strand, mapq, cigar, prim) in enumerate(test_data):
        candidate = Candidate(
            chrom=chrom,
            pos=pos,
            end=end,
            svtype=svtype,
            svlen=svlen,
            read=read,
            strand=strand,
            mapq=mapq,
            cigar=cigar,
            prim=prim,
            flags={'test': True, 'index': i}
        )
        candidates.append(candidate)
    
    return candidates


def test_interval_tree():
    """Test the interval tree functionality."""
    print("Creating test candidates...")
    candidates = create_test_candidates()
    
    print(f"Created {len(candidates)} test candidates")
    
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
    
    # Test interval queries
    print("\n--- Testing interval queries ---")
    
    # Test region query
    chr21_dels = structurer.query_region("chr21", "DEL", 1000, 1100)
    print(f"chr21 DEL candidates in region 1000-1100: {len(chr21_dels)}")
    for c in chr21_dels:
        print(f"  {c.chrom}:{c.pos}-{c.end} {c.svtype} (read={c.read})")
    
    # Test point query
    chr21_ins_at_2000 = structurer.query_point("chr21", "INS", 2000)
    print(f"\nchr21 INS candidates at position 2000: {len(chr21_ins_at_2000)}")
    for c in chr21_ins_at_2000:
        print(f"  {c.chrom}:{c.pos}-{c.end} {c.svtype} (read={c.read})")
    
    # Test overlap query
    chr21_dels_overlap_3000 = structurer.query_overlap("chr21", "DEL", 3000, 3100)
    print(f"\nchr21 DEL candidates overlapping 3000-3100: {len(chr21_dels_overlap_3000)}")
    for c in chr21_dels_overlap_3000:
        print(f"  {c.chrom}:{c.pos}-{c.end} {c.svtype} (read={c.read})")
    
    # Test overlapping candidates with window
    print("\n--- Testing overlapping candidates with window ---")
    test_candidate = candidates[0]  # chr21:1000-1050 DEL
    overlapping = structurer.get_overlapping_candidates(test_candidate, window=10)
    print(f"Candidates overlapping with {test_candidate.chrom}:{test_candidate.pos}-{test_candidate.end} "
          f"{test_candidate.svtype} (window=10): {len(overlapping)}")
    for c in overlapping:
        print(f"  {c.chrom}:{c.pos}-{c.end} {c.svtype} (read={c.read})")
    
    # Test tree visualization
    print("\n--- Tree visualization ---")
    print(structurer.get_tree_visualization("chr21", "DEL"))
    
    print("\nTest completed!")


if __name__ == "__main__":
    test_interval_tree()

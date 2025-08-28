"""
Test script for the interval tree functionality using the intervaltree package.
Also provides a converter to transform test_data/chr21_SVs.txt into
clustered-interval format with columns:

chrom	start	end	type	support	median_sv_len	reads	genes
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), 'isosv'))

from isosv.struct.intervals import Candidate, CandidateStructurer
import argparse
import pandas as pd


def _add_chr_prefix(val: str) -> str:
    s = str(val).strip()
    return s if s.startswith("chr") else f"chr{s}"


def load_clustered_candidates_from_converted(tsv_path: str):
    """Load candidates from the converted clustered TSV.

    Expects columns: chrom, start, end, type, support, median_sv_len, reads, genes
    Emits ONE Candidate per row (not expanded by reads/support) for testing trees.
    """
    df = pd.read_csv(tsv_path, sep='\t')
    required = ["chrom", "start", "end", "type", "support", "median_sv_len", "reads", "genes"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    candidates = []
    for i, row in df.iterrows():
        chrom = _add_chr_prefix(row["chrom"])  # ensure chr prefix for queries
        start = int(row["start"])
        end = int(row["end"]) if pd.notna(row["end"]) else start
        svtype = str(row["type"]).upper()
        svlen = int(abs(row["median_sv_len"])) if pd.notna(row["median_sv_len"]) else 0

        # Expand INS point to short span using svlen (>=1)
        if svtype == "INS" and start == end:
            end = start + max(1, svlen)

        candidates.append(
            Candidate(
                chrom=chrom,
                pos=start,
                end=end,
                svtype=svtype,
                svlen=svlen,
                read=f"row_{i+1}",
                strand="+",
                mapq=60,
                cigar="",
                prim=True,
                flags={"support": int(row["support"]) if pd.notna(row["support"]) else 1,
                       "genes": str(row["genes"]) if pd.notna(row["genes"]) else ""}
            )
        )
    return candidates


def test_interval_tree():
    """Test the interval tree functionality with converted clustered TSV."""
    project_root = os.path.dirname(os.path.dirname(__file__))
    converted = os.path.join(project_root, "test_data", "chr21_SVs_converted.tsv")
    print(f"Loading candidates from {converted} ...")
    candidates = load_clustered_candidates_from_converted(converted)
    
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
    
    # Test interval queries (driven by the data in the TSV)
    print("\n--- Testing interval queries (data-driven) ---")

    dels = [c for c in candidates if c.svtype == "DEL"]
    ins  = [c for c in candidates if c.svtype == "INS"]

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
    test_candidate = candidates[0]  # chr21:1000-1050 DEL
    overlapping = structurer.get_overlapping_candidates(test_candidate, window=10)
    print(f"Candidates overlapping with {test_candidate.chrom}:{test_candidate.pos}-{test_candidate.end} "
          f"{test_candidate.svtype} (window=10): {len(overlapping)}")
    for c in overlapping:
        print(f"  {c.chrom}:{c.pos}-{c.end} {c.svtype} (read={c.read})")
    
    # Test tree visualization
    print("\n--- Tree visualization ---")
    if dels:
        print(structurer.get_tree_visualization(dels[0].chrom, "DEL"))
    elif ins:
        print(structurer.get_tree_visualization(ins[0].chrom, "INS"))
    else:
        print("No candidates to visualize.")
    
    print("\nTest completed!")


def convert_chr21_to_clustered(in_path: str, out_path: str) -> None:
    """Convert chr21_SVs.txt (chr,start,end,svtype,svlen) to clustered format.

    Output columns: chrom,start,end,type,support,median_sv_len,reads,genes
    - chrom: strip 'chr' prefix if present (e.g., 'chr21' -> '21')
    - type: from svtype
    - support: default 1 (one row per input SV)
    - median_sv_len: abs(svlen) or 1 if missing/NA
    - reads: empty
    - genes: empty
    Rows with missing start/end are skipped.
    """
    df = pd.read_csv(in_path, sep="\t")
    required = ["chr", "start", "end", "svtype", "svlen"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Input is missing required columns: {missing}")

    def strip_chr(x: str) -> str:
        s = str(x)
        return s[3:] if s.startswith("chr") else s

    out = pd.DataFrame()
    out["chrom"] = df["chr"].map(strip_chr)
    # Drop rows with NA start/end
    valid = df["start"].notna() & df["end"].notna()
    df = df[valid].copy()
    out = out.loc[valid].copy()

    out["start"] = df["start"].astype(int)
    out["end"] = df["end"].astype(int)
    out["type"] = df["svtype"].astype(str)
    # Default support 1
    out["support"] = 1
    # median_sv_len: abs of svlen; coerce errors -> NaN then fill with 1
    svlen = pd.to_numeric(df["svlen"], errors="coerce").abs().fillna(1).astype(int)
    out["median_sv_len"] = svlen
    out["reads"] = ""
    out["genes"] = ""

    out.to_csv(out_path, sep="\t", index=False)
    print(f"Wrote converted clustered TSV: {out_path} (rows={len(out)})")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Interval tree tester and chr21 converter")
    parser.add_argument("--convert", action="store_true", help="Convert chr21_SVs.txt to clustered format")
    args = parser.parse_args()

    if args.convert:
        project_root = os.path.dirname(os.path.dirname(__file__))
        in_path = os.path.join(project_root, "test_data", "chr21_SVs.txt")
        out_path = os.path.join(project_root, "test_data", "chr21_SVs_converted.tsv")
        convert_chr21_to_clustered(in_path, out_path)
    else:
        test_interval_tree()

"""
IsoSV Integrated Pipeline

This script serves as the main entry point for the integrated IsoSV pipeline.
It combines three major steps into a single, streamlined workflow:
1.  **Parsing**: Reads a BAM/CRAM file to identify structural variant (SV)
    candidates, such as large indels, soft clips, and split reads.
2.  **Clustering**: Groups the identified SV candidates based on genomic
    proximity and type to form consolidated SV clusters.
3.  **Annotation**: Annotates the SV clusters against a gene reference
    (from a pickled IntervalTree) to provide functional context, such as
    identifying gene fusions or exonic deletions, and outputs a VCF file.

The pipeline is designed to be run from the command line and accepts various
arguments to customize the behavior of each step.
"""
import argparse
import sys
import pysam
import csv
# import all the functions from _utils.py, 
#   which contains the helper functions and dataclasses for parsing, clustering, and annotation
from _utils import (  
    IndelObs, ClipObs, SplitObs, GeneAnnot,
    walk_cigar_indels_and_clips, parse_sa_tag, ref_span_from_cigar,
    cluster_indels, cluster_clips, cluster_splits, load_tx_tree, annotate_candidates, write_to_vcf
)
import pandas as pd

def main():
    ap = argparse.ArgumentParser(description="Integrated pipeline for SV detection, clustering, and annotation.")
    ap.add_argument("in_bam", help="Input BAM/CRAM")
    ap.add_argument("-o", "--out-prefix", required=True, help="Output prefix for TSVs")
    ap.add_argument("--reference", help="Reference FASTA (required for CRAM)")
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-ins", type=int, default=20, help="Minimum insertion length to report")
    ap.add_argument("--min-del", type=int, default=20, help="Minimum deletion length to report")
    ap.add_argument("--min-clip", type=int, default=25, help="Minimum soft/hard clip length to report")
    ap.add_argument("--bp-window", type=int, default=10, help="Breakpoint clustering window (bp)")
    ap.add_argument("--gene-bed", help="BED with gene intervals to annotate outputs")
    ap.add_argument("--annotation-cache", help="Path to the transcript tree cache file for annotation")
    ap.add_argument("--max-readnames", type=int, default=3, help="Max example read names per cluster")
    
    args = ap.parse_args()

    # Step 1: Parsing
    samfile = pysam.AlignmentFile(args.in_bam, "rb")
    indel_obs, clip_obs, split_obs = [], [], []
    for aln in samfile.fetch(until_eof=True):
        if aln.is_unmapped or aln.is_secondary or aln.is_duplicate:
            continue
        if aln.mapping_quality < args.min_mapq:
            continue

        for obs in walk_cigar_indels_and_clips(aln, args.min_ins, args.min_del, args.min_clip):
            if isinstance(obs, IndelObs):
                indel_obs.append(obs)
            else:
                clip_obs.append(obs)
        
        sa = aln.get_tag("SA") if aln.has_tag("SA") else None
        if sa:
            try:
                sa_entries = parse_sa_tag(sa)
                primary_right = aln.reference_end
                for seg in sa_entries:
                    rname2, pos2 = seg["rname"], seg["pos"]
                    split_obs.append(SplitObs(
                        chrom1=aln.reference_name, pos1=primary_right,
                        chrom2=rname2, pos2=pos2,
                        mapq=aln.mapping_quality, rname=aln.query_name,
                        same_chrom=(aln.reference_name == rname2)
                    ))
            except Exception:
                pass # Ignore SA parsing errors
    samfile.close()

    gene_annot = None
    if args.gene_bed:
        gene_annot = GeneAnnot()
        gene_annot.load_bed(args.gene_bed)

    # Step 2: Clustering
    indel_clusters = cluster_indels(indel_obs, args.bp_window, gene_annot=gene_annot, max_readnames=args.max_readnames)
    clip_clusters = cluster_clips(clip_obs, args.bp_window, gene_annot=gene_annot, max_readnames=args.max_readnames)
    split_clusters = cluster_splits(split_obs, args.bp_window, gene_annot=gene_annot, max_readnames=args.max_readnames)

    # Step 3: Annotation
    if args.annotation_cache:
        tx_tree = load_tx_tree(args.annotation_cache)
        
        # Convert clusters to DataFrames for annotation
        indel_df = pd.DataFrame(indel_clusters)
        
        # Ensure correct columns for annotation function
        if not indel_df.empty:
            indel_df["support"] = indel_df["support"].astype(int)
            indel_df = indel_df.rename(columns={"pos": "start", "svtype": "type", "length": "median_sv_len"})
            indel_df["end"] = indel_df.apply(lambda row: row["start"] + row["median_sv_len"] - 1 if row["type"] == "DEL" else row["start"], axis=1)
            
            annotated_indels = annotate_candidates(indel_df, tx_tree)
            annotated_df = pd.DataFrame(annotated_indels, columns=["chr", "start", "stop", "SVTYPE", "ANN_TYPE", "SUPPORT", "SVLEN", "REGION", "TX_ALIAS", "BIOTYPE_GENE", "BIOTYPE_TX"])
            
            vcf_outpath = f"{args.out_prefix}.annotated.vcf"
            write_to_vcf(annotated_df, vcf_outpath)
            print(f"Annotated VCF written to {vcf_outpath}")

    # Step 3: Annotation and Output
    with open(f"{args.out_prefix}.indels.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="	")
        w.writerow(["chrom","pos","end","svtype","length","support","median_mapq","example_reads","example_ins_seq","genes"])
        for c in sorted(indel_clusters, key=lambda x: (x["chrom"], x["pos"], x["svtype"])):
            end = c["pos"] + c["length"] - 1 if c["svtype"] == "DEL" else c["pos"]
            w.writerow([c["chrom"], c["pos"], end, c["svtype"], c["length"], c["support"], c["median_mapq"], c["example_reads"], c["example_ins_seq"], c["genes"]])

    with open(f"{args.out_prefix}.softclips.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="	")
        w.writerow(["chrom","pos","side","support","median_mapq","median_cliplen","example_reads","genes"])
        for c in sorted(clip_clusters, key=lambda x: (x["chrom"], x["pos"], x["side"])):
            w.writerow([c["chrom"], c["pos"], c["side"], c["support"], c["median_mapq"], c["median_cliplen"], c["example_reads"], c["genes"]])

    with open(f"{args.out_prefix}.splits.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="	")
        w.writerow(["chrom1","pos1","chrom2","pos2","support","median_mapq","notes","example_reads","genes_left","genes_right"])
        for c in sorted(split_clusters, key=lambda x: (x["chrom1"], x["pos1"], x["chrom2"], x["pos2"])):
            w.writerow([c["chrom1"], c["pos1"], c["chrom2"], c["pos2"], c["support"], c["median_mapq"], c["notes"], c["example_reads"], c["genes_left"], c["genes_right"]])
            
    print(f"Found {len(indel_clusters)} indel clusters.")
    print(f"Found {len(clip_clusters)} clip clusters.")
    print(f"Found {len(split_clusters)} split clusters.")
    print("Pipeline finished.")

if __name__ == "__main__":
    main()

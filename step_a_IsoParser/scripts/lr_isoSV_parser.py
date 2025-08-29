#!/usr/bin/env python3
"""
Created on Wed August 28 2025
@author: Bharati Jadhav (Github:jadhab01)
Detect SV-sized indels in long-read RNA-seq (BAM/CRAM) using:
  1. CIGAR I(sertion) and D(deltion)
  2. Soft-clipped ends (S/H)
  3. Split reads (SA tag - supplementary alignment)

Coords used in TSV : 1-based POS; DEL end reported (half-open internally).

"""

import argparse
import csv
import statistics
from collections import defaultdict
from dataclasses import dataclass
import pysam
import sys
from bisect import bisect_left, bisect_right

### pysam CIGAR operator codes : #M-match, I-inertion, D-deletion, N-skipped region, S-soft clipping, H-hard clipping, P-Padding, X-mis-match bases, EQ-extented tag, identical bases (X and EQ - hel to detect higher error rate in LR)

M, I, D, N, S, H, P, EQ, X, B = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

###check poly-A tails to avoid inflated SV breakpoints
def is_poly_at(seq: str, min_len: int = 20, frac: float = 0.8) -> bool:
    if not seq or len(seq) < min_len:
        return False
    u = seq.upper()
    a = u.count("A"); t = u.count("T")
    return (a / len(u) >= frac) or (t / len(u) >= frac)

@dataclass
class IndelObs:
    chrom: str
    pos: int            # 1-based
    svtype: str         # "INS" or "DEL"
    length: int
    mapq: int
    rname: str
    seq: str | None = None  # for INS

@dataclass
class ClipObs:
    chrom: str
    pos: int            # 1-based
    side: str           # 'L' or 'R'
    cliplen: int
    mapq: int
    rname: str

@dataclass
class SplitObs:
    chrom1: str
    pos1: int           # 1-based
    chrom2: str
    pos2: int           # 1-based
    mapq: int
    rname: str
    same_chrom: bool
    note: str = ""

###check ref basesd and alignemnt bases based on the CIGAR string
def ref_span_from_cigar(cigar_tuples) -> int:
    span = 0
    for op, length in cigar_tuples or []:
        if op in (M, EQ, X, D, N):
            span += length
    return span



###scan CIGAR to find INDELS
def walk_cigar_indels_and_clips(aln, min_ins, min_del, min_clip):
    chrom = aln.reference_name
    mapq = aln.mapping_quality
    rname = aln.query_name
    try:
        read_seq = aln.query_sequence
    except Exception:
        read_seq = None

    cigar = aln.cigartuples or []
    ref_pos = aln.reference_start  # 0-based
    read_pos = 0

    # Left clip
    if cigar:
        op0, len0 = cigar[0]
        if op0 in (S, H) and len0 >= min_clip:
            clipseq = (read_seq[:len0] if (op0 == S and read_seq is not None) else "")
            if not is_poly_at(clipseq):
                yield ClipObs(chrom, aln.reference_start + 1, 'L', len0, mapq, rname)

    for op, length in cigar:
        if op in (M, EQ, X):
            ref_pos += length
            read_pos += length
        elif op == I:
            ins_seq = (read_seq[read_pos: read_pos + length] if (length >= min_ins and read_seq is not None) else None)
            if length >= min_ins:
                yield IndelObs(chrom, ref_pos + 1, "INS", length, mapq, rname, ins_seq)
            read_pos += length
        elif op == D:
            if length >= min_del:
                yield IndelObs(chrom, ref_pos + 1, "DEL", length, mapq, rname, None)
            ref_pos += length
        elif op == N:
            ref_pos += length
        elif op == S:
            read_pos += length
        elif op in (H, P, B):
            pass

    # Right clip
    if cigar:
        opn, lenn = cigar[-1]
        if opn in (S, H) and lenn >= min_clip:
            clipseq = (read_seq[-lenn:] if (opn == S and read_seq is not None) else "")
            if not is_poly_at(clipseq):
                yield ClipObs(chrom, aln.reference_end, 'R', lenn, mapq, rname)

###scan for supplenatry alignement for split-read evidence (SV breakpoints)
def parse_sa_tag(sa: str):
    segs = []
    if not sa:
        return segs
    entries = [e for e in sa.strip().split(';') if e]
    for e in entries:
        rname, pos, strand, cig, mapq, nm = e.split(',')
        segs.append({
            "rname": rname,
            "pos": int(pos),
            "strand": strand,
            "cigar": cig,
            "mapq": int(mapq),
            "nm": int(nm)
        })
    return segs

###helper fun
def cigarstr_to_tuples(cigar_str: str):
    out = []
    num = ""
    opmap = {'M': M, 'I': I, 'D': D, 'N': N, 'S': S, 'H': H, 'P': P, '=': EQ, 'X': X}
    for ch in cigar_str:
        if ch.isdigit():
            num += ch
        else:
            if ch not in opmap:
                raise ValueError(f"Unknown CIGAR op in SA: {ch}")
            out.append((opmap[ch], int(num)))
            num = ""
    if num:
        raise ValueError("Trailing digits in CIGAR string")
    return out

### for clsutering to merge signals
def within(a, b, w):
    return abs(a - b) <= w

###
def length_similar(a, b, abs_slop=5, frac=0.10):
    return abs(a - b) <= max(abs_slop, int(frac * max(a, b)))

class IntervalIndex:
    """Simple per-chrom interval index with bisect."""
    def __init__(self):
        self.data = {}  # chrom -> list of (start, end, name)

    def add(self, chrom, start, end, name=None):
        self.data.setdefault(chrom, []).append((start, end, name))

    def finalize(self):
        for chrom in self.data:
            self.data[chrom].sort(key=lambda x: (x[0], x[1]))

    def overlaps(self, chrom, pos):
        """pos: 1-based; intervals stored 1-based inclusive."""
        arr = self.data.get(chrom, [])
        if not arr:
            return False, []
        starts = [s for s, _, _ in arr]
        i = bisect_right(starts, pos)
        hits = []
        for k in (i-1, i, i+1):
            if 0 <= k < len(arr):
                s, e, nm = arr[k]
                if s <= pos <= e:
                    hits.append((s, e, nm))
        return (len(hits) > 0), hits

### Gene annotation
class GeneAnnot:
    def __init__(self):
        self.idx = IntervalIndex()

    def load_bed(self, path):
        with open(path) as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 3:
                    continue
                chrom = p[0]
                s1 = int(p[1]) + 1  # 0-based -> 1-based
                e1 = int(p[2])      # inclusive end in 1-based
                name = p[3] if len(p) >= 4 and p[3] else f"gene:{chrom}:{s1}-{e1}"
                self.idx.add(chrom, s1, e1, name)
        self.idx.finalize()

    def genes_at(self, chrom, pos):
        hit, arr = self.idx.overlaps(chrom, pos)
        if not hit:
            return []
        return sorted({nm for _, _, nm in arr if nm})

def main():
    ap = argparse.ArgumentParser(description="SV-Large Indels detection in long-read RNA-seq")
    ap.add_argument("in_bam", help="Input BAM/CRAM ")
    ap.add_argument("-o", "--out-prefix", required=True, help="Output prefix for TSVs")
    ap.add_argument("--reference", help="Reference FASTA (required for CRAM)")
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-ins", type=int, default=20, help="Minimum insertion length to report")
    ap.add_argument("--min-del", type=int, default=20, help="Minimum deletion length to report")
    ap.add_argument("--min-clip", type=int, default=25, help="Minimum soft/hard clip length to report")
    ap.add_argument("--bp-window", type=int, default=10, help="Breakpoint clustering window (bp)")

    # New options
   # ap.add_argument("--splice-bed", help="BED of splice sites (will mask soft-clip clusters within ±--splice-window bp)")
    #ap.add_argument("--splice-window", type=int, default=5, help="Window (bp) around splice sites to mask")
    #ap.add_argument("--vcf", help="Path to write VCF of indel clusters")
    ap.add_argument("--gene-bed", help="BED with gene intervals (name in column 4) to annotate outputs")
    ap.add_argument("--max-readnames", type=int, default=3, help="Max example read names per cluster")
    args = ap.parse_args()

    # Load inputs
    samfile = pysam.AlignmentFile(args.in_bam, "rb", reference_filename=args.reference) if args.reference else pysam.AlignmentFile(args.in_bam, "rb")

    #splice_idx = None
    #if args.splice-bed if False else False:
    #    pass  # placate linter
    #if args.splice_bed:
    #    splice_idx = load_splice_sites_bed(args.splice_bed, args.splice_window)

    gene_annot = None
    if args.gene_bed:
        gene_annot = GeneAnnot()
        gene_annot.load_bed(args.gene_bed)

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

        # Split reads via SA
        sa = aln.get_tag("SA") if aln.has_tag("SA") else None
        if sa:
            try:
                sa_entries = parse_sa_tag(sa)
                primary_left = aln.reference_start + 1
                primary_right = aln.reference_end
                for seg in sa_entries:
                    rname2 = seg["rname"]
                    pos2 = seg["pos"]
                    cigar2 = cigarstr_to_tuples(seg["cigar"])
                    span2 = ref_span_from_cigar(cigar2)
                    split_obs.append(SplitObs(
                        chrom1=aln.reference_name, pos1=primary_right,
                        chrom2=rname2, pos2=pos2,
                        mapq=aln.mapping_quality, rname=aln.query_name,
                        same_chrom=(aln.reference_name == rname2),
                        note="primaryR->SAstart"
                    ))
                    split_obs.append(SplitObs(
                        chrom1=aln.reference_name, pos1=primary_left,
                        chrom2=rname2, pos2=(pos2 + span2 - 1 if span2 > 0 else pos2),
                        mapq=aln.mapping_quality, rname=aln.query_name,
                        same_chrom=(aln.reference_name == rname2),
                        note="primaryL->SAend"
                    ))
            except Exception:
                pass

    samfile.close()

### Indel Clustering
    def cluster_indels(obs_list, bpw):
        clusters = []
        by_chr_type = defaultdict(list)
        for o in obs_list:
            by_chr_type[(o.chrom, o.svtype)].append(o)
        for (chrom, svtype), lst in by_chr_type.items():
            lst.sort(key=lambda x: (x.pos, x.length))
            used = [False]*len(lst)
            for i, oi in enumerate(lst):
                if used[i]: continue
                cluster = [i]
                for j in range(i+1, len(lst)):
                    if used[j]: continue
                    oj = lst[j]
                    if within(oi.pos, oj.pos, bpw) and length_similar(oi.length, oj.length):
                        cluster.append(j); used[j] = True
                pos_med = int(statistics.median([lst[k].pos for k in cluster]))
                len_med = int(statistics.median([lst[k].length for k in cluster]))
                mapqs = [lst[k].mapq for k in cluster]
                rnames = [lst[k].rname for k in cluster][:args.max_readnames]
                seqs = [lst[k].seq for k in cluster if (svtype == "INS" and lst[k].seq)]
                genes = []
                if gene_annot:
                    genes = sorted(set(g for g in gene_annot.genes_at(chrom, pos_med)))
                clusters.append({
                    "chrom": chrom, "svtype": svtype, "pos": pos_med, "length": len_med,
                    "support": len(cluster), "median_mapq": int(statistics.median(mapqs)),
                    "example_reads": ",".join(rnames),
                    "example_ins_seq": (seqs[0] if seqs else ""),
                    "genes": ",".join(genes) if genes else ""
                })
        return clusters

###Soft-clip clustering
    def cluster_clips(obs_list, bpw):
        clusters = []
        by_chr_side = defaultdict(list)
        for o in obs_list:
            by_chr_side[(o.chrom, o.side)].append(o)
        for (chrom, side), lst in by_chr_side.items():
            lst.sort(key=lambda x: x.pos)
            used = [False]*len(lst)
            for i, oi in enumerate(lst):
                if used[i]: continue
                cluster = [i]
                for j in range(i+1, len(lst)):
                    if used[j]: continue
                    oj = lst[j]
                    if within(oi.pos, oj.pos, bpw):
                        cluster.append(j); used[j] = True
                pos_med = int(statistics.median([lst[k].pos for k in cluster]))
                mapqs = [lst[k].mapq for k in cluster]
                clips = [lst[k].cliplen for k in cluster]
                rnames = [lst[k].rname for k in cluster][:args.max_readnames]
                genes = []
                if gene_annot:
                    genes = sorted(set(g for g in gene_annot.genes_at(chrom, pos_med)))
                clusters.append({
                    "chrom": chrom, "pos": pos_med, "side": side,
                    "support": len(cluster),
                    "median_mapq": int(statistics.median(mapqs)),
                    "median_cliplen": int(statistics.median(clips)),
                    "example_reads": ",".join(rnames),
                    "genes": ",".join(genes) if genes else ""
                })
        return clusters

###Split read clustering

    def cluster_splits(obs_list, bpw):
        clusters = []
        by_pair = defaultdict(list)
        for o in obs_list:
            # Keep only inter-chromosomal
            if o.same_chrom:
                continue
            by_pair[(o.chrom1, o.chrom2)].append(o)
        for (c1, c2), lst in by_pair.items():
            lst.sort(key=lambda x: (x.pos1, x.pos2))
            used = [False]*len(lst)
            for i, oi in enumerate(lst):
                if used[i]: continue
                cluster = [i]
                for j in range(i+1, len(lst)):
                    if used[j]: continue
                    oj = lst[j]
                    if within(oi.pos1, oj.pos1, bpw) and within(oi.pos2, oj.pos2, bpw):
                        cluster.append(j); used[j] = True
                pos1_med = int(statistics.median([lst[k].pos1 for k in cluster]))
                pos2_med = int(statistics.median([lst[k].pos2 for k in cluster]))
                mapqs = [lst[k].mapq for k in cluster]
                rnames = [lst[k].rname for k in cluster][:args.max_readnames]
                notes = sorted(set(lst[k].note for k in cluster))
                genes1 = []; genes2 = []
                if gene_annot:
                    genes1 = sorted(set(g for g in gene_annot.genes_at(c1, pos1_med)))
                    genes2 = sorted(set(g for g in gene_annot.genes_at(c2, pos2_med)))
                clusters.append({
                    "chrom1": c1, "pos1": pos1_med, "chrom2": c2, "pos2": pos2_med,
                    "support": len(cluster), "median_mapq": int(statistics.median(mapqs)),
                    "notes": "|".join(notes),
                    "example_reads": ",".join(rnames),
                    "genes_left": ",".join(genes1) if genes1 else "",
                    "genes_right": ",".join(genes2) if genes2 else ""
                })
        return clusters

    indel_clusters = cluster_indels(indel_obs, args.bp_window)
    clip_clusters  = cluster_clips(clip_obs, args.bp_window)
    split_clusters = cluster_splits(split_obs, args.bp_window)  # already inter-chr only


    ### Write Outputs to TSVs
    with open(f"{args.out_prefix}.indels.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["chrom","pos","end","svtype","length","support","median_mapq","example_reads","example_ins_seq","genes"])
        for c in sorted(indel_clusters, key=lambda x: (x["chrom"], x["pos"], x["svtype"])):
            end = c["pos"] + c["length"] - 1 if c["svtype"] == "DEL" else c["pos"]
            w.writerow([c["chrom"], c["pos"], end, c["svtype"], c["length"], c["support"], c["median_mapq"], c["example_reads"], c["example_ins_seq"], c["genes"]])

    with open(f"{args.out_prefix}.softclips.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["chrom","pos","side","support","median_mapq","median_cliplen","example_reads","genes"])
        for c in sorted(clip_clusters, key=lambda x: (x["chrom"], x["pos"], x["side"])):
            w.writerow([c["chrom"], c["pos"], c["side"], c["support"], c["median_mapq"], c["median_cliplen"], c["example_reads"], c["genes"]])

    with open(f"{args.out_prefix}.splits.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["chrom1","pos1","chrom2","pos2","support","median_mapq","notes","example_reads","genes_left","genes_right"])
        for c in sorted(split_clusters, key=lambda x: (x["chrom1"], x["pos1"], x["chrom2"], x["pos2"])):
            w.writerow([c["chrom1"], c["pos1"], c["chrom2"], c["pos2"], c["support"], c["median_mapq"], c["notes"], c["example_reads"], c["genes_left"], c["genes_right"]])

    ###Print Summary
    print(f"[indels] clusters: {len(indel_clusters)}", file=sys.stderr)
    print(f"[softclips] clusters: {len(clip_clusters)}", file=sys.stderr)
    print(f"[splits inter-chr] clusters: {len(split_clusters)}", file=sys.stderr)
if __name__ == "__main__":
    main()

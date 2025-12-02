#!/usr/bin/env python3
"""
@author    : Bharati Jadhav (github: bharatij)
RNA-seq SV caller: isoSV_parser.py

Modes:
  lr  – long-read RNA-seq (PacBio / ONT)
  sr  – paired-end short-read RNA-seq (Illumina)

Outputs (by mode & prefix):
  lr and sr: <prefix>.<mode>.indels.tsv, <prefix>.<mode>.softclips.tsv, <prefix>.<mode>.splits.tsv
  sr sepcific:  <prefix>.<mode>.pairs.tsv, <prefix>.pairs_unmapped.tsv, optional VCF
"""

import argparse
import csv
import statistics
import sys
import re
import time
from collections import defaultdict
from dataclasses import dataclass
from bisect import bisect_right

import pysam

# CIGAR op codes (pysam)
M, I, D, N, S, H, P, EQ, X, B = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

# ---------- Shared utilities ----------
def sample_read_stats(bam, region=None, max_reads=20000):
    # Quick check of BAM characteristics:no. of primary reads examined, median read length ans fracition of paired reads
    
    n = 0
    n_paired = 0
    n_unpaired = 0
    lengths = []

    if region:
        it = bam.fetch(region=region)
    else:
        it = bam.fetch(until_eof=True)

    for aln in it:
        if aln.is_secondary or aln.is_supplementary:
            continue
        if aln.query_length is None:
            continue
        lengths.append(aln.query_length)
        if aln.is_paired:
            n_paired += 1
        else:
            n_unpaired += 1
        n += 1
        if n >= max_reads:
            break

    if not lengths:
        return {"n": 0, "med_len": 0, "frac_paired": 0.0}

    med_len = int(statistics.median(lengths))
    denom = (n_paired + n_unpaired) or 1
    frac_paired = n_paired / denom

    return {
        "n": n,
        "med_len": med_len,
        "frac_paired": frac_paired,
    }


def is_poly_at(seq: str, min_len: int = 20, frac: float = 0.8) -> bool:
    if not seq or len(seq) < min_len:
        return False
    u = seq.upper()
    a = u.count("A")
    t = u.count("T")
    return (a / len(u) >= frac) or (t / len(u) >= frac)


def ref_span_from_cigar(cigs) -> int:
    span = 0
    for op, l in cigs or []:
        if op in (M, EQ, X, D, N):
            span += l
    return span


def cigarstr_to_tuples(s: str):
    m = {'M': M, 'I': I, 'D': D, 'N': N, 'S': S, 'H': H, 'P': P, '=': EQ, 'X': X}
    out = []
    num = ""
    for ch in s:
        if ch.isdigit():
            num += ch
        else:
            if ch not in m:
                raise ValueError(f"Unknown CIGAR op: {ch}")
            out.append((m[ch], int(num)))
            num = ""
    if num:
        raise ValueError("Trailing digits in CIGAR")
    return out


def parse_sa_tag(sa: str):
    segs = []
    if not sa:
        return segs
    for e in [x for x in sa.strip().split(';') if x]:
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


def within(a, b, w):
    return abs(a - b) <= w


def length_similar(a, b, abs_slop=5, frac=0.10):
    return abs(a - b) <= max(abs_slop, int(frac * max(a, b)))


# ---------- Interval index & gene/ splice annotation ----------

class IntervalIndex:
    """Per-chrom 1-based inclusive interval index with cached starts."""
    def __init__(self):
        self.data = {}
        self.starts = {}

    def add(self, chrom, start, end, name=None):
        self.data.setdefault(chrom, []).append((start, end, name))

    def finalize(self):
        for chrom, arr in self.data.items():
            arr.sort(key=lambda x: (x[0], x[1]))
            self.starts[chrom] = [s for s, _, _ in arr]

    def overlaps(self, chrom, pos):
        arr = self.data.get(chrom, [])
        if not arr:
            return False, []
        starts = self.starts.get(chrom, [])
        i = bisect_right(starts, pos)
        hits = []
        for k in (i - 1, i, i + 1):
            if 0 <= k < len(arr):
                s, e, nm = arr[k]
                if s <= pos <= e:
                    hits.append((s, e, nm))
        return (len(hits) > 0), hits


def load_splice_sites_bed(path, window_bp):
    idx = IntervalIndex()
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            s1 = int(p[1]) + 1 - window_bp
            e1 = int(p[2]) + window_bp
            if s1 < 1:
                s1 = 1
            name = p[3] if len(p) >= 4 else None
            idx.add(p[0], s1, e1, name)
    idx.finalize()
    return idx


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
                e1 = int(p[2])      # inclusive
                name = p[3] if len(p) >= 4 and p[3] else f"gene:{chrom}:{s1}-{e1}"
                self.idx.add(chrom, s1, e1, name)
        self.idx.finalize()

    def genes_at(self, chrom, pos):
        hit, arr = self.idx.overlaps(chrom, pos)
        if not hit:
            return []
        return sorted({nm for _, _, nm in arr if nm})


def maj_strand(strs):
    return '+' if sum(s == '+' for s in strs) >= sum(s == '-' for s in strs) else '-'


# ---------- Shared dataclasses ----------

@dataclass
class IndelObs:
    chrom: str
    pos: int
    svtype: str
    length: int
    mapq: int
    rname: str
    seq: str | None = None
    strand: str = "."


@dataclass
class ClipObs:
    chrom: str
    pos: int
    side: str
    cliplen: int
    mapq: int
    rname: str
    strand: str = "."


@dataclass
class SplitObs:
    chrom1: str
    pos1: int
    chrom2: str
    pos2: int
    mapq: int
    rname: str
    same_chrom: bool
    note: str = ""
    strand1: str = "."
    strand2: str = "."


### ---------- LONG-READ  ----------

def walk_cigar_indels_and_clips_lr(aln, min_ins, min_del, min_clip):
    """Long-read CIGAR walker (polyA filter min_len=20)."""
    # --- scan CIGAR to emit indels and end-clips (with strand) ---
    chrom = aln.reference_name
    mapq = aln.mapping_quality
    rname = aln.query_name
    strand = '-' if aln.is_reverse else '+'
    try:
        read_seq = aln.query_sequence
    except Exception:
        read_seq = None

    cigar = aln.cigartuples or []
    ref_pos = aln.reference_start
    read_pos = 0

    # Left clip
    if cigar:
        op0, len0 = cigar[0]
        if op0 in (S, H) and len0 >= min_clip:
            clipseq = (read_seq[:len0] if (op0 == S and read_seq is not None) else "")
            if not is_poly_at(clipseq, min_len=20):
                yield ClipObs(chrom, aln.reference_start + 1, 'L', len0, mapq, rname, strand)

    for op, length in cigar:
        if op in (M, EQ, X):
            ref_pos += length
            read_pos += length
        elif op == I:
            if length >= min_ins:
                ins_seq = (read_seq[read_pos: read_pos + length] if read_seq is not None else None)
                yield IndelObs(chrom, ref_pos + 1, "INS", length, mapq, rname, ins_seq, strand)
            read_pos += length
        elif op == D:
            if length >= min_del:
                yield IndelObs(chrom, ref_pos + 1, "DEL", length, mapq, rname, None, strand)
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
            if not is_poly_at(clipseq, min_len=20):
                yield ClipObs(chrom, aln.reference_end, 'R', lenn, mapq, rname, strand)

# --- Indel clustering ---
def cluster_indels_lr(obs_list, bpw, gene_annot, max_readnames):
    clusters = []
    by_chr_type = defaultdict(list)
    for o in obs_list:
        by_chr_type[(o.chrom, o.svtype)].append(o)
    for (chrom, svtype), lst in by_chr_type.items():
        lst.sort(key=lambda x: (x.pos, x.length))
        used = [False] * len(lst)
        for i, oi in enumerate(lst):
            if used[i]:
                continue
            cluster_idx = [i]
            for j in range(i + 1, len(lst)):
                if used[j]:
                    continue
                oj = lst[j]
                if within(oi.pos, oj.pos, bpw) and length_similar(oi.length, oj.length):
                    cluster_idx.append(j)
                    used[j] = True
            pos_med = int(statistics.median([lst[k].pos for k in cluster_idx]))
            len_med = int(statistics.median([lst[k].length for k in cluster_idx]))
            mapqs = [lst[k].mapq for k in cluster_idx]
            rnames = [lst[k].rname for k in cluster_idx][:max_readnames]
            seqs = [lst[k].seq for k in cluster_idx if (svtype == "INS" and lst[k].seq)]
            strands = [lst[k].strand for k in cluster_idx]
            scons = maj_strand(strands)
            genes = []
            if gene_annot:
                genes = sorted(set(g for g in gene_annot.genes_at(chrom, pos_med)))
            clusters.append({
                "chrom": chrom,
                "svtype": svtype,
                "pos": pos_med,
                "length": len_med,
                "support": len(cluster_idx),
                "median_mapq": int(statistics.median(mapqs)),
                "example_reads": ",".join(rnames),
                "example_ins_seq": (seqs[0] if seqs else ""),
                "genes": ",".join(genes) if genes else "",
                "strand": scons,
            })
    return clusters

# --- Clip read clustering ---
def cluster_clips_lr(obs_list, bpw, gene_annot, max_readnames):
    clusters = []
    by_chr_side = defaultdict(list)
    for o in obs_list:
        by_chr_side[(o.chrom, o.side)].append(o)
    for (chrom, side), lst in by_chr_side.items():
        lst.sort(key=lambda x: x.pos)
        used = [False] * len(lst)
        for i, oi in enumerate(lst):
            if used[i]:
                continue
            cluster_idx = [i]
            for j in range(i + 1, len(lst)):
                if used[j]:
                    continue
                oj = lst[j]
                if within(oi.pos, oj.pos, bpw):
                    cluster_idx.append(j)
                    used[j] = True
            pos_med = int(statistics.median([lst[k].pos for k in cluster_idx]))
            mapqs = [lst[k].mapq for k in cluster_idx]
            clips = [lst[k].cliplen for k in cluster_idx]
            rnames = [lst[k].rname for k in cluster_idx][:max_readnames]
            strands = [lst[k].strand for k in cluster_idx]
            scons = maj_strand(strands)
            genes = []
            if gene_annot:
                genes = sorted(set(g for g in gene_annot.genes_at(chrom, pos_med)))
            clusters.append({
                "chrom": chrom,
                "pos": pos_med,
                "side": side,
                "support": len(cluster_idx),
                "median_mapq": int(statistics.median(mapqs)),
                "median_cliplen": int(statistics.median(clips)),
                "example_reads": ",".join(rnames),
                "genes": ",".join(genes) if genes else "",
                "strand": scons,
            })
    return clusters

# --- split read clustering ---
def cluster_splits_lr(obs_list, bpw, gene_annot, max_readnames):
    clusters = []
    by_pair = defaultdict(list)
    for o in obs_list:
        # keep only inter-chromosomal
        if o.same_chrom:
            continue
        by_pair[(o.chrom1, o.chrom2)].append(o)
    for (c1, c2), lst in by_pair.items():
        lst.sort(key=lambda x: (x.pos1, x.pos2))
        used = [False] * len(lst)
        for i, oi in enumerate(lst):
            if used[i]:
                continue
            cluster_idx = [i]
            for j in range(i + 1, len(lst)):
                if used[j]:
                    continue
                oj = lst[j]
                if within(oi.pos1, oj.pos1, bpw) and within(oi.pos2, oj.pos2, bpw):
                    cluster_idx.append(j)
                    used[j] = True
            pos1_med = int(statistics.median([lst[k].pos1 for k in cluster_idx]))
            pos2_med = int(statistics.median([lst[k].pos2 for k in cluster_idx]))
            mapqs = [lst[k].mapq for k in cluster_idx]
            rnames = [lst[k].rname for k in cluster_idx][:max_readnames]
            notes = sorted(set(lst[k].note for k in cluster_idx))
            spairs = [f"{lst[k].strand1}>{lst[k].strand2}" for k in cluster_idx]
            pair_mode = max(set(spairs), key=spairs.count)
            genes1 = []
            genes2 = []
            if gene_annot:
                genes1 = sorted(set(g for g in gene_annot.genes_at(c1, pos1_med)))
                genes2 = sorted(set(g for g in gene_annot.genes_at(c2, pos2_med)))
            clusters.append({
                "chrom1": c1,
                "pos1": pos1_med,
                "chrom2": c2,
                "pos2": pos2_med,
                "support": len(cluster_idx),
                "median_mapq": int(statistics.median(mapqs)),
                "notes": "|".join(notes),
                "example_reads": ",".join(rnames),
                "genes_left": ",".join(genes1) if genes1 else "",
                "genes_right": ",".join(genes2) if genes2 else "",
                "strand_pair": pair_mode,
            })
    return clusters


def write_lr_tsvs(out_prefix, indel_clusters, clip_clusters, split_clusters):
    with open(f"{out_prefix}.lr.indels.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "chrom", "pos", "end", "svtype", "length", "support",
            "median_mapq", "example_reads", "example_ins_seq", "genes", "strand"
        ])
        for c in sorted(indel_clusters, key=lambda x: (x["chrom"], x["pos"], x["svtype"], x["length"])):
            end = c["pos"] + c["length"] - 1 if c["svtype"] == "DEL" else c["pos"]
            w.writerow([
                c["chrom"], c["pos"], end, c["svtype"], c["length"], c["support"],
                c["median_mapq"], c["example_reads"], c["example_ins_seq"],
                c["genes"], c["strand"]
            ])

    with open(f"{out_prefix}.lr.softclips.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "chrom", "pos", "side", "support", "median_mapq",
            "median_cliplen", "example_reads", "genes", "strand"
        ])
        for c in sorted(clip_clusters, key=lambda x: (x["chrom"], x["pos"], x["side"])):
            w.writerow([
                c["chrom"], c["pos"], c["side"], c["support"], c["median_mapq"],
                c["median_cliplen"], c["example_reads"], c["genes"], c["strand"]
            ])

    with open(f"{out_prefix}.lr.splits.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "chrom1", "pos1", "chrom2", "pos2", "support", "median_mapq",
            "notes", "example_reads", "genes_left", "genes_right", "strand_pair"
        ])
        for c in sorted(split_clusters, key=lambda x: (x["chrom1"], x["pos1"], x["chrom2"], x["pos2"])):
            w.writerow([
                c["chrom1"], c["pos1"], c["chrom2"], c["pos2"], c["support"],
                c["median_mapq"], c["notes"], c["example_reads"],
                c["genes_left"], c["genes_right"], c["strand_pair"]
            ])


def run_lr(args):
    samfile = pysam.AlignmentFile(
        args.in_bam,
        "rb",
        reference_filename=args.reference if args.reference else None,
    )
    # --- MODE CHECK: is really long-read? ---
    stats = sample_read_stats(samfile, max_reads=20000)
    if stats["n"] == 0:
        sys.exit("[error] No primary reads found in BAM/CRAM")
    #frac_paired small and insert size is big in lr
    if stats["frac_paired"] > 0.5 and stats["med_len"] < 400:
        # Looks like paired-end short-read
        sys.exit( f"[ERROR!] BAM looks like paired short-read data "
    f"(median read length ~{stats['med_len']} bp, "
    f"frac_paired={stats['frac_paired']:.2f}), "
    f"but you selected mode 'lr'.")

    samfile.close()

    samfile = pysam.AlignmentFile(
        args.in_bam,
        "rb",
        reference_filename=args.reference if args.reference else None,
    )

    splice_idx = load_splice_sites_bed(args.splice_bed, args.splice_window) if args.splice_bed else None

    gene_annot = None
    if args.gene_bed:
        gene_annot = GeneAnnot()
        gene_annot.load_bed(args.gene_bed)

    indel_obs, clip_obs, split_obs = [], [], []
    if args.region:
        it = samfile.fetch(region=args.region)
    else:
        it = samfile.fetch(until_eof=True)

    for aln in it:
    #for aln in samfile.fetch(until_eof=True):
        if aln.is_unmapped or aln.is_secondary or aln.is_duplicate:
            continue
        if aln.mapping_quality < args.min_mapq:
            continue

        for obs in walk_cigar_indels_and_clips_lr(aln, args.min_ins, args.min_del, args.min_clip):
            if isinstance(obs, IndelObs):
                indel_obs.append(obs)
            else:
                clip_obs.append(obs)

        # split reads via SA, primary-only, inter-chr only
        if (not aln.is_supplementary) and aln.has_tag("SA"):
            sa_entries = []
            try:
                sa_entries = parse_sa_tag(aln.get_tag("SA"))
            except Exception:
                pass

            primary_left = aln.reference_start + 1
            primary_right = aln.reference_end
            strand1 = '-' if aln.is_reverse else '+'
            rname1 = aln.reference_name
            seen_edges = set()

            for seg in sa_entries:
                try:
                    r2 = seg["rname"]
                    p2 = seg["pos"]
                    s2 = seg["strand"]
                    c2 = cigarstr_to_tuples(seg["cigar"])
                    span2 = ref_span_from_cigar(c2)
                except Exception:
                    continue
                if r2 == rname1:
                    continue  # inter-chr only

                # primaryR -> SAstart
                k1 = ("R", r2, p2)
                if k1 not in seen_edges:
                    seen_edges.add(k1)
                    split_obs.append(SplitObs(
                        chrom1=rname1, pos1=primary_right,
                        chrom2=r2, pos2=p2,
                        mapq=aln.mapping_quality, rname=aln.query_name,
                        same_chrom=False, note="primaryR->SAstart",
                        strand1=strand1, strand2=s2
                    ))

                # primaryL -> SAend
                sa_end = p2 + span2 - 1 if span2 > 0 else p2
                k2 = ("L", r2, sa_end)
                if k2 not in seen_edges:
                    seen_edges.add(k2)
                    split_obs.append(SplitObs(
                        chrom1=rname1, pos1=primary_left,
                        chrom2=r2, pos2=sa_end,
                        mapq=aln.mapping_quality, rname=aln.query_name,
                        same_chrom=False, note="primaryL->SAend",
                        strand1=strand1, strand2=s2
                    ))

    samfile.close()

    indel_clusters = cluster_indels_lr(indel_obs, args.bp_window, gene_annot, args.max_readnames)
    clip_clusters = cluster_clips_lr(clip_obs, args.bp_window, gene_annot, args.max_readnames)
    split_clusters = cluster_splits_lr(split_obs, args.bp_window, gene_annot, args.max_readnames)

    # splice masking on clip clusters
    if splice_idx:
        before = len(clip_clusters)
        clip_clusters = [
            c for c in clip_clusters if not splice_idx.overlaps(c["chrom"], c["pos"])[0]
        ]
        removed = before - len(clip_clusters)
        print(f"[lr splice-mask] softclip clusters removed: {removed} (kept {len(clip_clusters)})", file=sys.stderr)

    write_lr_tsvs(args.out_prefix, indel_clusters, clip_clusters, split_clusters)
    print(f"[summary lr] indels={len(indel_clusters)} clips={len(clip_clusters)} splits_interchr={len(split_clusters)}",
          file=sys.stderr)


### ---------- SHORT-READ PAIR-END ----------

@dataclass
class PairObs:
    chrom: str
    pos: int
    chrom2: str
    pos2: int
    tlen: int
    mapq: int
    rname: str
    note: str
    orientation: str


def pair_orientation(aln):
    r = '-' if aln.is_reverse else '+'
    m = '-' if aln.mate_is_reverse else '+'
    if r == '+' and m == '-':
        return 'FR'
    elif r == '-' and m == '+':
        return 'RF'
    elif r == '+' and m == '+':
        return 'FF'
    else:
        return 'RR'


def read_outer_edge(aln):
    return (aln.reference_start + 1) if aln.is_reverse else aln.reference_end


def iter_reads(bam, region: str | None):
    if region:
        yield from bam.fetch(region=region)
    else:
        yield from bam.fetch(until_eof=True)


def robust_insert_stats(bam, region=None, sample_reads=500000):
    vals = []
    n = 0
    for aln in iter_reads(bam, region):
        if n > sample_reads:
            break
        n += 1
        if not aln.is_proper_pair:
            continue
        if aln.is_unmapped or aln.mate_is_unmapped:
            continue
        if aln.reference_id != aln.next_reference_id:
            continue
        t = abs(aln.template_length)
        if 0 < t < 200000:
            vals.append(t)
    if len(vals) < 100:
        return 300, 80
    m = statistics.median(vals)
    mad = statistics.median([abs(x - m) for x in vals]) or 1
    sd = 1.4826 * mad
    return int(m), int(sd)


def walk_cigar_indels_and_clips_pe(aln, min_ins, min_del, min_clip):
    """Short-read CIGAR walker (polyA filter min_len=16)."""
    chrom = aln.reference_name
    mapq = aln.mapping_quality
    rname = aln.query_name
    strand = '-' if aln.is_reverse else '+'
    try:
        seq = aln.query_sequence
    except Exception:
        seq = None
    cigs = aln.cigartuples or []
    ref_pos = aln.reference_start
    read_pos = 0

    # left end clip
    if cigs:
        op0, l0 = cigs[0]
        if op0 in (S, H) and l0 >= min_clip:
            clipseq = (seq[:l0] if (op0 == S and seq) else "")
            if not is_poly_at(clipseq, min_len=16):
                yield ClipObs(chrom, aln.reference_start + 1, 'L', l0, mapq, rname, strand)

    for op, l in cigs:
        if op in (M, EQ, X):
            ref_pos += l
            read_pos += l
        elif op == I:
            if l >= min_ins:
                ins_seq = (seq[read_pos:read_pos + l] if seq else None)
                yield IndelObs(chrom, ref_pos + 1, "INS", l, mapq, rname, ins_seq, strand)
            read_pos += l
        elif op == D:
            if l >= min_del:
                yield IndelObs(chrom, ref_pos + 1, "DEL", l, mapq, rname, None, strand)
            ref_pos += l
        elif op == N:
            ref_pos += l
        elif op == S:
            read_pos += l
        elif op in (H, P, B):
            pass

    # right end clip
    if cigs:
        opn, ln = cigs[-1]
        if opn in (S, H) and ln >= min_clip:
            clipseq = (seq[-ln:] if (opn == S and seq) else "")
            if not is_poly_at(clipseq, min_len=16):
                yield ClipObs(chrom, aln.reference_end, 'R', ln, mapq, rname, strand)


# ---- PE clustering ----

def cluster_clips_linear(obs, bpw, gene, max_readnames, splice_idx=None):
    out = []
    by = defaultdict(list)
    for o in obs:
        by[(o.chrom, o.side)].append(o)

    for (chrom, side), L in by.items():
        L.sort(key=lambda x: x.pos)
        i = 0
        while i < len(L):
            j = i
            while j + 1 < len(L) and (L[j + 1].pos - L[i].pos) <= bpw:
                j += 1
            cluster = L[i:j + 1]
            pos_med = int(statistics.median([x.pos for x in cluster]))
            if splice_idx and splice_idx.overlaps(chrom, pos_med)[0]:
                i = j + 1
                continue
            out.append({
                "chrom": chrom,
                "pos": pos_med,
                "side": side,
                "support": len(cluster),
                "median_mapq": int(statistics.median([x.mapq for x in cluster])),
                "median_cliplen": int(statistics.median([x.cliplen for x in cluster])),
                "example_reads": ",".join([x.rname for x in cluster][:max_readnames]),
                "genes": (",".join(sorted(set(gene.genes_at(chrom, pos_med)))) if gene else ""),
                "strand": (
                    '+' if sum(x.strand == '+' for x in cluster) >= sum(x.strand == '-' for x in cluster)
                    else '-'
                ),
            })
            i = j + 1
    return out


def cluster_indels_nearlinear(obs, bpw, gene, max_readnames, len_bin=10, abs_slop=5, frac=0.10):
    out = []
    by = defaultdict(list)
    for o in obs:
        key = (o.chrom, o.svtype, int(o.length / len_bin))
        by[key].append(o)
    for (chrom, svt, _), L in by.items():
        L.sort(key=lambda x: x.pos)
        i = 0
        while i < len(L):
            j = i
            while j + 1 < len(L):
                nxt = L[j + 1]
                next_ok = ((nxt.pos - L[i].pos) <= bpw and
                           abs(nxt.length - L[i].length) <= max(abs_slop, int(frac * max(nxt.length, L[i].length))))
                if not next_ok:
                    break
                j += 1
            cluster = L[i:j + 1]
            pos_med = int(statistics.median([x.pos for x in cluster]))
            len_med = int(statistics.median([x.length for x in cluster]))
            out.append({
                "chrom": chrom,
                "svtype": svt,
                "pos": pos_med,
                "length": len_med,
                "support": len(cluster),
                "median_mapq": int(statistics.median([x.mapq for x in cluster])),
                "example_reads": ",".join([x.rname for x in cluster][:max_readnames]),
                "example_ins_seq": (next((x.seq for x in cluster if svt == "INS" and x.seq), "")),
                "genes": (",".join(sorted(set(gene.genes_at(chrom, pos_med)))) if gene else ""),
                "strand": (
                    '+' if sum(x.strand == '+' for x in cluster) >= sum(x.strand == '-' for x in cluster)
                    else '-'
                ),
            })
            i = j + 1
    return out


def sa_has_interchr(sa: str, primary_rname: str) -> bool:
    if not sa:
        return False
    for e in sa.split(';'):
        if not e:
            continue
        rname = e.split(',', 1)[0]
        if rname != primary_rname:
            return True
    return False


def cluster_splits_pe(obs, bpw, gene, max_readnames, interchr_only=True):
    by = defaultdict(list)
    for o in obs:
        if interchr_only and o.same_chrom:
            continue
        by[(o.chrom1, o.chrom2)].append(o)
    out = []
    for (c1, c2), L in by.items():
        L.sort(key=lambda x: (x.pos1, x.pos2))
        used = [False] * len(L)
        for i, oi in enumerate(L):
            if used[i]:
                continue
            idx = [i]
            for j in range(i + 1, len(L)):
                if used[j]:
                    continue
                if within(oi.pos1, L[j].pos1, bpw) and within(oi.pos2, L[j].pos2, bpw):
                    idx.append(j)
                    used[j] = True
            pos1 = int(statistics.median([L[k].pos1 for k in idx]))
            pos2 = int(statistics.median([L[k].pos2 for k in idx]))
            mapq = int(statistics.median([L[k].mapq for k in idx]))
            reads = ",".join([L[k].rname for k in idx][:max_readnames])
            notes = "|".join(sorted(set(L[k].note for k in idx)))
            spairs = [f"{L[k].strand1}>{L[k].strand2}" for k in idx]
            strand_pair = max(set(spairs), key=spairs.count)
            g1 = g2 = ""
            if gene:
                g1 = ",".join(sorted(set(gene.genes_at(c1, pos1))))
                g2 = ",".join(sorted(set(gene.genes_at(c2, pos2))))
            out.append({
                "chrom1": c1,
                "pos1": pos1,
                "chrom2": c2,
                "pos2": pos2,
                "support": len(idx),
                "median_mapq": mapq,
                "notes": notes,
                "example_reads": reads,
                "genes_left": g1,
                "genes_right": g2,
                "strand_pair": strand_pair,
            })
    return out


def summarize_pair_cluster(sub, c1, c2, bpw, gene, max_readnames):
    pos1 = int(statistics.median([x.pos for x in sub]))
    pos2 = int(statistics.median([x.pos2 for x in sub]))
    med_mapq = int(statistics.median([x.mapq for x in sub])) if sub else 0
    tlens = [x.tlen for x in sub if x.tlen > 0]
    med_tlen = int(statistics.median(tlens)) if tlens else 0
    reads = ",".join(list({x.rname for x in sub})[:max_readnames])
    notes = "|".join(sorted(set(x.note for x in sub)))
    orients = ",".join(sorted(set(x.orientation for x in sub)))
    g1 = g2 = ""
    if gene:
        g1 = ",".join(sorted(set(gene.genes_at(c1, pos1))))
        g2 = ",".join(sorted(set(gene.genes_at(c2, pos2))))
    return {
        "chrom1": c1,
        "pos1": pos1,
        "chrom2": c2,
        "pos2": pos2,
        "support": len(sub),
        "median_mapq": med_mapq,
        "median_tlen": med_tlen,
        "notes": notes,
        "orientations": orients,
        "example_reads": reads,
        "genes_left": g1,
        "genes_right": g2,
    }


def cluster_pairs_grid_refine(obs, bpw, gene, max_readnames,
                              refine_min=20, include_neighbors=True):
    from collections import defaultdict as _dd
    grid = _dd(list)
    for o in obs:
        key = (o.chrom, o.pos // bpw, o.chrom2, o.pos2 // bpw)
        grid[key].append(o)

    heavy = [k for k, v in grid.items() if len(v) >= refine_min]
    heavy.sort(key=lambda k: -len(grid[k]))

    assigned = set()
    clusters = []
    processed = set()

    for key in heavy:
        if key in processed:
            continue
        c1, b1, c2, b2 = key
        neigh = {(c1, b1, c2, b2)}
        if include_neighbors:
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    k2 = (c1, b1 + dx, c2, b2 + dy)
                    if k2 in grid:
                        neigh.add(k2)

        local = []
        for kk in neigh:
            for o in grid[kk]:
                if id(o) not in assigned:
                    local.append(o)
        if not local:
            processed.update(neigh)
            continue

        local.sort(key=lambda x: (x.pos, x.pos2))
        used = [False] * len(local)
        for i, oi in enumerate(local):
            if used[i]:
                continue
            idx = [i]
            for j in range(i + 1, len(local)):
                if used[j]:
                    continue
                if within(oi.pos, local[j].pos, bpw) and within(oi.pos2, local[j].pos2, bpw):
                    idx.append(j)
                    used[j] = True
            sub = [local[k] for k in idx]
            for o in sub:
                assigned.add(id(o))
            clusters.append(summarize_pair_cluster(sub, c1, c2, bpw, gene, max_readnames))

        processed.update(neigh)

    # aggregate remaining per cell
    for (c1, b1, c2, b2), L in grid.items():
        rest = [o for o in L if id(o) not in assigned]
        if not rest:
            continue
        pos1 = b1 * bpw + bpw // 2
        pos2 = b2 * bpw + bpw // 2
        med_mapq = int(statistics.median([o.mapq for o in rest]))
        tlens = [o.tlen for o in rest if o.tlen > 0]
        med_tlen = int(statistics.median(tlens)) if tlens else 0
        reads = ",".join(list({o.rname for o in rest})[:max_readnames])
        notes = "|".join(sorted(set(o.note for o in rest)))
        orients = ",".join(sorted(set(o.orientation for o in rest)))
        g1 = g2 = ""
        if gene:
            g1 = ",".join(sorted(set(gene.genes_at(c1, pos1))))
            g2 = ",".join(sorted(set(gene.genes_at(c2, pos2))))
        clusters.append({
            "chrom1": c1,
            "pos1": pos1,
            "chrom2": c2,
            "pos2": pos2,
            "support": len(rest),
            "median_mapq": med_mapq,
            "median_tlen": med_tlen,
            "notes": notes,
            "orientations": orients,
            "example_reads": reads,
            "genes_left": g1,
            "genes_right": g2,
        })
    return clusters


def cluster_unmapped_pairs_1d(obs, bpw, gene, max_readnames):
    out = []
    by = defaultdict(list)
    for o in obs:
        by[o.chrom].append(o)
    for chrom, L in by.items():
        L.sort(key=lambda x: x.pos)
        i = 0
        while i < len(L):
            j = i
            while j + 1 < len(L) and (L[j + 1].pos - L[i].pos) <= bpw:
                j += 1
            sub = L[i:j + 1]
            pos = int(statistics.median([x.pos for x in sub]))
            med_mapq = int(statistics.median([x.mapq for x in sub]))
            reads = ",".join(list({x.rname for x in sub})[:max_readnames])
            notes = "|".join(sorted(set(x.note for x in sub)))
            orients = ",".join(sorted(set(x.orientation for x in sub)))
            g = ""
            if gene:
                g = ",".join(sorted(set(gene.genes_at(chrom, pos))))
            out.append({
                "chrom": chrom,
                "pos": pos,
                "support": len(sub),
                "median_mapq": med_mapq,
                "notes": notes,
                "orientations": orients,
                "example_reads": reads,
                "genes": g,
            })
            i = j + 1
    return out


def write_pe_tsvs(out_prefix, indel_clusters, clip_clusters,
                  split_clusters, pair_clusters, unmapped_clusters):
    with open(f"{out_prefix}.sr.indels.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "chrom", "pos", "end", "svtype", "length", "support",
            "median_mapq", "example_reads", "example_ins_seq", "genes", "strand"
        ])
        for c in sorted(indel_clusters, key=lambda x: (x["chrom"], x["pos"], x["svtype"], x["length"])):
            end = c["pos"] + c["length"] - 1 if c["svtype"] == "DEL" else c["pos"]
            w.writerow([
                c["chrom"], c["pos"], end, c["svtype"], c["length"], c["support"],
                c["median_mapq"], c["example_reads"], c["example_ins_seq"],
                c["genes"], c["strand"]
            ])

    with open(f"{out_prefix}.sr.softclips.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "chrom", "pos", "side", "support", "median_mapq",
            "median_cliplen", "example_reads", "genes", "strand"
        ])
        for c in sorted(clip_clusters, key=lambda x: (x["chrom"], x["pos"], x["side"])):
            w.writerow([
                c["chrom"], c["pos"], c["side"], c["support"],
                c["median_mapq"], c["median_cliplen"],
                c["example_reads"], c["genes"], c["strand"]
            ])

    with open(f"{out_prefix}.sr.splits.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "chrom1", "pos1", "chrom2", "pos2", "support", "median_mapq",
            "notes", "example_reads", "genes_left", "genes_right", "strand_pair"
        ])
        for c in sorted(split_clusters, key=lambda x: (x["chrom1"], x["pos1"], x["chrom2"], x["pos2"])):
            w.writerow([
                c["chrom1"], c["pos1"], c["chrom2"], c["pos2"], c["support"],
                c["median_mapq"], c["notes"], c["example_reads"],
                c["genes_left"], c["genes_right"], c["strand_pair"]
            ])

    with open(f"{out_prefix}.sr.pairs.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "chrom1", "pos1", "chrom2", "pos2", "support", "median_mapq",
            "median_tlen", "notes", "orientations", "example_reads",
            "genes_left", "genes_right"
        ])
        for c in sorted(pair_clusters, key=lambda x: (x["chrom1"], x["pos1"], x["chrom2"], x["pos2"])):
            w.writerow([
                c["chrom1"], c["pos1"], c["chrom2"], c["pos2"], c["support"],
                c["median_mapq"], c["median_tlen"], c["notes"],
                c["orientations"], c["example_reads"],
                c["genes_left"], c["genes_right"]
            ])

    with open(f"{out_prefix}.sr.pairs_unmapped.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "chrom", "pos", "support", "median_mapq",
            "notes", "orientations", "example_reads", "genes"
        ])
        for c in sorted(unmapped_clusters, key=lambda x: (x["chrom"], x["pos"])):
            w.writerow([
                c["chrom"], c["pos"], c["support"], c["median_mapq"],
                c["notes"], c["orientations"],
                c["example_reads"], c["genes"]
            ])


def run_pe(args):
    print("Starting BAM scan for preliminary checks .....", file=sys.stderr)
    bam = pysam.AlignmentFile(
        args.in_bam, "rb",
        reference_filename=args.reference if args.reference else None,
        threads=args.threads,
    )

    if not bam.has_index():
        sys.exit("[error] BAM/CRAM is not indexed (.bai/.crai missing). Run: samtools index <file>")

    if args.region:
        contig = args.region.split(':', 1)[0]
        if contig not in set(bam.references):
            ex = ", ".join(list(bam.references[:5]))
            sys.exit(f"[error] Contig '{contig}' not found in BAM header. Example contigs: {ex}")

    # --- MODE CHECK: is paired short-read? ---
    stats = sample_read_stats(bam, region=args.region, max_reads=20000)
    if stats["n"] == 0:
        sys.exit("[error] No primary reads found in BAM/CRAM")

    if stats["frac_paired"] < 0.5 or stats["med_len"] > 600:
        # Looks more like long-read or single-end
        sys.exit( f"[ERROR!] BAM does not look like paired short-read data "
            f"(median read length ~{stats['med_len']} bp, "
            f"frac_paired={stats['frac_paired']:.2f}); "
            f"but you selected mode 'sr'.")

    print("Generating splice index if bed provided .....", file=sys.stderr)
    splice_idx = load_splice_sites_bed(args.splice_bed, args.splice_window) if args.splice_bed else None
    print("Generating gene annotation index if bed provided .....", file=sys.stderr)
    gene = None
    if args.gene_bed:
        gene = GeneAnnot()
        gene.load_bed(args.gene_bed)

    if args.tlen_thresh is not None:
        del_like_thresh = int(args.tlen_thresh)
        print(f"[pe] TLEN threshold (override): {del_like_thresh}", file=sys.stderr)
    else:
        mean_is, sd_is = robust_insert_stats(bam, region=args.region)
        del_like_thresh = int(mean_is + args.pair_z * max(sd_is, 10))
        print(f"[pe] insert mean={mean_is}, sd={sd_is}, DEL-like TLEN≥{del_like_thresh}", file=sys.stderr)
        bam.close()
        print("Preliminary BAM scan finished!", file=sys.stderr)
        bam = pysam.AlignmentFile(
            args.in_bam, "rb",
            reference_filename=args.reference if args.reference else None,
            threads=args.threads,
        )

    sa_allow_pat = re.compile(args.sa_allow_contigs)

    indel_obs = []
    clip_obs = []
    split_obs = []
    pair_obs_main = []
    pair_obs_unmapped = []
    seen_tpl_cats = set()

    count = 0
    t0 = time.time()
    print("Starting BAM scan for SV detection .....", file=sys.stderr)

    for aln in iter_reads(bam, args.region):
        count += 1
        if args.progress_reads and (count % args.progress_reads == 0):
            elapsed = time.time() - t0
            print(f"[progress] {count:,} reads processed in {args.region or 'ALL'}  ({elapsed:.1f}s)",
                  file=sys.stderr, flush=True)

        if aln.is_unmapped or aln.is_secondary or aln.is_duplicate:
            continue
        if aln.mapping_quality < args.min_mapq:
            continue

        # --- pair evidence ---
        if aln.is_paired:
            qn = aln.query_name
            if aln.mate_is_unmapped:
                cat = ("unmapped")
                if (qn, cat) not in seen_tpl_cats:
                    seen_tpl_cats.add((qn, cat))
                    note = "one_end_unmapped"
                    cigs = aln.cigartuples or []
                    lclip = cigs[0][1] if cigs and cigs[0][0] in (S, H) else 0
                    rclip = cigs[-1][1] if cigs and cigs[-1][0] in (S, H) else 0
                    if max(lclip, rclip) >= args.min_clip:
                        note += "|clip_support"
                    edge = read_outer_edge(aln)
                    o = PairObs(aln.reference_name, edge, aln.reference_name, edge,
                                0, aln.mapping_quality, qn, note, pair_orientation(aln))
                    pair_obs_unmapped.append(o)
                    if args.pair_include_unmapped_in_main:
                        pair_obs_main.append(o)
            else:
                if aln.reference_id != aln.next_reference_id:
                    cat = ("interchr")
                    if (qn, cat) not in seen_tpl_cats:
                        seen_tpl_cats.add((qn, cat))
                        pair_obs_main.append(PairObs(
                            aln.reference_name, read_outer_edge(aln),
                            bam.get_reference_name(aln.next_reference_id), aln.next_reference_start + 1,
                            abs(aln.template_length), aln.mapping_quality, qn,
                            "interchr", pair_orientation(aln)))
                else:
                    orient = pair_orientation(aln)
                    tlen = abs(aln.template_length)
                    left = min(aln.reference_start, aln.next_reference_start) + 1
                    right = max(aln.reference_end, aln.next_reference_start + 1)
                    if tlen >= del_like_thresh:
                        cat = ("tlen_large")
                        if (qn, cat) not in seen_tpl_cats:
                            seen_tpl_cats.add((qn, cat))
                            pair_obs_main.append(PairObs(
                                aln.reference_name, left, aln.reference_name, right,
                                tlen, aln.mapping_quality, qn, "DEL_like_large_TLEN", orient))
                    if orient in ("RF", "FF", "RR"):
                        cat = (f"orient_{orient}")
                        if (qn, cat) not in seen_tpl_cats:
                            seen_tpl_cats.add((qn, cat))
                            pair_obs_main.append(PairObs(
                                aln.reference_name, left, aln.reference_name, right,
                                tlen, aln.mapping_quality, qn, f"orientation_{orient}", orient))

        # --- pre-check for local evidence ---
        cigs = aln.cigartuples or []
        has_big_endclip = (
            (cigs and (cigs[0][0] in (S, H) and cigs[0][1] >= args.min_clip)) or
            (cigs and (cigs[-1][0] in (S, H) and cigs[-1][1] >= args.min_clip))
        )
        has_big_indel = any(
            (op == I and l >= args.min_ins) or (op == D and l >= args.min_del)
            for op, l in cigs
        )
        has_sa = (not aln.is_supplementary) and aln.has_tag("SA")
        if has_sa and args.interchr_only_splits:
            sa_str = aln.get_tag("SA")
            if not sa_has_interchr(sa_str, aln.reference_name):
                has_sa = False

        if not (has_big_endclip or has_big_indel or has_sa):
            continue

        # --- local evidence ---
        for o in walk_cigar_indels_and_clips_pe(aln, args.min_ins, args.min_del, args.min_clip):
            if isinstance(o, IndelObs):
                indel_obs.append(o)
            else:
                clip_obs.append(o)

        if has_sa:
            try:
                sa_entries = parse_sa_tag(aln.get_tag("SA"))
            except Exception:
                sa_entries = []
            sa_entries = [seg for seg in sa_entries
                          if seg["mapq"] >= args.sa_min_mapq and sa_allow_pat.match(seg["rname"])]
            sa_entries.sort(key=lambda s: (-s["mapq"], s["nm"]))
            sa_entries = sa_entries[:max(0, args.sa_max_per_read)]

            primary_left = aln.reference_start + 1
            primary_right = aln.reference_end
            strand1 = '-' if aln.is_reverse else '+'
            seen_edges = set()
            for seg in sa_entries:
                try:
                    r2 = seg["rname"]
                    p2 = seg["pos"]
                    s2 = seg["strand"]
                    c2 = cigarstr_to_tuples(seg["cigar"])
                    span2 = ref_span_from_cigar(c2)
                except Exception:
                    continue
                k1 = ("R", r2, p2)
                if k1 not in seen_edges:
                    seen_edges.add(k1)
                    split_obs.append(SplitObs(
                        chrom1=aln.reference_name, pos1=primary_right,
                        chrom2=r2, pos2=p2,
                        mapq=aln.mapping_quality, rname=aln.query_name,
                        same_chrom=(aln.reference_name == r2),
                        note="primaryR->SAstart",
                        strand1=strand1, strand2=s2
                    ))
                sa_end = p2 + span2 - 1 if span2 > 0 else p2
                k2 = ("L", r2, sa_end)
                if k2 not in seen_edges:
                    seen_edges.add(k2)
                    split_obs.append(SplitObs(
                        chrom1=aln.reference_name, pos1=primary_left,
                        chrom2=r2, pos2=sa_end,
                        mapq=aln.mapping_quality, rname=aln.query_name,
                        same_chrom=(aln.reference_name == r2),
                        note="primaryL->SAend",
                        strand1=strand1, strand2=s2
                    ))

    bam.close()
    print("Finish main BAM scan!", file=sys.stderr)

    if args.debug_pairs:
        by = defaultdict(int)
        for o in pair_obs_main:
            by[(o.chrom, o.chrom2)] += 1
        top = sorted(by.items(), key=lambda x: -x[1])[:5]
        print(f"[debug] pair_obs_main={len(pair_obs_main)} buckets(top): {top}", file=sys.stderr)
        print(f"[debug] pair_obs_unmapped={len(pair_obs_unmapped)}", file=sys.stderr)

    print("Performing Indel Clustering .....", file=sys.stderr)
    indel_clusters = cluster_indels_nearlinear(indel_obs, args.bp_window, gene, args.max_readnames)
    print("Performing Clip Clustering .....", file=sys.stderr)
    clip_clusters = cluster_clips_linear(clip_obs, args.bp_window, gene, args.max_readnames, splice_idx)
    print("Performing Split read Clustering .....", file=sys.stderr)
    split_clusters = cluster_splits_pe(split_obs, args.bp_window, gene, args.max_readnames, args.interchr_only_splits)
    print("Performing pair read Clustering .....", file=sys.stderr)
    pair_clusters = cluster_pairs_grid_refine(
        pair_obs_main, args.bp_window, gene, args.max_readnames,
        refine_min=args.pair_refine_min, include_neighbors=True
    )
    unmapped_clusters = cluster_unmapped_pairs_1d(
        pair_obs_unmapped, args.bp_window, gene, args.max_readnames
    )

    print("Writing to output files .....", file=sys.stderr)
    write_pe_tsvs(args.out_prefix, indel_clusters, clip_clusters,
                  split_clusters, pair_clusters, unmapped_clusters)

    # Optional VCF (indels only)
    if args.vcf:
        with open(args.vcf, "w") as v:
            v.write("##fileformat=VCFv4.3\n")
            v.write("##ALT=<ID=DEL,Description=\"Deletion\">\n")
            v.write("##ALT=<ID=INS,Description=\"Insertion\">\n")
            v.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"SV type\">\n")
            v.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"SV length (DEL negative)\">\n")
            v.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End coordinate\">\n")
            v.write("##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"CI around POS\">\n")
            v.write("##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"CI around END\">\n")
            v.write("##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description=\"Cluster support\">\n")
            v.write("##INFO=<ID=MEDIAN_MAPQ,Number=1,Type=Integer,Description=\"Median MAPQ\">\n")
            v.write("##INFO=<ID=RNAMES,Number=.,Type=String,Description=\"Example read names\">\n")
            v.write("##INFO=<ID=GENE,Number=.,Type=String,Description=\"Genes overlapping POS\">\n")
            v.write("##INFO=<ID=STRAND,Number=1,Type=String,Description=\"Consensus strand of cluster\">\n")
            v.write("##source=rnaseq_sv.py:pe\n")
            v.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for idx, c in enumerate(sorted(indel_clusters, key=lambda x: (x["chrom"], x["pos"], x["svtype"], x["length"]))):
                chrom = c["chrom"]
                pos = c["pos"]
                L = c["length"]
                svt = c["svtype"]
                if svt == "DEL":
                    alt = "<DEL>"
                    end = pos + L - 1
                    svlen = -abs(L)
                else:
                    alt = "<INS>"
                    end = pos
                    svlen = abs(L)
                info = [
                    f"SVTYPE={svt}",
                    f"SVLEN={svlen}",
                    f"END={end}",
                    f"CIPOS=-{args.bp_window},{args.bp_window}",
                    f"CIEND=-{args.bp_window},{args.bp_window}",
                    f"SUPPORT={c['support']}",
                    f"MEDIAN_MAPQ={c['median_mapq']}",
                ]
                if c.get("example_reads"):
                    info.append(f"RNAMES={c['example_reads']}")
                if c.get("genes"):
                    info.append(f"GENE={c['genes']}")
                info.append(f"STRAND={c.get('strand', '.')}")
                rid = f"{svt}_{chrom}_{pos}_{L}_{idx+1}"
                v.write(f"{chrom}\t{pos}\t{rid}\tN\t{alt}\t.\tPASS\t{';'.join(info)}\n")

    print(f"[summary pe] region={args.region or 'ALL'} "
          f"indels={len(indel_clusters)} clips={len(clip_clusters)} "
          f"splits={len(split_clusters)} pairs={len(pair_clusters)} "
          f"unmapped_pairs={len(unmapped_clusters)}", file=sys.stderr)

def build_parser():
    p = argparse.ArgumentParser(description="RNA-seq SV caller (long-read (PacBio/ONT) + paired-end short-read (Illumina))")
    sub = p.add_subparsers(dest="mode", required=True)

    # Long-read subcommand
    pl = sub.add_parser("lr", help="Long-read RNA-seq SV caller")
    pl.add_argument("in_bam", help="Input BAM/CRAM")
    pl.add_argument("-o", "--out-prefix", required=True, help="Output prefix")
    pl.add_argument("--reference", help="Reference FASTA (for CRAM)")
    pl.add_argument("--region", help="Region like chr1 or chr1:100000-200000")
    pl.add_argument("--min-mapq", type=int, default=20, help="Minimum mapq to report, default=20")
    pl.add_argument("--min-ins", type=int, default=30, help="Minimum insertion length to report, default=30")
    pl.add_argument("--min-del", type=int, default=30, help="Minimum deletion length to report, default=30")
    pl.add_argument("--min-clip", type=int, default=25, help="Minimum soft/hard clip length to report, default=25")
    pl.add_argument("--bp-window", type=int, default=10, help="Breakpoint clustering window (bp), default=10")
    pl.add_argument("--splice-bed", help="Splice-site BED for soft-clip masking" )
    pl.add_argument("--splice-window", type=int, default=5, help="Window (bp) around splice sites to mask, default=5")
    pl.add_argument("--gene-bed", help="Gene BED3/4+ (name in col4) to annotate")
    pl.add_argument("--max-readnames", type=int, default=3, help="Max example read names per cluster, default=3")
    # Short-read PE subcommand
    pp = sub.add_parser("sr", help="Paired-end short-read RNA-seq SV caller")
    pp.add_argument("in_bam", help="Input BAM/CRAM")
    pp.add_argument("-o", "--out-prefix", required=True, help="Output prefix")
    pp.add_argument("--reference", help="Reference FASTA (for CRAM)")
    pp.add_argument("--region", help="Region like chr1 or chr1:100000-200000")
    pp.add_argument("--threads", type=int, default=4, help="BGZF decompression threads")
    pp.add_argument("--min-mapq", type=int, default=20, help="Minimum mapq to report, default=20")
    pp.add_argument("--min-ins", type=int, default=30, help="Minimum insertion length to report, default=30")
    pp.add_argument("--min-del", type=int, default=30, help="Minimum deletion length to report, default=30")
    pp.add_argument("--min-clip", type=int, default=10, help="Minimum soft/hard clip length to report, default=10")
    pp.add_argument("--bp-window", type=int, default=10, help="Breakpoint clustering window (bp), default=10")
    pp.add_argument("--splice-bed", help="Splice-site BED to mask soft-clip clusters")
    pp.add_argument("--splice-window", type=int, default=10, help="Window (bp) around splice sites to mask, default=10")
    pp.add_argument("--gene-bed", help="Gene BED3/4+ (name in col4) to annotate")
    pp.add_argument("--interchr-only-splits", action="store_true", default=True,
                    help="Emit only inter-chrom splits (recommended)")
    pp.add_argument("--include-intrachr-splits", dest="interchr_only_splits",
                    action="store_false", help="Also emit same-chrom splits")
    pp.add_argument("--pair-z", type=float, default=4.0, help="Z threshold for TLEN outliers")
    pp.add_argument("--tlen-thresh", type=int, help="Override TLEN threshold; skip sampling")
    pp.add_argument("--vcf", help="Write VCF 4.3 of indel clusters")
    pp.add_argument("--max-readnames", type=int, default=3, help="Max example read names per cluster, default=3")
    pp.add_argument("--progress-reads", type=int, default=1_000_000, help="Print progress every N reads")
    # SA throttles
    pp.add_argument("--sa-min-mapq", type=int, default=20)
    pp.add_argument("--sa-max-per-read", type=int, default=2)
    pp.add_argument("--sa-allow-contigs",
                    default=r"^(chr([0-9]{1,2}|1[0-9]|2[0-2]|X|Y|M)|[0-9]{1,2}|X|Y|MT|chrM)$")
    # Pair clustering knobs
    pp.add_argument("--pair-refine-min", type=int, default=20, help="Grid cell support threshold for local refine")
    pp.add_argument("--pair-include-unmapped-in-main", action="store_true", help="Also include one-end-unmapped pairs in main pairs clustering")
    pp.add_argument("--debug-pairs", action="store_true", help="Print pair bucket sizes (debug)")

    return p


def main():
    parser = build_parser()
    args = parser.parse_args()
    if args.mode == "lr":
        run_lr(args)
    elif args.mode == "sr":
        run_pe(args)
    else:
        parser.error("Unknown mode")


if __name__ == "__main__":
    main()


#!/usr/bin/env python3
"""
sr_isoSV_parser.py - Short-read parser that produces the same three TSV outputs
as the long-read script: <prefix>.indels.tsv, <prefix>.softclips.tsv, <prefix>.splits.tsv

- Coordinates in TSV: 1-based POS; DEL end reported (POS + length - 1)
- Gene annotation: expects a BED-like file; use --gene-name-col to pick the column with gene name.
"""
from __future__ import annotations
import argparse
import csv
import statistics
import pysam
import sys
import os
from collections import defaultdict
from dataclasses import dataclass
from bisect import bisect_right

# pysam CIGAR codes
M, I, D, N, S, H, P, EQ, X = 0, 1, 2, 3, 4, 5, 6, 7, 8


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
    seq: str | None = None


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


def ref_span_from_cigar(cigar_tuples):
    span = 0
    for op, length in cigar_tuples or []:
        if op in (M, EQ, X, D, N):
            span += length
    return span


def cigarstr_to_tuples(cigar_str: str):
    out = []
    num = ""
    opmap = {'M': M, 'I': I, 'D': D, 'N': N, 'S': S, 'H': H, 'P': P, '=': EQ, 'X': X}
    for ch in cigar_str:
        if ch.isdigit():
            num += ch
            continue
        if ch not in opmap:
            raise ValueError(f"Unknown CIGAR op: {ch}")
        out.append((opmap[ch], int(num)))
        num = ""
    if num:
        raise ValueError("Trailing digits in CIGAR")
    return out


def parse_sa_tag(sa: str):
    segs = []
    if not sa:
        return segs
    entries = [e for e in sa.strip().split(';') if e]
    for e in entries:
        try:
            rname, pos, strand, cig, mapq, nm = e.split(',')
            segs.append({
                "rname": rname,
                "pos": int(pos),
                "strand": strand,
                "cigar": cig,
                "mapq": int(mapq),
                "nm": int(nm)
            })
        except Exception:
            continue
    return segs


def walk_cigar_indels_and_clips(aln, min_ins, min_del, min_clip):
    chrom = aln.reference_name
    mapq = aln.mapping_quality or 0
    rname = aln.query_name
    try:
        read_seq = aln.query_sequence
    except Exception:
        read_seq = None

    cigar = aln.cigartuples or []
    ref_pos = aln.reference_start  # 0-based
    read_pos = 0

    # left clip
    if cigar:
        op0, len0 = cigar[0]
        if op0 in (S, H) and len0 >= min_clip:
            clipseq = (read_seq[:len0] if (op0 == S and read_seq) else "")
            if not is_poly_at(clipseq, min_len=min(20, len0)):
                yield ClipObs(chrom, aln.reference_start + 1, 'L', len0, mapq, rname)

    for op, length in cigar:
        if op in (M, EQ, X):
            ref_pos += length
            read_pos += length
        elif op == I:
            if length >= min_ins:
                ins_seq = (read_seq[read_pos: read_pos + length] if read_seq is not None else None)
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
        elif op in (H, P):
            pass

    # right clip
    if cigar:
        opn, lenn = cigar[-1]
        if opn in (S, H) and lenn >= min_clip:
            clipseq = (read_seq[-lenn:] if (opn == S and read_seq) else "")
            if not is_poly_at(clipseq, min_len=min(20, lenn)):
                yield ClipObs(chrom, aln.reference_end, 'R', lenn, mapq, rname)


### interval index for gene BED annotation ###
class IntervalIndex:
    def __init__(self):
        self.data = {}
    def add(self, chrom, start, end, name=None):
        self.data.setdefault(chrom, []).append((start, end, name))
    def finalize(self):
        for chrom in self.data:
            self.data[chrom].sort(key=lambda x: (x[0], x[1]))
    def overlaps(self, chrom, pos):
        arr = self.data.get(chrom, [])
        if not arr:
            return False, []
        starts = [s for s, _, _ in arr]
        i = bisect_right(starts, pos)
        hits = []
        for k in (i-1, i, i+1):
            if 0 <= k < len(arr):
                s,e,nm = arr[k]
                if s <= pos <= e:
                    hits.append((s,e,nm))
        return (len(hits) > 0), hits


class GeneAnnot:
    def __init__(self, name_col=4):
        # name_col: 1-based column index in BED that contains the gene name (user-specified)
        self.name_col = int(name_col)
        self.idx = IntervalIndex()

    def _store_for_chrom_variants(self, chrom, s1, e1, name):
        """
        Store the interval under multiple chromosome name variants to maximize the chance of match:
          - original chrom
          - if chrom startswith 'chr' also store without 'chr'
          - if chrom does not start with 'chr' also store with 'chr' prefix
        This avoids mismatches between BED using `chrN` and BAM/SA sometimes using `N`.
        """
        to_add = {chrom}
        if chrom.startswith("chr"):
            to_add.add(chrom[3:])
        else:
            to_add.add("chr" + chrom)
        for c in to_add:
            self.idx.add(c, s1, e1, name)

    def load_bed(self, path):
        """
        Load a BED-like file. Uses self.name_col (1-based) if present and not just '+'/'-'.
        Heuristic fallback: look for a non-strand token in columns 4..10 (prefer non-ENSG
        token as gene symbol; fall back to ENSG if needed). If nothing found, store empty.
        """
        count = 0
        with open(path) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 3:
                    continue
                chrom = p[0]
                try:
                    s1 = int(p[1]) + 1  # 0-based -> 1-based
                except Exception:
                    continue
                try:
                    e1 = int(p[2])
                except Exception:
                    continue

                # Try user-requested column first (if available)
                name = ""
                if len(p) >= self.name_col:
                    cand = p[self.name_col - 1].strip()
                    if cand and cand not in ("+", "-"):
                        name = cand

                # If that didn't look like a gene name, search heuristically in a small window
                if not name:
                    ens_candidate = ""
                    for idx in range(3, min(len(p), 10)):  # examine cols 4..10 (0-based indexes 3..9)
                        tok = p[idx].strip()
                        if not tok:
                            continue
                        if tok in ("+", "-"):
                            continue
                        # prefer short-ish tokens that look like gene symbols, not just ENS IDs
                        if not tok.upper().startswith("ENSG") and len(tok) <= 50:
                            name = tok
                            break
                        # keep ENSG as fallback if nothing better found
                        if tok.upper().startswith("ENSG") and not ens_candidate:
                            ens_candidate = tok
                    if not name and ens_candidate:
                        name = ens_candidate

                # final fallback: empty string (user requested empty if not found)
                if not name:
                    name = ""

                # store interval (and alternate chr variants)
                self._store_for_chrom_variants(chrom, s1, e1, name)
                count += 1
        self.idx.finalize()
        # small informational message
        print(f"[INFO] Loaded {count} intervals from gene BED (indexed chrom variants)", file=sys.stderr)

    def genes_at(self, chrom, pos):
        # try the chrom as-is first
        hit, arr = self.idx.overlaps(chrom, pos)
        if hit:
            return sorted({nm for _, _, nm in arr if nm})
        # fallback: try alt 'chr' / no-'chr' form
        alt = (chrom[3:] if chrom.startswith("chr") else ("chr" + chrom))
        hit2, arr2 = self.idx.overlaps(alt, pos)
        if hit2:
            return sorted({nm for _, _, nm in arr2 if nm})
        return []


### clustering helpers ###
def within(a, b, w):
    return abs(a - b) <= w

def length_similar(a, b, abs_slop=5, frac=0.10):
    return abs(a - b) <= max(abs_slop, int(frac * max(a, b)))


def cluster_indels(obs_list, bpw, max_readnames=3):
    clusters = []
    by_chr_type = defaultdict(list)
    for o in obs_list:
        by_chr_type[(o.chrom, o.svtype)].append(o)
    for (chrom, svtype), lst in by_chr_type.items():
        lst.sort(key=lambda x: (x.pos, x.length))
        used = [False]*len(lst)
        for i, oi in enumerate(lst):
            if used[i]:
                continue
            cluster = [i]
            for j in range(i+1, len(lst)):
                if used[j]:
                    continue
                oj = lst[j]
                if within(oi.pos, oj.pos, bpw) and length_similar(oi.length, oj.length):
                    cluster.append(j); used[j] = True
            pos_med = int(statistics.median([lst[k].pos for k in cluster]))
            len_med = int(statistics.median([lst[k].length for k in cluster]))
            mapqs = [lst[k].mapq for k in cluster]
            rnames = [lst[k].rname for k in cluster][:max_readnames]
            seqs = [lst[k].seq for k in cluster if (svtype == "INS" and lst[k].seq)]
            clusters.append({
                "chrom": chrom, "svtype": svtype, "pos": pos_med, "length": len_med,
                "support": len(cluster), "median_mapq": int(statistics.median(mapqs)) if mapqs else 0,
                "example_reads": ",".join(rnames),
                "example_ins_seq": (seqs[0] if seqs else "")
            })
    return clusters


def cluster_clips(obs_list, bpw, max_readnames=3):
    clusters = []
    by_chr_side = defaultdict(list)
    for o in obs_list:
        by_chr_side[(o.chrom, o.side)].append(o)
    for (chrom, side), lst in by_chr_side.items():
        lst.sort(key=lambda x: x.pos)
        used = [False]*len(lst)
        for i, oi in enumerate(lst):
            if used[i]:
                continue
            cluster = [i]
            for j in range(i+1, len(lst)):
                if used[j]:
                    continue
                oj = lst[j]
                if within(oi.pos, oj.pos, bpw):
                    cluster.append(j); used[j] = True
            pos_med = int(statistics.median([lst[k].pos for k in cluster]))
            mapqs = [lst[k].mapq for k in cluster]
            clips = [lst[k].cliplen for k in cluster]
            rnames = [lst[k].rname for k in cluster][:max_readnames]
            clusters.append({
                "chrom": chrom, "pos": pos_med, "side": side,
                "support": len(cluster),
                "median_mapq": int(statistics.median(mapqs)) if mapqs else 0,
                "median_cliplen": int(statistics.median(clips)) if clips else 0,
                "example_reads": ",".join(rnames)
            })
    return clusters


def cluster_splits(obs_list, bpw, max_readnames=3):
    clusters = []
    by_pair = defaultdict(list)
    for o in obs_list:
        if o.same_chrom:
            continue
        by_pair[(o.chrom1, o.chrom2)].append(o)
    for (c1, c2), lst in by_pair.items():
        lst.sort(key=lambda x: (x.pos1, x.pos2))
        used = [False]*len(lst)
        for i, oi in enumerate(lst):
            if used[i]:
                continue
            cluster = [i]
            for j in range(i+1, len(lst)):
                if used[j]:
                    continue
                oj = lst[j]
                if within(oi.pos1, oj.pos1, bpw) and within(oi.pos2, oj.pos2, bpw):
                    cluster.append(j); used[j] = True
            pos1_med = int(statistics.median([lst[k].pos1 for k in cluster]))
            pos2_med = int(statistics.median([lst[k].pos2 for k in cluster]))
            mapqs = [lst[k].mapq for k in cluster]
            rnames = [lst[k].rname for k in cluster][:max_readnames]
            notes = sorted(set(lst[k].note for k in cluster))
            clusters.append({
                "chrom1": c1, "pos1": pos1_med, "chrom2": c2, "pos2": pos2_med,
                "support": len(cluster), "median_mapq": int(statistics.median(mapqs)) if mapqs else 0,
                "notes": "|".join(notes),
                "example_reads": ",".join(rnames)
            })
    return clusters


def main():
    ap = argparse.ArgumentParser(description="Short-read parser that writes indels/softclips/splits TSVs")
    ap.add_argument("in_bam", help="Input sorted BAM (indexed recommended)")
    ap.add_argument("-o", "--out-prefix", required=True, help="Output prefix for TSVs (no extension)")
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--min-ins", type=int, default=20, help="Minimum insertion length to report")
    ap.add_argument("--min-del", type=int, default=20, help="Minimum deletion length to report")
    ap.add_argument("--min-clip", type=int, default=20, help="Minimum soft/hard clip length to report")
    ap.add_argument("--bp-window", type=int, default=10, help="Breakpoint clustering window (bp)")
    ap.add_argument("--sample-parse", type=int, default=None, help="Only parse this many reads")
    ap.add_argument("--gene-bed", help="Optional gene BED for simple annotation (name col default=4)")
    ap.add_argument("--gene-name-col", type=int, default=4, help="1-based BED column index to use as gene name (default 4)")
    ap.add_argument("--max-readnames", type=int, default=3, help="Max example read names per cluster")
    args = ap.parse_args()

    if not os.path.exists(args.in_bam):
        print(f"[ERROR] BAM not found: {args.in_bam}", file=sys.stderr)
        sys.exit(2)

    gene_annot = None
    if args.gene_bed:
        gene_annot = GeneAnnot(name_col=args.gene_name_col)
        gene_annot.load_bed(args.gene_bed)

    sam = pysam.AlignmentFile(args.in_bam, "rb")
    indel_obs = []
    clip_obs = []
    split_obs = []

    n_examined = 0
    for aln in sam.fetch(until_eof=True):
        n_examined += 1
        if args.sample_parse and n_examined > args.sample_parse:
            break
        if aln.is_unmapped or aln.is_secondary or aln.is_duplicate:
            continue
        if aln.mapping_quality is None or aln.mapping_quality < args.min_mapq:
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
                        mapq=aln.mapping_quality or 0, rname=aln.query_name,
                        same_chrom=(aln.reference_name == rname2),
                        note="primaryR->SAstart"
                    ))
                    split_obs.append(SplitObs(
                        chrom1=aln.reference_name, pos1=primary_left,
                        chrom2=rname2, pos2=(pos2 + span2 - 1 if span2 > 0 else pos2),
                        mapq=aln.mapping_quality or 0, rname=aln.query_name,
                        same_chrom=(aln.reference_name == rname2),
                        note="primaryL->SAend"
                    ))
            except Exception:
                pass

    sam.close()

    indel_clusters = cluster_indels(indel_obs, args.bp_window, max_readnames=args.max_readnames)
    clip_clusters  = cluster_clips(clip_obs, args.bp_window, max_readnames=args.max_readnames)
    split_clusters = cluster_splits(split_obs, args.bp_window, max_readnames=args.max_readnames)

    # annotate gene names if provided (exact overlap only)
    if gene_annot:
        for c in indel_clusters:
            genes = gene_annot.genes_at(c["chrom"], c["pos"])
            c["genes"] = ",".join(genes) if genes else ""
        for c in clip_clusters:
            genes = gene_annot.genes_at(c["chrom"], c["pos"])
            c["genes"] = ",".join(genes) if genes else ""
        for c in split_clusters:
            genes_l = gene_annot.genes_at(c["chrom1"], c["pos1"])
            genes_r = gene_annot.genes_at(c["chrom2"], c["pos2"])
            c["genes_left"] = ",".join(genes_l) if genes_l else ""
            c["genes_right"] = ",".join(genes_r) if genes_r else ""
    else:
        for c in indel_clusters:
            c["genes"] = ""
        for c in clip_clusters:
            c["genes"] = ""
        for c in split_clusters:
            c["genes_left"] = ""
            c["genes_right"] = ""

    out_pre = args.out_prefix
    with open(f"{out_pre}.indels.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["chrom","pos","end","svtype","length","support","median_mapq","example_reads","example_ins_seq","genes"])
        for c in sorted(indel_clusters, key=lambda x: (x["chrom"], x["pos"], x["svtype"])):
            end = c["pos"] + c["length"] - 1 if c["svtype"] == "DEL" else c["pos"]
            w.writerow([c["chrom"], c["pos"], end, c["svtype"], c["length"], c["support"], c["median_mapq"], c["example_reads"], c.get("example_ins_seq",""), c.get("genes","")])

    with open(f"{out_pre}.softclips.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["chrom","pos","side","support","median_mapq","median_cliplen","example_reads","genes"])
        for c in sorted(clip_clusters, key=lambda x: (x["chrom"], x["pos"], x["side"])):
            w.writerow([c["chrom"], c["pos"], c["side"], c["support"], c["median_mapq"], c["median_cliplen"], c["example_reads"], c.get("genes","")])

    with open(f"{out_pre}.splits.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["chrom1","pos1","chrom2","pos2","support","median_mapq","notes","example_reads","genes_left","genes_right"])
        for c in sorted(split_clusters, key=lambda x: (x["chrom1"], x["pos1"], x["chrom2"], x["pos2"])):
            w.writerow([c["chrom1"], c["pos1"], c["chrom2"], c["pos2"], c["support"], c["median_mapq"], c["notes"], c["example_reads"], c.get("genes_left",""), c.get("genes_right","")])

    print(f"[indels] clusters: {len(indel_clusters)}", file=sys.stderr)
    print(f"[softclips] clusters: {len(clip_clusters)}", file=sys.stderr)
    print(f"[splits inter-chr] clusters: {len(split_clusters)}", file=sys.stderr)


if __name__ == "__main__":
    main()

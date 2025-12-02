"""
Utility functions and data structures for the IsoSV integrated pipeline.

This module contains the core logic and helper functions used by the main
pipeline script (`run_pipeline.py`). It includes:
- Data classes for representing different types of SV observations (Indels,
  Clips, Splits).
- Functions for parsing pysam alignment records to extract SV candidates from
  CIGAR strings and SA tags.
- Clustering algorithms to group raw SV observations into meaningful events.
- Annotation logic to intersect SVs with gene models and produce functional
  classifications.
- VCF writing utilities to format the final output.
"""
import argparse
import csv
import statistics
from collections import defaultdict
from dataclasses import dataclass
import pysam
import sys
from bisect import bisect_left, bisect_right
import os
import pandas as pd

# pysam CIGAR operator codes
M, I, D, N, S, H, P, EQ, X, B = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

# define dataclasses for the three types of observations
@dataclass
class IndelObs:
    """Represents an InDel observation from a single read."""
    chrom: str
    pos: int
    svtype: str
    length: int
    mapq: int
    rname: str
    seq: str | None = None

@dataclass
class ClipObs:
    """Represents a soft/hard clip observation from a single read."""
    chrom: str
    pos: int
    side: str
    cliplen: int
    mapq: int
    rname: str

@dataclass
class SplitObs:
    """Represents a split-read observation from an SA tag."""
    chrom1: str
    pos1: int
    chrom2: str
    pos2: int
    mapq: int
    rname: str
    same_chrom: bool
    note: str = ""

# define a class for gene annotation, similar to the one in step_a_IsoParser/scripts/lr_isoSV_parser.py
class GeneAnnot:
    """Handles loading and querying a gene annotation BED file."""
    def __init__(self):
        self.idx = IntervalIndex()

    def load_bed(self, path):
        with open(path) as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                p = line.rstrip("").split("	")
                if len(p) < 3:
                    continue
                chrom, s1, e1 = p[0], int(p[1]) + 1, int(p[2])
                name = p[3] if len(p) >= 4 and p[3] else f"gene:{chrom}:{s1}-{e1}"
                self.idx.add(chrom, s1, e1, name)
        self.idx.finalize()

    def genes_at(self, chrom, pos):
        hit, arr = self.idx.overlaps(chrom, pos)
        return sorted({nm for _, _, nm in arr if nm}) if hit else []

def is_poly_at(seq: str, min_len: int = 20, frac: float = 0.8) -> bool:
    """
    Checks if a sequence is likely a poly-A or poly-T tail.

    Args:
        seq: The sequence to check.
        min_len: The minimum length for the check to apply.
        frac: The minimum fraction of A's or T's to be considered poly-A/T.

    Returns:
        True if the sequence is likely a poly-A/T tail, False otherwise.
    """
    if not seq or len(seq) < min_len:
        return False
    u = seq.upper()
    a = u.count("A"); t = u.count("T")
    return (a / len(u) >= frac) or (t / len(u) >= frac)

def ref_span_from_cigar(cigar_tuples) -> int:
    """Calculates the span of a CIGAR string on the reference genome."""
    span = 0
    for op, length in cigar_tuples or []:
        if op in (M, EQ, X, D, N):
            span += length
    return span

def walk_cigar_indels_and_clips(aln, min_ins, min_del, min_clip):
    """
    Scan a read's CIGAR to yield Indel and Clip observations.

    This is adapted from the long-read parser in `step_a_IsoParser` so that
    the integrated pipeline sees the same indel / soft-clip events.
    """
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
                # report the left breakpoint as 1-based start
                yield ClipObs(chrom, aln.reference_start + 1, 'L', len0, mapq, rname)

    for op, length in cigar:
        if op in (M, EQ, X):
            ref_pos += length
            read_pos += length
        elif op == I:
            # insertion relative to reference at ref_pos
            ins_seq = None
            if length >= min_ins and read_seq is not None:
                ins_seq = read_seq[read_pos: read_pos + length]
            if length >= min_ins:
                yield IndelObs(chrom, ref_pos + 1, "INS", length, mapq, rname, ins_seq)
            read_pos += length
        elif op == D:
            # deletion relative to reference starting at ref_pos
            if length >= min_del:
                yield IndelObs(chrom, ref_pos + 1, "DEL", length, mapq, rname, None)
            ref_pos += length
        elif op == N:
            # skipped region (intron)
            ref_pos += length
        elif op == S:
            # soft clip inside the read (already handled at ends)
            read_pos += length
        elif op in (H, P, B):
            # hard clip / padding / back
            pass

    # Right clip
    if cigar:
        opn, lenn = cigar[-1]
        if opn in (S, H) and lenn >= min_clip:
            clipseq = (read_seq[-lenn:] if (opn == S and read_seq is not None) else "")
            if not is_poly_at(clipseq):
                # report the right breakpoint as reference_end (1-based)
                yield ClipObs(chrom, aln.reference_end, 'R', lenn, mapq, rname)

def parse_sa_tag(sa: str):
    """
    Parse the SA tag into a list of dicts.

    Matches the format used in `lr_isoSV_parser.py` / `sr_isoSV_parser.py` so
    that downstream split-read handling can treat entries uniformly.
    """
    segs = []
    if not sa:
        return segs
    entries = [e for e in sa.strip().split(";") if e]
    for e in entries:
        try:
            rname, pos, strand, cig, mapq, nm = e.split(",")
        except ValueError:
            # malformed SA entry – skip
            continue
        segs.append({
            "rname": rname,
            "pos": int(pos),
            "strand": strand,
            "cigar": cig,
            "mapq": int(mapq),
            "nm": int(nm),
        })
    return segs


def cigarstr_to_tuples(cigar_str: str):
    """
    Convert a CIGAR string from an SA tag into (op, length) tuples
    using the same numeric op codes as pysam.
    """
    out = []
    num = ""
    opmap = {'M': M, 'I': I, 'D': D, 'N': N, 'S': S, 'H': H, 'P': P, '=': EQ, 'X': X}
    for ch in cigar_str:
        if ch.isdigit():
            num += ch
            continue
        if ch not in opmap:
            raise ValueError(f"Unknown CIGAR op in SA: {ch}")
        if not num:
            raise ValueError("Missing length before CIGAR op in SA")
        out.append((opmap[ch], int(num)))
        num = ""
    if num:
        raise ValueError("Trailing digits in CIGAR string")
    return out

def within(a, b, w):
    return abs(a - b) <= w

def length_similar(a, b, abs_slop=5, frac=0.10):
    return abs(a - b) <= max(abs_slop, int(frac * max(a, b)))

def cluster_indels(obs_list, bpw, gene_annot=None, max_readnames=3):
    """
    Clusters Indel observations based on position and length.

    Returns:
        A list of dictionaries, where each dictionary represents a cluster.
    """
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
            rnames = [lst[k].rname for k in cluster][:max_readnames]
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

def cluster_clips(obs_list, bpw, gene_annot=None, max_readnames=3):
    """
    Clusters soft/hard clip observations based on position.

    Returns:
        A list of dictionaries, where each dictionary represents a cluster.
    """
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
            rnames = [lst[k].rname for k in cluster][:max_readnames]
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
def cluster_splits(obs_list, bpw, gene_annot=None, max_readnames=3):
    """
    Clusters split-read observations for inter-chromosomal events.

    Returns:
        A list of dictionaries, where each dictionary represents a cluster.
    """
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
            rnames = [lst[k].rname for k in cluster][:max_readnames]
            notes = sorted(set(lst[k].note for k in cluster))
            genes1, genes2 = [], []
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

def load_tx_tree(filepath):
    """
    Load a pickled IntervalTree from disk.

    Args:
        filepath: Path to the .pkl file containing the IntervalTree.

    Returns:
        The loaded IntervalTree object.
    """
    import pickle
    if not os.path.exists(filepath):
        raise FileNotFoundError("Transcript tree not found! ")
    with open(filepath, 'rb') as f:
        return pickle.load(f)

def annotate_candidates(candidates: pd.DataFrame, tx_tree):
    """
    Annotates SV candidates based on their overlap with gene features.

    Args:
        candidates: A DataFrame of clustered SVs.
        tx_tree: A loaded IntervalTree containing gene and transcript information.

    Returns:
        A list of tuples, with each tuple containing the annotated fields for a VCF.
    """
    results = []
    for _, row in candidates.iterrows():
        chrom = row['chrom']
        start = int(row['start'])
        end = int(row['end'])
        svtype = row['type']
        support = int(row['support'])
        median_svlen = int(row['median_sv_len'])
        
        if not str(chrom).startswith("chr"):
            chrom = "chr" + str(chrom)
        
        if chrom not in tx_tree:
            continue
        
        query_start, query_end = (start, end) if svtype == 'DEL' else (start, start + 1)
        overlapping_hits = tx_tree[chrom][query_start:query_end]
        
        if not overlapping_hits:
            results.append((chrom, start, end, svtype, "Intergenic", support, median_svlen, "Intergenic", "na", "na", "na"))
            continue

        touched_genes, hit_exons, hit_introns = set(), [], []
        all_tx_aliases, all_gene_biotypes, all_tx_biotypes = set(), set(), set()

        for hit in overlapping_hits:
            region_str, tx_alias, gene_biotype, tx_biotype = hit.data
            try:
                gene_name, _, region_type = region_str.split(":", 2)
            except ValueError:
                continue

            touched_genes.add(gene_name)
            all_tx_aliases.add(tx_alias)
            all_gene_biotypes.add(gene_biotype)
            all_tx_biotypes.add(tx_biotype)

            if 'exon' in region_type:
                hit_exons.append(hit)
            elif 'intron' in region_type:
                hit_introns.append(hit)
        
        final_annotation = "Unclassified"
        if svtype == 'DEL':
            if len(touched_genes) > 1:
                final_annotation = "Gene_Fusion"
            elif len(touched_genes) == 1:
                is_exonic_deletion = any(start <= exon.begin and end >= exon.end for exon in hit_exons)
                if is_exonic_deletion:
                    final_annotation = "Exonic_Deletion"
                else:
                    is_canonical_splicing = any(abs(start - intron.begin) <= 5 and abs(end - intron.end) <= 5 for intron in hit_introns)
                    final_annotation = "Canonical_Splicing" if is_canonical_splicing else "Intronic/Novel_Splicing"
        elif svtype == 'INS':
            final_annotation = "Insertion_in_Fusion_Region" if len(touched_genes) > 1 else "Intragenic_Insertion"

        region_field = ",".join(sorted(list(touched_genes))) or "na"
        tx_alias_field = ",".join(sorted(list(all_tx_aliases))) or "na"
        gene_biotype_field = ",".join(sorted(list(all_gene_biotypes))) or "na"
        tx_biotype_field = ",".join(sorted(list(all_tx_biotypes))) or "na"

        results.append((
            chrom, start, end, svtype, final_annotation, support, median_svlen,
            region_field, tx_alias_field, gene_biotype_field, tx_biotype_field
        ))
            
    return results

def write_to_vcf(annotated_df, outpath):
    """
    Writes annotated SVs to a VCF file.

    Args:
        annotated_df: A DataFrame containing the annotated SVs.
        outpath: The path for the output VCF file.
    """
    with open(outpath, "w") as vcf_out:
        header_lines = [
            "##fileformat=VCFv4.2", "##source=IsoSV", "##FILTER=<ID=PASS,Description=\"All filters passed\">",
            "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">",
            "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">",
            "##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description=\"Support count for this SV\">",
            "##INFO=<ID=ANN_TYPE,Number=1,Type=String,Description=\"Functional annotation of the SV\">",
            "##INFO=<ID=REGION,Number=1,Type=String,Description=\"Gene region\">",
            "##INFO=<ID=TX_NAME,Number=1,Type=String,Description=\"Transcript name\">",
            "##INFO=<ID=GENE_BIOTYPE,Number=1,Type=String,Description=\"Gene biotype\">",
            "##INFO=<ID=TX_BIOTYPE,Number=1,Type=String,Description=\"Transcript biotype\">",
            "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO"
        ]
        vcf_out.write("".join(header_lines) + "\n")

        for _, row in annotated_df.iterrows():
            chrom, start, stop, svtype, ann_type, support, svlen, region, tx_alias, biotype_gene, biotype_tx = row
            record = f"{chrom}	{start}	.	N	<{svtype.upper()}>	.	PASS	END={stop};SVTYPE={svtype};SVLEN={svlen};SUPPORT={support};ANN_TYPE={ann_type};REGION={region}"
            if tx_alias != "na": record += f";TX_NAME={tx_alias}"
            if biotype_gene != "na": record += f";GENE_BIOTYPE={biotype_gene}"
            if biotype_tx != "na": record += f";TX_BIOTYPE={biotype_tx}"
            vcf_out.write(record + "\n")

class IntervalIndex:
    def __init__(self):
        self.data = {}  # chrom -> list of (start, end, name)

    def add(self, chrom, start, end, name=None):
        self.data.setdefault(chrom, []).append((start, end, name))

    def finalize(self):
        for chrom in self.data:
            self.data[chrom].sort(key=lambda x: (x[0], x[1]))
    def overlaps(self, chrom, pos):
        hits = []
        for start, end, name in self.data.get(chrom, []):
            if start <= pos <= end:
                hits.append((start, end, name))
        return (len(hits) > 0), hits

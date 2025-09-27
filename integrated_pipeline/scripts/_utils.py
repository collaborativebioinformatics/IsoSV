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

@dataclass
class IndelObs:
    chrom: str
    pos: int
    svtype: str
    length: int
    mapq: int
    rname: str
    seq: str | None = None

@dataclass
class ClipObs:
    chrom: str
    pos: int
    side: str
    cliplen: int
    mapq: int
    rname: str

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

class GeneAnnot:
    def __init__(self):
        self.idx = IntervalIndex()

    def load_bed(self, path):
        with open(path) as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                p = line.rstrip("
").split("	")
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
    if not seq or len(seq) < min_len:
        return False
    u = seq.upper()
    a = u.count("A"); t = u.count("T")
    return (a / len(u) >= frac) or (t / len(u) >= frac)

def ref_span_from_cigar(cigar_tuples) -> int:
    span = 0
    for op, length in cigar_tuples or []:
        if op in (M, EQ, X, D, N):
            span += length
    return span

def walk_cigar_indels_and_clips(aln, min_ins, min_del, min_clip):
    if not is_poly_at(clipseq):
        yield ClipObs(chrom, aln.reference_end, 'R', lenn, mapq, rname)

def parse_sa_tag(sa: str):
    segs = []
    for seg in sa.split(","):
        if ";" in seg:
            seg, qual = seg.split(";")
            segs.append((seg, qual))
        else:
            segs.append((seg, ""))
    return segs

def cigarstr_to_tuples(cigar_str: str):
    out = []
    num = ""
    for c in cigar_str:
        if c.isdigit():
            num += c
        else:
            if num:
                out.append((int(num), c))
                num = ""
            else:
                out.append((1, c))
    if num:
        raise ValueError("Trailing digits in CIGAR string")
    return out

def within(a, b, w):
    return abs(a - b) <= w

def length_similar(a, b, abs_slop=5, frac=0.10):
    return abs(a - b) <= max(abs_slop, int(frac * max(a, b)))

def cluster_indels(obs_list, bpw, gene_annot=None, max_readnames=3):
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
    """Load the IntervalTree from disk"""
    import pickle
    if not os.path.exists(filepath):
        raise FileNotFoundError("Transcript tree not found! ")
    with open(filepath, 'rb') as f:
        return pickle.load(f)

def annotate_candidates(candidates: pd.DataFrame, tx_tree):
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
        vcf_out.write("
".join(header_lines) + "
")

        for _, row in annotated_df.iterrows():
            chrom, start, stop, svtype, ann_type, support, svlen, region, tx_alias, biotype_gene, biotype_tx = row
            record = f"{chrom}	{start}	.	N	<{svtype.upper()}>	.	PASS	END={stop};SVTYPE={svtype};SVLEN={svlen};SUPPORT={support};ANN_TYPE={ann_type};REGION={region}"
            if tx_alias != "na": record += f";TX_NAME={tx_alias}"
            if biotype_gene != "na": record += f";GENE_BIOTYPE={biotype_gene}"
            if biotype_tx != "na": record += f";TX_BIOTYPE={biotype_tx}"
            vcf_out.write(record + "
")

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

"""
Lightweight clustering utilities for structural-variant evidence.

Provides three functions that group read-level observations into concise
summary clusters: `cluster_indels` (clusters indel observations by
chromosome, position and length), `cluster_clips` (clusters soft/hard clips
by chromosome, side and position), and `cluster_splits` (clusters inter-
chromosomal split-read pairs by paired positions). Each function returns a
list of dictionaries containing aggregated fields such as median position/
length, support count, median MAPQ, example read names, and optional gene
annotations.
"""

import statistics
from collections import defaultdict
from utils import within, length_similar

def cluster_indels(obs_list, bpw, trees=None, max_readnames=3):
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

            if svtype == "DEL":
                genes = set(x.data["gene_name"] for x in trees[chrom][oi.pos:oi.pos+oi.length])
            else:
                genes = set(x.data["gene_name"] for x in trees[chrom][oi.pos:oi.pos+1])
                
            clusters.append({
                "chrom": chrom, "svtype": svtype, "pos": pos_med, "length": len_med,
                "support": len(cluster), "median_mapq": int(statistics.median(mapqs)),
                "example_reads": ",".join(rnames),
                "example_ins_seq": (seqs[0] if seqs else ""),
                "genes": ",".join(genes) if genes else ""
            })
    return clusters

def cluster_clips(obs_list, bpw, trees=None, max_readnames=3):
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
            clips_med = int(statistics.median(clips))
            if trees:
                genes = set(x.data["gene_name"] for x in trees[chrom][pos_med-clips_med:pos_med+clips_med])
            clusters.append({
                "chrom": chrom, "pos": pos_med, "side": side,
                "support": len(cluster),
                "median_mapq": int(statistics.median(mapqs)),
                "median_cliplen": clips_med,
                "example_reads": ",".join(rnames),
                "genes": ",".join(genes) if genes else ""
            })
    return clusters


def cluster_splits(obs_list, bpw, trees=None, max_readnames=3):
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
            if trees:
                genes1 = [iv.data["gene_name"] for iv in trees[c1][pos1_med-1:pos1_med+1]]
                genes2 = [iv.data["gene_name"] for iv in trees[c2][pos2_med-1:pos2_med+1]]
            clusters.append({
                "chrom1": c1, "pos1": pos1_med, "chrom2": c2, "pos2": pos2_med,
                "support": len(cluster), "median_mapq": int(statistics.median(mapqs)),
                "notes": "|".join(notes),
                "example_reads": ",".join(rnames),
                "genes_left": ",".join(genes1) if genes1 else "",
                "genes_right": ",".join(genes2) if genes2 else ""
            })
    return clusters



#!/usr/bin/env python3
"""
iso_parser.py - parser + clustering + optional GTF annotation

This updated version adds:
- an `--out-dir` option so all outputs (candidates, clusters, logs) go into a single folder
- a short automatic summary printed to stdout and saved as `<out_prefix>.summary.txt`
  which contains per-type candidate counts and simple percentages
- richer docstrings and in-line comments so team members can understand each function

Features (unchanged):
- Robust BAM classification + parsing (skips & logs corrupted records)
- Per-read candidate extraction (INS/DEL/SOFTCLIP/LARGE_N/SPLIT)
- Clustering of per-read candidates into region-level events
- Optional GTF annotation of clusters (Gencode v43 compatible)

Usage highlights (example):
    python iso_parser.py input.bam results.tsv --out-dir results/ --insert 1 --delete 1 \
        --softclip 1 --sample-parse 20000 --cluster-window 50 --min-support 2 --salvage

Outputs (in the out-dir or next to `results.tsv` if --out-dir omitted):
  results.tsv                    - per-read candidate TSV (tab-separated)
  results.tsv.errors.log         - parse/salvage errors
  results.clusters.tsv           - cluster-level TSV
  results.clusters.txt           - cluster-level TXT (same content, tool-friendly ext.)
  results.summary.txt            - short summary (counts by type + basic stats)

"""
import pysam
import argparse
import statistics
import sys
import os
import gzip
import re
import json
from collections import defaultdict, Counter
from bisect import bisect_left

# CIGAR op codes mapping used in logic
# 0=M,1=I,2=D,3=N,4=S,5=H,6=P,7==,8=X
OP_CONSUME_REF = {0, 2, 3, 7, 8}


def infer_bam_type_from_header(bam_path):
    """
    Inspect the BAM header's @RG platform (PL) tags to guess short vs long read BAM.

    Returns:
        str or None: "short" if header indicates Illumina; "long" for PacBio/ONT; None if unknown.
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    rgs = bam.header.get('RG', [])
    platforms = set()
    for rg in (rgs or []):
        pl = rg.get('PL')
        if pl:
            platforms.add(pl.upper())
    bam.close()
    if not platforms:
        return None
    if any("ILLUMINA" in p for p in platforms):
        return "short"
    if any(p in ("PACBIO", "PACBIO_RS", "PACBIO_SMRT", "PACBIO_HIFI", "ONT", "OXFORD_NANOPORE") for p in platforms):
        return "long"
    return None


def infer_bam_type_by_readlen(bam_path, sample_reads=5000, readlen_cutoff=1000):
    """
    Robustly sample reads from the BAM and compute median/mean read lengths.
    Skips reads that raise exceptions so partially corrupted BAMs still yield a classification.

    Args:
        bam_path (str): path to BAM file
        sample_reads (int): how many reads to sample (stops early if EOF)
        readlen_cutoff (int): median read length threshold to call a BAM "long"

    Returns:
        (bam_type_or_None, median_len, mean_len, n_sampled, n_skipped)
    """
    lengths = []
    skipped = 0
    n_seen = 0
    bam = pysam.AlignmentFile(bam_path, "rb")
    it = bam.fetch(until_eof=True)
    while len(lengths) < sample_reads:
        try:
            r = next(it)
        except StopIteration:
            break
        except Exception:
            skipped += 1
            n_seen += 1
            continue
        n_seen += 1
        qlen = r.query_length
        if qlen is None:
            seq = r.query_sequence
            if seq is not None:
                qlen = len(seq)
        if qlen is not None and qlen > 0:
            lengths.append(qlen)
    bam.close()
    if not lengths:
        return None, 0, 0.0, len(lengths), skipped
    median_len = int(statistics.median(lengths))
    mean_len = float(statistics.mean(lengths))
    bam_type = "short" if median_len < readlen_cutoff else "long"
    return bam_type, median_len, mean_len, len(lengths), skipped


def choose_thresholds(bam_type, short_defaults=(30, 30, 30), long_defaults=(200, 200, 200)):
    """
    Choose default size thresholds (INS, DEL, SOFTCLIP) based on inferred BAM type.
    Override with command-line args if provided.

    Returns:
        tuple: (ins_default, del_default, sc_default)
    """
    return long_defaults if bam_type == "long" else short_defaults


def _get_read_length(read):
    """Return the read length using pysam's query_length or fallback to sequence length."""
    qlen = read.query_length
    if qlen is None:
        seq = read.query_sequence
        if seq is not None:
            qlen = len(seq)
    if qlen is None:
        qlen = 0
    return int(qlen)


def _write_candidate_line(fout, chrom, start, end, qname, sv_type, sv_len, note, read_len, cigar):
    """Helper to write a tab-delimited candidate line to the open file handle."""
    fout.write(f"{chrom}\t{start}\t{end}\t{qname}\t{sv_type}\t{sv_len}\t{note}\t{read_len}\t{cigar}\n")


# ---------------- GTF helpers ----------------

def parse_gtf_exons(gtf_path):
    """
    Parse a GTF/GFF file (optionally gzipped) and return a dict mapping chrom -> list of (start, end, gene_name).
    We index by start position (sorted list) to allow fast interval scanning later.
    """
    exons = defaultdict(list)
    openf = gzip.open if gtf_path.endswith(".gz") else open
    # attribute regex: key "value" pairs
    attr_re = re.compile(r'(\S+)\s+"([^"]+)"')
    with openf(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = parts[:9]
            if feature.lower() not in ("exon", "cds", "transcript", "gene"):
                continue
            gene_name = None
            for m in re.finditer(attr_re, attrs):
                key, val = m.group(1), m.group(2)
                if key == "gene_name":
                    gene_name = val
                    break
            if gene_name is None:
                m = re.search(r'gene_id\s+"([^"]+)"', attrs)
                if m:
                    gene_name = m.group(1)
            try:
                s = int(start) - 1
                e = int(end)
            except Exception:
                continue
            exons[chrom].append((s, e, gene_name or "NA"))
    for chrom in exons:
        exons[chrom].sort(key=lambda x: x[0])
    return exons


def gtf_find_genes_for_interval(exons_dict, chrom, qstart, qend):
    """
    Given a pre-parsed exons dict, return a set of gene names overlapping [qstart, qend).
    We use a small window around the bisect lookup to limit scanning.
    """
    genes = set()
    if chrom not in exons_dict:
        return genes
    arr = exons_dict[chrom]
    starts = [iv[0] for iv in arr]
    i = bisect_left(starts, qend)
    lo = max(0, i - 1000)
    hi = min(len(arr), i + 100)
    for j in range(lo, hi):
        s, e, g = arr[j]
        if e <= qstart:
            continue
        if s >= qend:
            continue
        genes.add(g)
    return genes


# ---------------- clustering helpers ----------------

def cluster_candidates(candidates, cluster_window=50, min_support=1):
    """
    Simple one-pass clustering grouped by (chrom, sv_type):
      - sort candidates by start
      - merge candidates if the next candidate start <= cur_end + cluster_window
      - record distinct read support and median sv_len

    Returns: list of cluster dicts with keys: chrom, sv_type, start, end, support, reads, median_sv_len
    """
    clusters = []
    by_ch_type = defaultdict(list)
    for c in candidates:
        key = (c['chrom'], c['sv_type'])
        by_ch_type[key].append(c)

    for (chrom, sv_type), items in by_ch_type.items():
        items.sort(key=lambda x: x['start'])
        cur_start = None
        cur_end = None
        cur_reads = set()
        cur_sv_lens = []
        for it in items:
            s = it['start']
            e = it['end'] if it['end'] is not None else it['start']
            if cur_start is None:
                cur_start = s
                cur_end = e
                cur_reads = {it['read_name']}
                cur_sv_lens = [it['sv_len']]
            else:
                if s <= cur_end + cluster_window:
                    cur_end = max(cur_end, e)
                    cur_reads.add(it['read_name'])
                    cur_sv_lens.append(it['sv_len'])
                else:
                    support = len(cur_reads)
                    if support >= min_support:
                        clusters.append({
                            'chrom': chrom,
                            'sv_type': sv_type,
                            'start': cur_start,
                            'end': cur_end,
                            'support': support,
                            'reads': sorted(cur_reads),
                            'median_sv_len': int(statistics.median(cur_sv_lens)) if cur_sv_lens else 0
                        })
                    cur_start = s
                    cur_end = e
                    cur_reads = {it['read_name']}
                    cur_sv_lens = [it['sv_len']]
        if cur_start is not None:
            support = len(cur_reads)
            if support >= min_support:
                clusters.append({
                    'chrom': chrom,
                    'sv_type': sv_type,
                    'start': cur_start,
                    'end': cur_end,
                    'support': support,
                    'reads': sorted(cur_reads),
                    'median_sv_len': int(statistics.median(cur_sv_lens)) if cur_sv_lens else 0
                })
    clusters.sort(key=lambda x: (x['chrom'], x['start'], x['sv_type']))
    return clusters


# ---------------- main parser (collect candidates then cluster) ----------------

def parse_bam_and_cluster(bam_path, out_tsv, final_type, median_len, mean_len,
                          mq_thresh=20, ins_thresh=30, del_thresh=30, sc_thresh=30,
                          large_N_thresh=None, sample_limit=None, salvage=False,
                          cluster_window=None, min_support=1, gtf_exons=None, out_dir=None):
    """
    Parse BAM and extract per-read SV-like candidates. Optionally cluster them and annotate.

    Behavior:
      - writes per-read TSV at `out_tsv` (or into out_dir if provided)
      - writes clusters files if clustering requested
      - logs parse errors to `<out_prefix>.errors.log`
      - writes a short summary file `<out_prefix>.summary.txt`

    Returns:
        (candidates_list, clusters_list)
    """
    # prepare output folder and consistent file prefixes
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        base = os.path.basename(out_tsv)
        out_prefix = os.path.join(out_dir, os.path.splitext(base)[0])
        out_tsv = os.path.join(out_dir, base)
    else:
        out_prefix = os.path.splitext(out_tsv)[0]

    bam = pysam.AlignmentFile(bam_path, "rb")
    fout = open(out_tsv, "w")
    fout.write(f"#BAM_CLASSIFICATION={final_type}\n")
    fout.write(f"#SAMPLED_MEDIAN_READLEN={median_len}\n")
    fout.write(f"#SAMPLED_MEAN_READLEN={mean_len:.1f}\n")
    fout.write(f"#THRESHOLDS_INS={ins_thresh}\tDEL={del_thresh}\tSOFTCLIP={sc_thresh}\tMQ={mq_thresh}\n")
    fout.write("chrom\tpos_start\tpos_end\tread_name\tsv_type\tsv_len\tnote\tread_len\tcigar\n")

    err_log = out_prefix + ".errors.log"
    elog = open(err_log, "a")

    total_examined = 0
    skipped = 0
    candidates = []
    counts = Counter()

    for read in bam.fetch(until_eof=True):
        total_examined += 1
        if sample_limit and total_examined > sample_limit:
            break

        try:
            if read.is_unmapped:
                continue
            if read.mapping_quality is None or read.mapping_quality < mq_thresh:
                continue
            if read.reference_id is None or read.reference_id < 0:
                continue

            chrom = bam.get_reference_name(read.reference_id)
            ref_pos = read.reference_start
            try:
                sa = read.get_tag("SA")
            except KeyError:
                sa = None
            except Exception as e_tag:
                sa = None
                elog.write(f"{getattr(read,'query_name','UNKNOWN')}\tSA_read_error\t{str(e_tag)}\n")

            cigarstr = read.cigarstring or "."
            read_len = _get_read_length(read)
            if read.cigartuples is None:
                continue

            ref_cursor = ref_pos
            for idx, (op, length) in enumerate(read.cigartuples):
                # insertion in CIGAR: I
                if op == 1 and length >= ins_thresh:
                    _write_candidate_line(fout, chrom, ref_cursor, ref_cursor, read.query_name, "INS", length, "CIGAR_I", read_len, cigarstr)
                    candidates.append({'chrom': chrom, 'start': ref_cursor, 'end': ref_cursor, 'read_name': read.query_name, 'sv_type': 'INS', 'sv_len': length})
                    counts['INS'] += 1
                # deletion on reference: D
                if op == 2 and length >= del_thresh:
                    _write_candidate_line(fout, chrom, ref_cursor, ref_cursor + length - 1, read.query_name, "DEL", length, "CIGAR_D", read_len, cigarstr)
                    candidates.append({'chrom': chrom, 'start': ref_cursor, 'end': ref_cursor + length - 1, 'read_name': read.query_name, 'sv_type': 'DEL', 'sv_len': length})
                    counts['DEL'] += 1
                # soft-clip: S
                if op == 4 and length >= sc_thresh:
                    note = "CIGAR_S_left" if idx == 0 else ("CIGAR_S_right" if idx == len(read.cigartuples) - 1 else "CIGAR_S_internal")
                    _write_candidate_line(fout, chrom, ref_cursor, ref_cursor, read.query_name, "SOFTCLIP", length, note, read_len, cigarstr)
                    candidates.append({'chrom': chrom, 'start': ref_cursor, 'end': ref_cursor, 'read_name': read.query_name, 'sv_type': 'SOFTCLIP', 'sv_len': length})
                    counts['SOFTCLIP'] += 1
                # N (skipped region / intron in RNA-seq); optional large_N_thresh
                if op == 3 and large_N_thresh and length >= large_N_thresh:
                    _write_candidate_line(fout, chrom, ref_cursor, ref_cursor + length - 1, read.query_name, "LARGE_N", length, "CIGAR_N", read_len, cigarstr)
                    candidates.append({'chrom': chrom, 'start': ref_cursor, 'end': ref_cursor + length - 1, 'read_name': read.query_name, 'sv_type': 'LARGE_N', 'sv_len': length})
                    counts['LARGE_N'] += 1
                if op in OP_CONSUME_REF:
                    ref_cursor += length

            # SA tag (supplementary/split alignment) - treat as split-support
            if sa is not None:
                _write_candidate_line(fout, chrom, read.reference_start, read.reference_end, read.query_name, "SPLIT", 0, "SA_tag", read_len, cigarstr)
                candidates.append({'chrom': chrom, 'start': read.reference_start, 'end': read.reference_end, 'read_name': read.query_name, 'sv_type': 'SPLIT', 'sv_len': 0})
                counts['SPLIT'] += 1

        except Exception as e:
            skipped += 1
            qn = getattr(read, "query_name", "UNKNOWN")
            elog.write(f"{qn}\tparse_error\t{str(e)}\n")
            if salvage:
                try:
                    # conservative salvage: inspect CIGAR tuples if present and write using a minimal chrom placeholder
                    if getattr(read, 'cigartuples', None):
                        chrom = None
                        if getattr(read, 'reference_id', None) is not None and read.reference_id >= 0:
                            try:
                                chrom = bam.get_reference_name(read.reference_id)
                            except Exception:
                                chrom = None
                        ref_cursor = getattr(read, 'reference_start', 0) or 0
                        cigarstr = getattr(read, 'cigarstring', None) or "."
                        read_len = _get_read_length(read)
                        for (op, length) in getattr(read, 'cigartuples', []):
                            if op == 1 and length >= ins_thresh:
                                _write_candidate_line(fout, chrom or ".", ref_cursor, ref_cursor, qn, "INS", length, "SALVAGED_CIGAR_I", read_len, cigarstr)
                                candidates.append({'chrom': chrom or ".", 'start': ref_cursor, 'end': ref_cursor, 'read_name': qn, 'sv_type': 'INS', 'sv_len': length})
                                counts['INS'] += 1
                            if op == 2 and length >= del_thresh:
                                _write_candidate_line(fout, chrom or ".", ref_cursor, ref_cursor + length - 1, qn, "DEL", length, "SALVAGED_CIGAR_D", read_len, cigarstr)
                                candidates.append({'chrom': chrom or ".", 'start': ref_cursor, 'end': ref_cursor + length - 1, 'read_name': qn, 'sv_type': 'DEL', 'sv_len': length})
                                counts['DEL'] += 1
                            if op == 4 and length >= sc_thresh:
                                _write_candidate_line(fout, chrom or ".", ref_cursor, ref_cursor, qn, "SOFTCLIP", length, "SALVAGED_CIGAR_S", read_len, cigarstr)
                                candidates.append({'chrom': chrom or ".", 'start': ref_cursor, 'end': ref_cursor, 'read_name': qn, 'sv_type': 'SOFTCLIP', 'sv_len': length})
                                counts['SOFTCLIP'] += 1
                            if op in OP_CONSUME_REF:
                                ref_cursor += length
                except Exception as e2:
                    elog.write(f"{qn}\tsalvage_failed\t{str(e2)}\n")
            continue

    # write tiny header summary inside the candidates TSV for quick inspection
    fout.write(f"#SUMMARY_total_examined={total_examined}\n")
    fout.write(f"#SUMMARY_total_candidates={len(candidates)}\n")
    fout.write(f"#SUMMARY_total_skipped={skipped}\n")
    for k in ['INS', 'DEL', 'SOFTCLIP', 'LARGE_N', 'SPLIT']:
        fout.write(f"#SUMMARY_{k}={counts.get(k,0)}\n")

    fout.close()
    elog.close()
    bam.close()

    # clustering step (optional)
    clusters = []
    if cluster_window and len(candidates) > 0:
        clusters = cluster_candidates(candidates, cluster_window=cluster_window, min_support=min_support)

    # annotate clusters with GTF if available
    if gtf_exons and clusters:
        for c in clusters:
            genes = gtf_find_genes_for_interval(gtf_exons, c['chrom'], c['start'], c['end'])
            c['genes'] = sorted(list(genes)) if genes else []

    # write clusters to files
    clusters_tsv = out_prefix + ".clusters.tsv"
    clusters_txt = out_prefix + ".clusters.txt"
    header = "chrom\tstart\tend\ttype\tsupport\tmedian_sv_len\treads\tgenes\n"
    with open(clusters_tsv, "w") as cf:
        cf.write(header)
        for c in clusters:
            reads_str = ",".join(c['reads'])
            genes_str = ",".join(c.get('genes', []))
            cf.write(f"{c['chrom']}\t{c['start']}\t{c['end']}\t{c['sv_type']}\t{c['support']}\t{c['median_sv_len']}\t{reads_str}\t{genes_str}\n")
    with open(clusters_txt, "w") as cf2:
        cf2.write(header)
        for c in clusters:
            reads_str = ",".join(c['reads'])
            genes_str = ",".join(c.get('genes', []))
            cf2.write(f"{c['chrom']}\t{c['start']}\t{c['end']}\t{c['sv_type']}\t{c['support']}\t{c['median_sv_len']}\t{reads_str}\t{genes_str}\n")

    # write a simple machine- and human-readable summary file
    summary_path = out_prefix + ".summary.txt"
    summary_obj = {
        'total_examined': total_examined,
        'total_candidates': len(candidates),
        'total_skipped': skipped,
        'counts': dict(counts),
        'clusters_found': len(clusters)
    }
    with open(summary_path, "w") as sf:
        sf.write(json.dumps(summary_obj, indent=2))
        sf.write("\n\nHuman-readable counts by type:\n")
        total_cand = len(candidates) or 1
        for k, v in counts.most_common():
            pct = float(v) / total_cand * 100.0 if total_cand else 0.0
            sf.write(f"{k}: {v} ({pct:.1f}%)\n")

    # console summary (friendly)
    print(f"[SUMMARY] total_examined={total_examined}, per-read_candidates={len(candidates)}, clusters_found={len(clusters)}, total_skipped={skipped}")
    print("[SUMMARY] counts by type:")
    for k, v in counts.most_common():
        pct = float(v) / (len(candidates) or 1) * 100.0
        print(f"   {k}: {v} ({pct:.1f}%)")

    if skipped:
        print(f"[WARN] See {err_log} for parse/skipped read errors.")
    if clusters:
        print(f"[INFO] clusters written to: {clusters_tsv} and {clusters_txt}")
    print(f"[INFO] summary written to: {summary_path}")

    return candidates, clusters


def main():
    p = argparse.ArgumentParser(description="iso_parser.py: parse BAM for SV-like signals, cluster, and optionally annotate with GTF")
    p.add_argument("bam", help="input sorted BAM (indexed recommended)")
    p.add_argument("out", help="output TSV of per-read candidates (filename; placed into --out-dir if given)")
    p.add_argument("--out-dir", help="optional directory to place all outputs (created if missing)")
    p.add_argument("--mq", type=int, default=20, help="mapping quality threshold")
    p.add_argument("--sample", type=int, default=5000, help="how many reads to sample for type inference")
    p.add_argument("--readlen-cutoff", type=int, default=1000, help="median read length cutoff to separate short vs long")
    p.add_argument("--insert", type=int, default=None, help="min insertion size to report (override auto)")
    p.add_argument("--delete", type=int, default=None, help="min deletion size to report (override auto)")
    p.add_argument("--softclip", type=int, default=None, help="min softclip length to report (override auto)")
    p.add_argument("--large-N", type=int, default=None, help="report very large N (skipped region) as candidate if >= value")
    p.add_argument("--sample-parse", type=int, default=None, help="only parse this many reads (for quick tests)")
    p.add_argument("--salvage", action="store_true", help="try to salvage minimal info from reads that raised parse errors (conservative)")
    p.add_argument("--cluster-window", type=int, default=None, help="bp window to cluster per-read candidates (e.g. 50)")
    p.add_argument("--min-support", type=int, default=1, help="min distinct supporting reads for a cluster to be kept")
    p.add_argument("--gtf", type=str, default=None, help="optional GTF (Gencode v43) to annotate clusters (can be .gz)")
    args = p.parse_args()

    if not os.path.exists(args.bam):
        print(f"[ERROR] BAM file not found: {args.bam}", file=sys.stderr)
        sys.exit(2)

    header_type = infer_bam_type_from_header(args.bam)
    if header_type:
        print(f"[INFO] Header RG platform suggests BAM type = {header_type}")
    else:
        print("[INFO] No clear RG PL tag found in header.")

    sampled_type, median_len, mean_len, n_sampled, n_skipped_sampling = infer_bam_type_by_readlen(
        args.bam, sample_reads=args.sample, readlen_cutoff=args.readlen_cutoff)
    if sampled_type:
        print(f"[INFO] Sampled reads: n_sampled={n_sampled}, median_length={median_len}, mean_length={mean_len:.1f}. Inferred type = {sampled_type}")
    else:
        print(f"[WARN] Could not sample reads (n_sampled={n_sampled}); sampling skipped or failed. n_skipped_sampling={n_skipped_sampling}")

    final_type = header_type or sampled_type or "unknown"
    print(f"[INFO] Final BAM classification: {final_type}")

    short_defaults = (30, 30, 30)
    long_defaults = (200, 200, 200)
    ins_def, del_def, sc_def = choose_thresholds("long" if final_type == "long" else "short", short_defaults, long_defaults)

    ins_thresh = args.insert if args.insert is not None else ins_def
    del_thresh = args.delete if args.delete is not None else del_def
    sc_thresh = args.softclip if args.softclip is not None else sc_def

    print(f"[INFO] Using thresholds: INS={ins_thresh}, DEL={del_thresh}, SOFTCLIP={sc_thresh}, MQ={args.mq}, SALVAGE={args.salvage}")

    gtf_exons = None
    if args.gtf:
        if not os.path.exists(args.gtf):
            print(f"[ERROR] GTF file not found: {args.gtf}", file=sys.stderr)
            sys.exit(2)
        print(f"[INFO] Loading GTF (may take a while): {args.gtf}")
        gtf_exons = parse_gtf_exons(args.gtf)
        print("[INFO] GTF loaded (exons parsed)")

    parse_bam_and_cluster(args.bam, args.out, final_type, median_len, mean_len,
                          mq_thresh=args.mq, ins_thresh=ins_thresh, del_thresh=del_thresh,
                          sc_thresh=sc_thresh, large_N_thresh=args.large_N, sample_limit=args.sample_parse,
                          salvage=args.salvage, cluster_window=args.cluster_window, min_support=args.min_support,
                          gtf_exons=gtf_exons, out_dir=args.out_dir)


if __name__ == "__main__":
    main()

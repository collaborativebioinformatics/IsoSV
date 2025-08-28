#!/usr/bin/env python3

## This script builds BED files for candidate SV search, and builds intervaltree for transcript-level SV annotation queries.
## Contact: Louis SHE (snakesch@connect.hku.hk)

import pandas as pd
import re
import pickle

from intervaltree import Interval, IntervalTree

## Helpers 
def expand_gtf_attributes(df):
    """
    Expand a GTF attribute column into separate columns, assuming no repeated keys.
    """
    pair_re = re.compile(r'\s*([^;\s]+)\s+("([^"]*)"|[^;"]+)\s*;')

    def parse_row(attr: str) -> dict:
        if not isinstance(attr, str):
            return {}
        d: dict[str, str] = {}
        for m in pair_re.finditer(attr):
            k = m.group(1)
            raw = m.group(2).strip()
            v = raw[1:-1] if len(raw) >= 2 and raw[0] == '"' and raw[-1] == '"' else raw
            # keep only the first occurrence by spec/assumption
            if k not in d:
                d[k] = v
        return d

    attrs_df = df["annotations"].apply(parse_row).apply(pd.Series)

    attrs_df["level"] = pd.to_numeric(attrs_df["level"], errors="coerce").astype("Int64")

    out = pd.concat(
        [df.drop(columns=["annotations"]).reset_index(drop=True),
         attrs_df.reset_index(drop=True)],
        axis=1
    )
    return out

cleaned_gtf = "../resources/cleaned_annotations.gtf"

gtf = pd.read_csv(cleaned_gtf, sep="\t", header=None)
gtf.columns = ["chrom", "type", "start", "end", "strand", "annotations"]

## - 1. Build gene BED for SV candidate search - ##

## First construct a gene table
df = gtf.loc[gtf["type"] == "gene", ["chrom", "start", "end", "strand", "annotations"]]

## HGNC ID can be null, so we use ENSG gene ID as hashes
gene_table = expand_gtf_attributes(df).loc[:, ['chrom', 'start', 'end', 'strand', 'gene_id', 'gene_type', 'gene_name']]

## Write gene regions for candidate SVs
gene_table.to_csv("../resources/gencode.v43.annotation.genes.bed", sep="\t", index=False, header=False)

## - 2. Build intervaltree for transcripts - ##
df = gtf.loc[gtf["type"] != "gene", ["chrom", "start", "end", "strand", "type", "annotations"]]
tx_df = expand_gtf_attributes(df).loc[:, ['chrom', 'start', 'end', 'strand', 'type', 'gene_name', 'transcript_name', 'transcript_id', 'gene_type', 'transcript_type', 'tag', 'exon_number']]

tx_tree = {}
for tx_id, row_df in tx_df.groupby('transcript_id', sort=False):
    # Get the single transcript row as a Series
    tx_row = row_df.loc[row_df["type"] == "transcript"].iloc[0]

    chrom = str(tx_row["chrom"])
    tx_start = int(tx_row["start"])
    tx_end = int(tx_row["end"])
    gene_name = str(tx_row["gene_name"])
    transcript_name = str(tx_row["transcript_name"])
    transcript_id = str(tx_row["transcript_id"])
    gene_type = (None if "gene_type" not in tx_row or pd.isna(tx_row["gene_type"]) else str(tx_row["gene_type"]))
    transcript_type = (None if "transcript_type" not in tx_row or pd.isna(tx_row["transcript_type"]) else str(tx_row["transcript_type"]))
    strand = (None if "strand" not in tx_row or pd.isna(tx_row["strand"]) else str(tx_row["strand"]))
    tag = (None if "tag" not in tx_row or pd.isna(tx_row["tag"]) else str(tx_row["tag"]))

    if chrom not in tx_tree:
        tx_tree[chrom] = IntervalTree()

    # Build a pure-Python payload (no pandas objects)
    payload = [
        {
            "region": "intronic",
            "gene": gene_name,
            "gene_type": gene_type,
            "transcript_name": transcript_name,
            "transcript_id": transcript_id,
            "transcript_type": transcript_type,
            "strand": strand,
            "tag": tag,
        }
    ]

    # Prefer addi to avoid any special __setitem__ slicing edge cases; both are equivalent once types are clean
    tx_tree[chrom].addi(tx_start, tx_end + 1, payload)

    # Append elements into the same payload list
    for _, row in row_df.iterrows():
        rtype = row["type"]
        if rtype == "transcript":
            continue

        elem_start = int(row["start"])
        elem_end = int(row["end"])
        exon_num_val = row.get("exon_number", None)
        exon_number = None if pd.isna(exon_num_val) else int(exon_num_val)

        # Find the transcript interval by a point lookup and matching bounds
        holder = tx_tree[chrom][tx_start]
        for iv in holder:
            if iv.begin == tx_start and iv.end == tx_end + 1:
                iv.data.append({
                    "region": str(rtype) + str(exon_number),
                    "start": elem_start,
                    "end": elem_end,
                    "gene": str(row.get("gene_name", gene_name)) if not pd.isna(row.get("gene_name", gene_name)) else gene_name,
                    "transcript_name": str(row.get("transcript_name", transcript_name)) if not pd.isna(row.get("transcript_name", transcript_name)) else transcript_name,
                    "transcript_id": str(row.get("transcript_id", transcript_id)) if not pd.isna(row.get("transcript_id", transcript_id)) else transcript_id,
                    "strand": (None if pd.isna(row.get("strand", strand)) else str(row.get("strand", strand))),
                    "gene_type": (None if pd.isna(row.get("gene_type", gene_type)) else str(row.get("gene_type", gene_type))),
                    "transcript_type": (None if pd.isna(row.get("transcript_type", transcript_type)) else str(row.get("transcript_type", transcript_type))),
                })
                break
            
outpath = "../resources/tx_tree_cache.pkl"
with open(outpath, 'wb') as f:
    pickle.dump(tx_tree, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"Saved tx_tree to {outpath}")
#!/usr/bin/env python3

## This script builds BED files for candidate SV search, and builds intervaltree for transcript-level SV annotation queries.
## Contact: Louis SHE (snakesch@connect.hku.hk)

import pandas as pd
import argparse
import os

## Transcript tree usage:
## tx_tree["chr21"][39217043:39217093]
## {Interval(39215363, 39218152, ['BRWD1:ENST00000380800.7:intron_31', 'BRWD1-204', 'protein_coding', 'protein_coding']),
##  Interval(39216784, 39217401, ['TIMM9P2:ENST00000419324.1:intron_1', 'TIMM9P2-201', 'processed_pseudogene', 'processed_pseudogene'])}

def expand_gtf_attributes(df):
    """
    Expand a GTF attribute column into separate columns, assuming no repeated keys.
    """
    import re
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

def preprocess_gtf(gtf_path) -> pd.DataFrame:
    '''Preprocess GTF downloaded from GENCODE and returns a filtered GTF as pandas df'''
    ## 2025-08-29: no longer generates cleaned_annotations.gtf, original file kept for debugging
    import subprocess
    import os
    
    if not os.path.exists(gtf_path):
        raise FileNotFoundError("GTF file not found! ")
    
    cmd = f"grep -wv CDS {gtf_path} | grep -w 'level 1' | grep -vw TEC | cut -f1,4,5,3,7,9- "
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    df = pd.read_csv(p.stdout, sep="\t", header=None, names = ["chrom", "type", "start", "end", "strand", "annotations"])

    if df.empty:
        raise ValueError("Input GTF file may be corrupt. Check input. ")
    
    ret = p.wait()
    if ret != 0:
        raise RuntimeError(f"pipeline failed ({ret}): {p.stderr.read()}")
    
    return df 

def build_gene_bed(gtf: pd.DataFrame, outpath: str = "../resources/gencode.v43.annotation.genes.bed") -> None:
    '''Takes a preprocessed GTF and extracts gene coords'''
    ## First construct a gene table
    df = gtf.loc[gtf["type"] == "gene", ["chrom", "start", "end", "strand", "annotations"]]

    ## HGNC ID can be null, so we use ENSG gene ID as hashes
    gene_table = expand_gtf_attributes(df).loc[:, ['chrom', 'start', 'end', 'strand', 'gene_id', 'gene_type', 'gene_name']]

    ## Write gene regions for candidate SVs
    gene_table.to_csv(outpath, sep="\t", index=False, header=False)

def build_transcript_tree(full_gtf) -> dict:
    '''Takes the complete preprocessed GTF and returns a transcript tree (tx_tree) for later queries'''
    
    from intervaltree import Interval, IntervalTree
    
    df = full_gtf.loc[full_gtf["type"] != "gene", ["chrom", "start", "end", "strand", "type", "annotations"]]
    tx_df = expand_gtf_attributes(df).loc[:, ['chrom', 'start', 'end', 'strand', 'type', 'gene_name', 'transcript_name', 'transcript_id', 'gene_type', 'transcript_type', 'tag', 'exon_number']]

    tx_tree = {}
    for transcript_id, row_df in tx_df.groupby('transcript_id', sort=False):
        # Get the single transcript row as a Series
        tx_row = row_df.loc[row_df["type"] == "transcript"].iloc[0]

        chrom = str(tx_row["chrom"])
        tx_start = int(tx_row["start"])
        tx_end = int(tx_row["end"])
        gene_name = str(tx_row["gene_name"])
        transcript_name = str(tx_row["transcript_name"])
        gene_type = (None if "gene_type" not in tx_row or pd.isna(tx_row["gene_type"]) else str(tx_row["gene_type"]))
        transcript_type = (None if "transcript_type" not in tx_row or pd.isna(tx_row["transcript_type"]) else str(tx_row["transcript_type"]))
        strand = tx_row["strand"]

        if chrom not in tx_tree:
            tx_tree[chrom] = IntervalTree()

        row_df = row_df[~row_df["type"].isin(["transcript", "start_codon", "stop_codon", "Selenocysteine"])] ## Make sure only exon and UTR in type column
        
        ## make sure exon_df is sorted in ascending order of start coords
        exon_df = row_df[row_df["type"] == "exon"].sort_values("start")
        exon_df["intron_start"] = exon_df["end"] + 1
        exon_df["intron_end"] = ((exon_df["start"] - 1).shift(-1).fillna(-1)).astype(int)

        for _, row in exon_df.iterrows():
            rtype = row["type"]
            elem_start = row["start"]
            elem_end = row["end"]
            exon_number = int(row["exon_number"])
            
            ## Exonic regions
            tx_tree[chrom][elem_start:elem_end+1] = [
                "{}:{}:{}".format(gene_name, transcript_id, str(rtype) + "_" + str(exon_number)), transcript_name, gene_type, transcript_type
            ]
            
            ## Intronic regions
            intron_start = row["intron_start"]
            intron_end = row["intron_end"]
            
            if intron_end > 0 and intron_start > 0:
                intron_number = exon_number if strand == "+" else exon_number - 1
                tx_tree[chrom][intron_start:intron_end+1] = [
                "{}:{}:{}".format(gene_name, transcript_id, "intron_" + str(intron_number)), transcript_name, gene_type, transcript_type
                ]

        
        for _, row in row_df[row_df["type"] == "UTR"].iterrows():
            rtype = row["type"]
            elem_start = row["start"]
            elem_end = row["end"]
            number = int(row["exon_number"])
            
            ## UTR
            tx_tree[chrom][elem_start:elem_end+1] = [
                "{}:{}:{}".format(gene_name, transcript_id, str(rtype) + "_" + str(number)), transcript_name, gene_type, transcript_type
            ]
                
    return tx_tree


def pickle_tree(tree, outpath = "../resources/tx_tree_cache.pkl"):
    import pickle
    with open(outpath, 'wb') as f:
        pickle.dump(tree, f, protocol=pickle.HIGHEST_PROTOCOL)
        print(f"Saved tx_tree to {outpath}")

def main():
    parser = argparse.ArgumentParser(description="Prepare gene coordinates and builds transcript tree. ")
    parser.add_argument("--gtf", required=True, help="Path to downloaded GTF file")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory for gene BED and cache")
    args = parser.parse_args()
    
    ## gtf_path = "../resources/gencode.v43.annotation.gtf" ## This file is large ! 
    gtf = preprocess_gtf(gtf_path = args.gtf)
    print(f"Done preprocessing GTF from {args.gtf} ")

    build_gene_bed(gtf, outpath=os.path.join(args.outdir, 'gencode.v43.annotation.genes.bed'))

    tx_tree = build_transcript_tree(gtf)
    pickle_tree(tx_tree, outpath=os.path.join(args.outdir, 'tx_tree_cache.pkl'))
    print(f"Gene BED and transcript tree built successfully. ")

if __name__ == "__main__":
    main()
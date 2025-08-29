#!/usr/bin/env python3

## The plan is to define a single interval query and then iterate over all candidates from step b.

# Test file from IsoSV/step_b_IsoClustering/test_data/chr21_SVs_converted.tsv:
# chrom	start	end	    type	support	median_sv_len	reads	genes
# 21	5227407	5227477	DEL	    1	    70		
# 21	5252520	5255771	DEL	    1	    3251		
# 21	5276291	5276291	INS	    1	    2560	

import os
import argparse
import pandas as pd
from intervaltree import Interval, IntervalTree

cache = "/data/bnf/dev/saile/other/prj/hacathon_RNASV/resources/tx_tree_cache.pkl"
test_data = "../../step_b_IsoClustering/test_data/chr21_SVs_converted.tsv"

def load_tx_tree(filepath):
    """Load the IntervalTree from disk"""
    import pickle
    if not os.path.exists(filepath):
        raise FileNotFoundError("Transcript tree not found! ")
    with open(filepath, 'rb') as f:
        return pickle.load(f)

def annotate_candidates(candidates: pd.DataFrame, tx_tree):
    
    results = []

    for _, candidate in candidates.iterrows():
        chrom, start, end, svtype, support, median_svlen, *misc = candidate
        
        if not str(chrom).startswith("chr"):
            chrom = "chr" + str(chrom)
        
        overlapping_transcripts = tx_tree[chrom][start:end]
        
        if not overlapping_transcripts: ## Case: purely intergenic
            feature = "intergenic"
            results.append((chrom, start, end, svtype, support, median_svlen, "intergenic", "na", "na", "na"))
        elif len(overlapping_transcripts) == 1:
            ## simple SVs
            feature = next(iter(overlapping_transcripts)).data
            results.append((chrom, start, end, svtype, support, median_svlen, feature[0], feature[1], feature[2], feature[3]))
        else:
            ## fusion?
            feature0, feature1, feature2, feature3 = "", "", "", ""
            
            ## If a single transcript is involved, no fusion induced by deletion
            involved_tx = [iv.data[0].split(":")[1] for iv in overlapping_transcripts]
            
            if len(involved_tx) == 1:
                for tx in overlapping_transcripts:
                    feature0 += tx.data[0] + ";"
                    feature1 += tx.data[1] + ";"
                    feature2 += tx.data[2] + ";"
                    feature3 += tx.data[3] + ";"
                results.append((chrom, start, end, svtype, support, median_svlen, feature0, feature1, feature2, feature3))
            else:
                ## candidate fusion induced by long deletion
                for tx in overlapping_transcripts:
                    feature0 += tx.data[0] + ";"
                    feature1 += tx.data[1] + ";"
                    feature2 += tx.data[2] + ";"
                    feature3 += tx.data[3] + ";"
                results.append((chrom, start, end, "FUSION", support, median_svlen, feature0, feature1, feature2, feature3))
    return results

def write_to_vcf(annotated_df, outpath):
    
    vcf_out = open(outpath, "w")

    ## Need to add back the contig headers (VCF is still a bit problematic here ...)
    header_lines = [
        "##fileformat=VCFv4.2",
        "##source=IsoSV",
        "##FILTER=<ID=PASS,Description=\"All filters passed\">",
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">",
        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">",
        "##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description=\"Support count for this SV\">",
        "##INFO=<ID=REGION,Number=1,Type=String,Description=\"Gene region\">",
        "##INFO=<ID=TX_NAME,Number=1,Type=String,Description=\"Transcript name\">",
        "##INFO=<ID=GENE_BIOTYPE,Number=1,Type=String,Description=\"Gene biotype\">",
        "##INFO=<ID=TX_BIOTYPE,Number=1,Type=String,Description=\"Transcript biotype\">",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    ]

    for header in header_lines:
        vcf_out.write(header)
        vcf_out.write("\n")

    for _, row in annotated_df.iterrows():
        chrom, start, stop, svtype, support, svlen, region, tx_alias, biotype_gene, biotype_tx = row

        record = f"{chrom}\t{start}\t.\tN\t<{svtype.upper()}>\t.\tPASS\tEND={stop};SVTYPE={svtype};SVLEN={svlen};SUPPORT={support};REGION={region}"
        
        if tx_alias != "na":
            record += fr";TX_NAME={tx_alias}"
        if biotype_gene != "na":
            record += fr";GENE_BIOTYPE={biotype_gene}"
        if biotype_tx != "na":
            record += fr";TX_BIOTYPE={biotype_tx}"

        record += "\n"
        vcf_out.write(record)
        
    vcf_out.close()

def main():
    parser = argparse.ArgumentParser(description="Annotate SV candidates using prebuilt transcript trees. ")
    parser.add_argument("--candidates", required=True, help="Path to SV candidates") ## currently adapted to step b chr21 test output
    parser.add_argument("--cache", required=True, help="Transcript tree file path")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory for annotated VCF")
    args = parser.parse_args()
    
    
    tx_tree = load_tx_tree(args.cache)

    candidates = pd.read_csv(args.candidates, sep="\t")
    annotated = annotate_candidates(candidates, tx_tree)

    ## Write to VCF
    out = pd.DataFrame(annotated, columns=["chr", "start", "stop", "SVTYPE", "SUPPORT", "SVLEN", "REGION", "TX_ALIAS", "BIOTYPE_GENE", "BIOTYPE_TX"])
    write_to_vcf(out, os.path.join(args.outdir, "sv_candidates.annotated.vcf"))
    print(f"Result VCF written to {os.path.join(args.outdir, 'sv_candidates.annotated.vcf')}")
    
if __name__ == "__main__":
    main()
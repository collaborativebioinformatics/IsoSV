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

    for _, row in candidates.iterrows():
        # Read columns by name to be robust to input file format changes
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
        
        # The query region is the full span for DEL, and a single point for INS
        query_start = start
        query_end = end if svtype == 'DEL' else start + 1
        overlapping_hits = tx_tree[chrom][query_start:query_end]
        
        # 1. Handle Intergenic case
        if not overlapping_hits:
            results.append((chrom, start, end, svtype, "Intergenic", support, median_svlen, "Intergenic", "na", "na", "na"))
            continue

        # 2. Collect detailed information from all hits
        touched_genes = set()
        hit_exons = []
        hit_introns = []
        
        all_tx_aliases = set()
        all_gene_biotypes = set()
        all_tx_biotypes = set()

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

        # 3. Apply classification logic based on SVTYPE
        final_annotation = "Unclassified" # Default

        if svtype == 'DEL':
            if len(touched_genes) > 1:
                final_annotation = "Gene_Fusion"
            elif len(touched_genes) == 1:
                final_annotation = "Intronic/Novel_Splicing"
                is_exonic_deletion = any(start <= exon.begin and end >= exon.end for exon in hit_exons)
                
                if is_exonic_deletion:
                    final_annotation = "Exonic_Deletion"
                else:
                    is_canonical_splicing = any(abs(start - intron.begin) <= 5 and abs(end - intron.end) <= 5 for intron in hit_introns)
                    if is_canonical_splicing:
                        final_annotation = "Canonical_Splicing"
        elif svtype == 'INS':
            if len(touched_genes) > 1:
                final_annotation = "Insertion_in_Fusion_Region"
            elif len(touched_genes) == 1:
                final_annotation = "Intragenic_Insertion"

        # 4. Format results for VCF output
        region_field = ",".join(sorted(list(touched_genes))) if touched_genes else "na"
        tx_alias_field = ",".join(sorted(list(all_tx_aliases))) or "na"
        gene_biotype_field = ",".join(sorted(list(all_gene_biotypes))) or "na"
        tx_biotype_field = ",".join(sorted(list(all_tx_biotypes))) or "na"

        results.append((
            chrom, start, end, svtype,
            final_annotation,
            support, median_svlen,
            region_field, tx_alias_field,
            gene_biotype_field, tx_biotype_field
        ))
            
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
        "##INFO=<ID=ANN_TYPE,Number=1,Type=String,Description=\"Functional annotation of the SV (e.g., Gene_Fusion, Exonic_Deletion)\">",
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
        chrom, start, stop, svtype, ann_type, support, svlen, region, tx_alias, biotype_gene, biotype_tx = row

        record = f"{chrom}\t{start}\t.\tN\t<{svtype.upper()}>\t.\tPASS\tEND={stop};SVTYPE={svtype};SVLEN={svlen};SUPPORT={support};ANN_TYPE={ann_type};REGION={region}"
        
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
    out = pd.DataFrame(annotated, columns=["chr", "start", "stop", "SVTYPE", "ANN_TYPE", "SUPPORT", "SVLEN", "REGION", "TX_ALIAS", "BIOTYPE_GENE", "BIOTYPE_TX"])
    outpath = os.path.join(args.outdir, "GIAB002_chr22_region_LongReadSV.annotated.vcf")
    write_to_vcf(out, outpath)
    print(f"Result VCF written to {outpath}")
    
if __name__ == "__main__":
    main()
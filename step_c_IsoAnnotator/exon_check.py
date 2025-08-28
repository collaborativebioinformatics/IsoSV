#!/usr/bin/env python3
import argparse, csv, gzip, io, re, os


def open_text(path):
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"))
    return open(path, "r", encoding="utf-8")

def overlap(a1, a2, b1, b2):
    return not (a2 < b1 or a1 > b2)

def fully_within(a1, a2, b1, b2):
    return a1 >= b1 and a2 <= b2

# minimal regex for GTF attributes
ATTR_TX = re.compile(r'(?:^|;) *transcript_id *"([^"]+)"')
ATTR_GENE = re.compile(r'(?:^|;) *gene_id *"([^"]+)"')

def parse_attrs(attr_field):
    tx = ATTR_TX.search(attr_field)
    gene = ATTR_GENE.search(attr_field)
    return (tx.group(1) if tx else None, gene.group(1) if gene else None)

def load_exons(gtf_path):
    tx_exons = {}
    tx2gene = {}
    # go through the gtf fike and check and parse for exons and genes
    with open_text(gtf_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"): 
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9: 
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = fields
            if feature.lower() != "exon":
                continue
            try:
                start, end = int(start), int(end)
            except: 
                continue
            tx, gene = parse_attrs(attrs)
            if not tx:
                continue
            key = (chrom, strand, tx)
            tx_exons.setdefault(key, []).append((start,end))
            if gene: 
                tx2gene[tx] = gene
    # sort exons
    for k in tx_exons:
        tx_exons[k].sort()
    return tx_exons, tx2gene

## annotate whether the region regions are in intergenic, exonic or partial_exonic region
def classify_region(chrom, start, end, tx_exons, tx2gene):
    full, partial = False, False
    for (c, strand ,tx), exons in tx_exons.items():
        if c!=chrom: 
            continue
        for (s,e) in exons:
            if fully_within(start,end,s,e):
                full = True
            elif overlap(start,end,s,e):
                partial = True
    if full: 
        return "exon"
    if partial: 
        return "exon_partial"
    return "intergenic"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a","--annotation", required=True, help="GTF file (.gtf or .gtf.gz)")
    parser.add_argument("-i","--input", required=True, help="CSV with columns: Chr,Start,End,Readname,Type")
    parser.add_argument("-o","--output_file", type=str, required= True, help="Path to save the output")
    args = parser.parse_args()

    tx_exons, tx2gene = load_exons(args.annotation)
    reader = csv.DictReader(open(args.input))
    fieldnames = reader.fieldnames + ["Summary"]
    os.makedirs(os.path.dirname(args.output_file), exist_ok=True)

    with open(args.output_file, 'w') as output_file:
        writer = csv.DictWriter(output_file, fieldnames=fieldnames)
        writer.writeheader()

        for row in reader:
            chrom = row["Chr"].strip()
            start = int(row["Start"])
            end   = int(row["End"])
            if start > end: 
                start,end = end,start
            summary = classify_region(chrom, start, end, tx_exons, tx2gene)
            row["Summary"] = summary
            writer.writerow(row)

if __name__=="__main__":
    main()

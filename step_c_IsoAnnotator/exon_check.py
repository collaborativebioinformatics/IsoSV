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
    """
    Load exon coordinates from a GTF file

    Args:
        gtf_path (str): Path for the gtf file . If you already cleaned a file that is great ;)

    Returns:
        tuple:
            - dict: {(chrom, strand, transcript_id): [(start, end), ...]}  
              A dictionary mapping each transcript to its list of exon coordinates (sorted).
            - dict: {transcript_id: gene_id}  
              A mapping from transcript IDs to their associated gene IDs.
    """
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
    """
    Classify a query interval relative to exons.

    Args:
        chrom (str): Chromosome name of the query interval.
        start (int): Start coordinate of the query (1-based inclusive).
        end (int): End coordinate of the query (1-based inclusive).
        tx_exons (dict): Transcript → exon coordinates dictionary from load_exons().

    Returns:
        str: One of:
            - "exon"         → interval fully contained in at least one exon
            - "exon_partial" → interval overlaps an exon but is not fully contained
            - "intergenic"   → no overlap with any exon
    """
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

## annotate the splice sites 
def find_splice_hits(chrom, start, end, tx_exons, tx2gene, window):
    """
    Identify splice-site proximity hits around exon boundaries.
    """
    if window <= 0:
        return ("no", [])
    
    hits = []

    for (c, _strand, tx), exons in tx_exons.items():
        if c!= chrom:
            continue
        for i, (s,e) in enumerate(exons, start=1 ):
            # window around exon start
            ws1, we1 = s - window, s+ window
            if overlap(start, end, ws1, we1):
                hits.append(f"{tx}|{tx2gene.get(tx,'-')}|exon{i}|start|{s}")
            # window around exon end
            ws2, we2 = e - window, e + window
            if overlap(start, end, ws2, we2):
                hits.append(f"{tx}|{tx2gene.get(tx,'-')}|exon{i}|end|{e}")
    return ("yes" if hits else "no", hits)
            
        
## write results into vcf file
def write_vcf_header(out):
    out.write("##fileformat=VCFv4.3\n")
    out.write("##source=indel_annotation_splice\n")
    out.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
    out.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the structural variant">\n')
    out.write('##INFO=<ID=SUMMARY,Number=1,Type=String,Description="Exon overlap summary (exon|exon_partial|intergenic)">\n')
    out.write('##INFO=<ID=SPLICE,Number=1,Type=String,Description="yes/no if region overlaps ±N bp of any exon boundary">\n')
    out.write('##INFO=<ID=HITS,Number=.,Type=String,Description="Splice hits: TX|GENE|exonN|start/end|coord">\n')
    out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

## Calculate the SV length
def make_svlen(row_svtype, start, end, provided_svlen):
    """
    Prefer provided SVLen; otherwise compute from span with DEL negative, INS positive.
    """
    if provided_svlen not in (None, "", "NA", "NaN"):
        try:
            return int(provided_svlen)
        except:
            pass
    span = end - start + 1
    t = (row_svtype or "").upper()
    if t == "DEL":
        return -abs(span)
    if t == "INS":
        return abs(span)
    return span

def emit_vcf_record(out, row, summary, splice_flag, splice_hits, ref_base="N"):
    chrom = row["Chr"].strip()
    start = int(row["Start"])
    end   = int(row["End"])
    if start > end:
        start, end = end, start

    vid    = row.get("Readname(ID)") or row.get("ID") or "."
    svtype = (row.get("SVType") or "SV").upper()
    svlen  = make_svlen(svtype, start, end, row.get("SVLen"))
    alt    = f"<{svtype}>"

    info_parts = [
        f"SVTYPE={svtype}",
        f"SVLEN={svlen}",
        f"SUMMARY={summary}",
        f"SPLICE={splice_flag}",
    ]
    if splice_hits:
        info_parts.append("HITS=" + ",".join(splice_hits))
    info = ";".join(info_parts)

    out.write(f"{chrom}\t{start}\t{vid}\t{ref_base}\t{alt}\t.\tPASS\t{info}\n")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a","--annotation", required=True, help="GTF file (.gtf or .gtf.gz)")
    parser.add_argument("-i","--input", required=True, help="CSV with columns: Chr,Start,End,Readname,Type")
    parser.add_argument("-o","--output_file", type=str, required= True, help="Path to save the output")
    parser.add_argument("--splice-window", type=int, default=2,help="±N bp around exon edges to flag as splice proximity (default 2; 0 disables).")
    parser.add_argument("--ref-base", default="N", help="REF base for symbolic SVs (default N)")
    args = parser.parse_args()

    tx_exons, tx2gene = load_exons(args.annotation)
    

    out_dir = os.path.dirname(args.output_file)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

        
    with open(args.output_file, "w", encoding="utf-8") as output_file:
        write_vcf_header(output_file)
        reader = csv.DictReader(open(args.input, "r", encoding="utf-8"))
     

        for row in reader:
            if not row:
                continue
            chrom = row["Chr"].strip()
            start = int(row["Start"])
            end   = int(row["End"])
            if start > end: 
                start,end = end,start
            summary = classify_region(chrom, start, end, tx_exons, tx2gene)
            splice_flag, splice_hits = find_splice_hits(chrom, start, end, tx_exons, tx2gene, args.splice_window)


            row["Summary"] = summary
            row["SpliceSiteFlag"] = splice_flag
            row["SpliceHits"] = ";".join(splice_hits)
            emit_vcf_record(output_file, row, summary, splice_flag, splice_hits, ref_base=args.ref_base)

if __name__=="__main__":
    main()


### python exon_check.py -a path/to/resources/gencode.v43.annotation.gtf -i path/to/resources/mock.csv --splice-window 2 -o /path/to/results/ ./results/mock_annotated_regions.vcf
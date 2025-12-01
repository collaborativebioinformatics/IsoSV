
# pysam CIGAR operator codes
M, I, D, N, S, H, P, EQ, X, B = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

def within(a, b, w):
    return abs(a - b) <= w

def length_similar(a, b, abs_slop=5, frac=0.10):
    return abs(a - b) <= max(abs_slop, int(frac * max(a, b)))

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

def write_to_vcf(annotated_candidates, outpath):
    """
    Writes annotated SVs to a VCF file.

    Args:
        annotated_candidates: A list of tuples containing the annotated SVs.
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
            "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO"
        ]
        vcf_out.write("\n".join(header_lines) + "\n")

        for cand in annotated_candidates:
            chrom, start, stop, svtype, ann_type, support, svlen, region, tx_alias = cand
            record = f"{chrom}	{start}	.	N	<{svtype.upper()}>	.	PASS	END={stop};SVTYPE={svtype};SVLEN={svlen};SUPPORT={support};ANN_TYPE={ann_type};REGION={region}"
            if tx_alias != "na": record += f";TX_NAME={tx_alias}"
            vcf_out.write(record + "\n")
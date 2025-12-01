from variant_objects import IndelObs, ClipObs, SplitObs
from utils import is_poly_at

# pysam CIGAR operator codes
M, I, D, N, S, H, P, EQ, X, B = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9

def ref_span_from_cigar(cigar_tuples) -> int:
    """Calculates the span of a CIGAR string on the reference genome."""
    span = 0
    for op, length in cigar_tuples or []:
        if op in (M, EQ, X, D, N):
            span += length
    return span

def walk_cigar_indels_and_clips(aln, min_ins, min_del, min_clip):
    chrom = aln.reference_name
    mapq = aln.mapping_quality
    rname = aln.query_name
    try:
        read_seq = aln.query_sequence
    except Exception:
        read_seq = None

    cigar = aln.cigartuples or []
    ref_pos = aln.reference_start  # 0-based
    read_pos = 0

    # Left clip
    if cigar:
        op0, len0 = cigar[0]
        if op0 in (S, H) and len0 >= min_clip:
            clipseq = (read_seq[:len0] if (op0 == S and read_seq is not None) else "")
            if not is_poly_at(clipseq):
                yield ClipObs(chrom, aln.reference_start + 1, 'L', len0, mapq, rname)

    for op, length in cigar:
        if op in (M, EQ, X):
            ref_pos += length
            read_pos += length
        elif op == I:
            ins_seq = (read_seq[read_pos: read_pos + length] if (length >= min_ins and read_seq is not None) else None)
            if length >= min_ins:
                yield IndelObs(chrom, ref_pos + 1, "INS", length, mapq, rname, ins_seq)
            read_pos += length
        elif op == D:
            if length >= min_del:
                yield IndelObs(chrom, ref_pos + 1, "DEL", length, mapq, rname, None)
            ref_pos += length
        elif op == N:
            ref_pos += length
        elif op == S:
            read_pos += length
        elif op in (H, P, B):
            pass

    # Right clip
    if cigar:
        opn, lenn = cigar[-1]
        if opn in (S, H) and lenn >= min_clip:
            clipseq = (read_seq[-lenn:] if (opn == S and read_seq is not None) else "")
            if not is_poly_at(clipseq):
                yield ClipObs(chrom, aln.reference_end, 'R', lenn, mapq, rname)

###scan for supplenatry alignement for split-read evidence (SV breakpoints)
def parse_sa_tag(sa: str):
    segs = []
    if not sa:
        return segs
    entries = [e for e in sa.strip().split(';') if e]
    for e in entries:
        rname, pos, strand, cig, mapq, nm = e.split(',')
        segs.append({
            "rname": rname,
            "pos": int(pos),
            "strand": strand,
            "cigar": cig,
            "mapq": int(mapq),
            "nm": int(nm)
        })
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
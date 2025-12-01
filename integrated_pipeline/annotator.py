def annotate_candidate(candidate: dict, tx_tree):
    """
    Annotates SV candidates based on their overlap with gene features.

    Args:
        candidate: A dict with annotations
        tx_tree: An IntervalTree containing gene and transcript information.

    Returns:
        A tuple with each tuple containing the annotated fields for a VCF.
    """
    chrom = candidate['chrom']
    start = int(candidate['pos'])
    end = int(candidate['end'])
    svtype = candidate['svtype']
    support = int(candidate['support'])
    median_svlen = int(candidate['length'])
    
    if not str(chrom).startswith("chr"):
        chrom = "chr" + str(chrom)
    
    if chrom not in tx_tree:
        return None
    
    query_start, query_end = (start, end) if svtype == 'DEL' else (start, start + 1)
    overlapping_hits = tx_tree[chrom][query_start:query_end]
    
    if not overlapping_hits:
        result=(chrom, start, end, svtype, "Intergenic", support, median_svlen, "Intergenic", "na", "na", "na")
        return result

    touched_genes, hit_exons, hit_introns = set(), [], []
    transcripts = set()

    for hit in overlapping_hits:
        gene_name = hit.data["gene_name"]
        transcript_id = hit.data["transcript_id"]
        region_type = hit.data["type"]

        # strand = hit.data["strand"]
        # exon_number = hit.data.get("exon_number", "na")
        # intron_number = hit.data.get("intron_number", "na")

        touched_genes.add(gene_name)
        transcripts.add(transcript_id)

        if region_type == "exon":
            hit_exons.append(hit)
        elif region_type == "intron":
            hit_introns.append(hit)
    
    final_annotation = "Unclassified"
    if svtype == 'DEL':
        if len(touched_genes) > 1:
            final_annotation = "Gene_Fusion"
        elif len(touched_genes) == 1:
            is_exonic_deletion = any(start <= exon.begin and end >= exon.end for exon in hit_exons)
            if is_exonic_deletion:
                final_annotation = "Exonic_Deletion"
            else:
                is_canonical_splicing = any(abs(start - intron.begin) <= 5 and abs(end - intron.end) <= 5 for intron in hit_introns)
                final_annotation = "Canonical_Splicing" if is_canonical_splicing else "Intronic/Novel_Splicing"
    elif svtype == 'INS':
        final_annotation = "Insertion_in_Fusion_Region" if len(touched_genes) > 1 else "Intragenic_Insertion"

    region_field = ",".join(sorted(list(touched_genes))) or "na"
    transcript_field = ",".join(sorted(list(transcripts))) or "na"

    result = (
        chrom, start, end, svtype, final_annotation, support, median_svlen,
        region_field, transcript_field
    )

    return result

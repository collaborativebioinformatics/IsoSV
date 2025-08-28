import argparse
import pickle
import csv

# Load the pre-built IntervalTree 
def load_tx_tree(filepath):
    with open(filepath, 'rb') as f:
        return pickle.load(f)

def main():
    # Set up argument parsing 
    parser = argparse.ArgumentParser(description="Annotate SV candidates using a pre-built IntervalTree.")
    parser.add_argument("-i", "--input", required=True, help="Input CSV of SV candidates.")
    parser.add_argument("-t", "--tree", required=True, help="Path to the tx_tree_cache.pkl file.")
    parser.add_argument("-o", "--output", required=True, help="Path for the output annotated CSV file.")
    args = parser.parse_args()

    # Load the IntervalTree data structure
    tx_tree = load_tx_tree(args.tree)

    # Open the input and output files
    with open(args.input, 'r') as infile, open(args.output, 'w', newline='') as outfile:
        reader = csv.DictReader(infile, delimiter='\t')
        fieldnames = reader.fieldnames + ["Annotation", "Gene(s)"]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        # Loop through every SV candidate (the "multi-interval" part)
        for row in reader:
            chrom = "chr" + row["chrom"] # Ensure chrom format matches tree (e.g., 'chr21')
            start = int(row["start"])
            end = int(row["end"])

            # Perform the query for the current interval
            if chrom not in tx_tree:
                continue
            hits = tx_tree[chrom].overlap(start, end)

            # Apply classification logic
            touched_genes = set()
            for hit in hits:
                gene_name = hit.data[0].get("gene")
                if gene_name:
                    touched_genes.add(gene_name)
            
            annotation = "Intergenic" # Default classification if no gene is found
            gene_result = "N/A"
            
            # If the SV overlaps with two or more genes, it's a gene fusion
            if len(touched_genes) > 1:
                annotation = "Gene_Fusion"
                gene_result = ",".join(sorted(list(touched_genes)))
            # If the SV overlaps with one gene, check if it's an exonic deletion or canonical splicing
            elif len(touched_genes) == 1:
                gene_result = list(touched_genes)[0]
                is_exonic_deletion = False
                is_canonical_splicing = False
                
                # Check for exonic deletion
                for hit in hits:
                    for region in hit.data[1:]: # Skip the first summary region
                        if 'exon' in region['region']:
                            exon_start = region['start']
                            exon_end = region['end']
                            if start <= exon_start and end >= exon_end:
                                is_exonic_deletion = True
                                break
                    if is_exonic_deletion:
                        break
                
                # If not a full exonic deletion, check if it's a canonical splice event
                if not is_exonic_deletion:
                    for hit in hits:
                        # Sort exons by start position to identify introns between them
                        exons = sorted([r for r in hit.data[1:] if 'exon' in r['region']], key=lambda x: x['start'])
                        
                        # Check the space between consecutive exons (introns)
                        for i in range(len(exons) - 1):
                            exon1_end = exons[i]['end']
                            exon2_start = exons[i+1]['start']
                            
                            # Check if the SV perfectly matches the intron boundaries (+/- 5bp tolerance)
                            if abs(start - (exon1_end + 1)) <= 5 and abs(end - (exon2_start - 1)) <= 5:
                                is_canonical_splicing = True
                                break
                        if is_canonical_splicing:
                            break

                if is_exonic_deletion:
                    annotation = "Exonic_Deletion"
                elif is_canonical_splicing:
                    annotation = "Canonical_Splicing"
                else:
                    annotation = "Intronic/Novel_Splicing"

            # Write the annotated result to the output file
            row["Annotation"] = annotation
            row["Gene(s)"] = gene_result
            writer.writerow(row)

if __name__ == "__main__":
    main()
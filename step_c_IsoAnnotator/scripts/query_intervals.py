#!/usr/bin/env python3

## TODO: Add wrapper, argparse

import os
import pickle

cache = "../resources/tx_tree_cache.pkl"

def load_tx_tree(filepath):
    """Load the IntervalTree from disk"""
    if os.path.exists(filepath):
        with open(filepath, 'rb') as f:
            return pickle.load(f)
    return None

tx_tree = load_tx_tree(cache)

## Example query: chr1:130000-150000
# tx_tree["chr1"][130000-150000]

## Expected output: 
## {Interval(131025, 134837, [{'region': 'intronic', 'gene': 'CICP27', 'gene_type': 'processed_pseudogene', 'transcript_name': 'CICP27-201', 'transcript_id': 'ENST00000442987.3', 'transcript_type': 'processed_pseudogene', 'strand': '+', 'tag': 'pseudo_consens'}, {'region': 'exon1', 'start': 131025, 'end': 134836, 'gene': 'CICP27', 'transcript_name': 'CICP27-201', 'transcript_id': 'ENST00000442987.3', 'strand': '+', 'gene_type': 'processed_pseudogene', 'transcript_type': 'processed_pseudogene'}])}
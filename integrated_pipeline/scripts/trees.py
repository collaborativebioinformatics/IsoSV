"""
Optimized GENCODE GTF processor for gene annotation and interval tree construction.

Features:
- Filters by level and transcript support level
- Builds IntervalTrees per chromosome for fast lookups
- Memory-efficient streaming GTF parsing
"""

import re
import sys
from collections import defaultdict
from intervaltree import IntervalTree, Interval
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class GTFParser:
    """Memory-efficient GTF parser with streaming line processing."""
    
    # Pre-compiled regex for attribute parsing (faster than repeated compilation)
    ATTR_PATTERN = re.compile(r'(\w+)\s+"?([^";]+)"?')
    
    def __init__(self, gtf_file, max_level=2, max_tsl=3):
        """
        Initialize GTF parser.
        
        Args:
            gtf_file (str): Path to GTF file
            max_level (int): Maximum annotation level (2=high confidence)
            max_tsl (int): Maximum transcript support level
        """
        self.gtf_file = gtf_file
        self.max_level = max_level
        self.max_tsl = max_tsl
        self.genes = {}  # gene_id -> {chrom, start, end, strand, gene_name}
        self.transcripts = defaultdict(list)  # transcript_id -> [exons]
        self.transcript_info = {}  # transcript_id -> {gene_id, chrom, strand}
    
    def parse_attributes(self, attr_string):
        """
        Fast attribute parsing using pre-compiled regex.
        
        Args:
            attr_string (str): GTF attribute field (column 9)
            
        Returns:
            dict: Parsed attributes
        """
        return dict(self.ATTR_PATTERN.findall(attr_string))
    
    def passes_filters(self, attrs):
        """
        Check if feature passes quality filters.
        
        Args:
            attrs (dict): Parsed attributes
            
        Returns:
            bool: True if passes filters
        """
        # Check level (e.g., "level 2")
        level = attrs.get('level', '')
        if level.isnumeric():
            level = int(level)
            if level > self.max_level:
                return False
            
        # Check transcript support level (TSL)
        tsl = attrs.get('transcript_support_level', '')
        if tsl:
            # TSL format: "tsl:1" or just "1"
            if tsl.split(':')[-1].isnumeric():
                tsl_value = int(tsl.split(':')[-1])
                if tsl_value > self.max_tsl:
                    return False

        return True
    
    def stream_parse(self):
        """
        Stream parse GTF file line by line (memory efficient).
        Yields (feature_type, record) tuples.
        """
        try:
            with open(self.gtf_file, 'r') as f:
                for line in f:
                    # Skip header and comment lines
                    if line.startswith('#'):
                        continue
                    
                    line = line.rstrip('\n')
                    if not line:
                        continue
                    
                    fields = line.split('\t')
                    if len(fields) < 9:
                        logger.warning(f"Skipping malformed line: {line[:50]}...")
                        continue
                    chrom, source, feature_type, start, end, score, strand, frame, attributes = fields
                    
                    # Parse attributes
                    attrs = self.parse_attributes(attributes)
                    
                    # Apply filters (skip level 3 and high TSL)
                    if not self.passes_filters(attrs):
                        continue
                    
                    # Skip start/stop codons
                    if feature_type in ('start_codon', 'stop_codon'):
                        continue
                    
                    record = {
                        'chrom': chrom,
                        'start': int(start),
                        'end': int(end),
                        'strand': strand,
                        'feature_type': feature_type,
                        'attrs': attrs
                    }
                    
                    yield feature_type, record
        
        except FileNotFoundError:
            logger.error(f"GTF file not found: {self.gtf_file}")
            sys.exit(1)
        except Exception as e:
            logger.error(f"Error parsing GTF: {e}")
            sys.exit(1)
    
    def index_features(self):
        """
        Index all features from GTF file.
        Builds gene and transcript dictionaries for later use.
        """
        gene_count = 0
        transcript_count = 0
        exon_count = 0
        
        for feature_type, record in self.stream_parse():
            if feature_type == 'gene':
                gene_id = record['attrs'].get('gene_id')
                gene_name = record['attrs'].get('gene_name', gene_id)
                
                self.genes[gene_id] = {
                    'chrom': record['chrom'],
                    'start': record['start'],
                    'end': record['end'],
                    'strand': record['strand'],
                    'gene_name': gene_name
                }
                gene_count += 1
            
            elif feature_type == 'transcript':
                transcript_id = record['attrs'].get('transcript_id')
                gene_id = record['attrs'].get('gene_id')
                
                self.transcript_info[transcript_id] = {
                    'gene_id': gene_id,
                    'chrom': record['chrom'],
                    'strand': record['strand'],
                    'start': record['start'],
                    'end': record['end']
                }
                transcript_count += 1
            
            elif feature_type == 'exon':
                transcript_id = record['attrs'].get('transcript_id')
                exon_number = record['attrs'].get('exon_number')
                
                exon_data = {
                    'start': record['start'],
                    'end': record['end'],
                    'exon_number': exon_number,
                    'chrom': record['chrom'],
                    'strand': record['strand']
                }
                self.transcripts[transcript_id].append(exon_data)
                exon_count += 1
        
        logger.info(f"Indexed: {gene_count} genes, {transcript_count} transcripts, {exon_count} exons")
        return gene_count, transcript_count, exon_count


class IntervalTreeBuilder:
    """Construct IntervalTrees per chromosome."""
    
    @staticmethod
    def build_trees(parser_instance):
        """
        Build IntervalTree for each chromosome from transcript exons.
        
        Args:
            parser_instance (GTFParser): Indexed GTF parser
            
        Returns:
            dict: chrom -> IntervalTree
        """
        trees = defaultdict(IntervalTree)
        
        for transcript_id, exons in parser_instance.transcripts.items():
            if transcript_id not in parser_instance.transcript_info:
                continue
            
            tx_info = parser_instance.transcript_info[transcript_id]
            gene_id = tx_info['gene_id']
            gene_name = parser_instance.genes[gene_id]['gene_name']
            chrom = tx_info['chrom']
            strand = tx_info['strand']
            
            # Sort exons by exon_number
            sorted_exons = sorted(exons, key=lambda e: (int(e['start']), int(e['end'])))
            
            # Add exons to tree
            for exon_idx, exon in enumerate(sorted_exons):
                # Interval is [start, end) - convert to 0-based half-open
                interval = Interval(
                    begin=exon['start'],
                    end=exon['end'],
                    data={
                        'gene_id': gene_id,
                        'gene_name': gene_name,
                        'transcript_id': transcript_id,
                        'type': 'exon',
                        'exon_number': exon['exon_number'],
                        'strand': strand
                    }
                )
                if interval.begin >= interval.end:
                    continue
                trees[chrom].add(interval)
            
            # Add introns between exons
            for intron_idx in range(len(sorted_exons) - 1):
                intron_start = sorted_exons[intron_idx]['end']
                intron_end = sorted_exons[intron_idx + 1]['start']
                
                interval = Interval(
                    begin=intron_start,
                    end=intron_end,
                    data={
                        'gene_id': gene_id,
                        'gene_name': gene_name,
                        'transcript_id': transcript_id,
                        'type': 'intron',
                        'intron_number': intron_idx + 1,
                        'strand': strand
                    }
                )
                trees[chrom].add(interval)
        
        logger.info(f"Built IntervalTrees for {len(trees)} chromosomes")
        return dict(trees)


def build_tx_trees(gtf_file, max_level=2, max_tsl=3):
    """
    Main processing pipeline.
    
    Args:
        gtf_file (str): Input GTF file
        bed_output (str): Output BED file
        tree_output (str): Optional pickle file for IntervalTrees
        max_level (int): Maximum annotation level
        max_tsl (int): Maximum transcript support level
    """
    logger.info(f"Processing GTF: {gtf_file}")
    logger.info(f"Filters: level <= {max_level}, TSL <= {max_tsl}")
    
    # Parse and index GTF
    parser = GTFParser(gtf_file, max_level=max_level, max_tsl=max_tsl)
    parser.index_features()
    
    # Build IntervalTrees
    logger.info("Building IntervalTrees...")
    trees = IntervalTreeBuilder.build_trees(parser)
    
    # Summary
    total_intervals = sum(len(tree) for tree in trees.values())
    logger.info(f"Total intervals in trees: {total_intervals}")
    
    return trees

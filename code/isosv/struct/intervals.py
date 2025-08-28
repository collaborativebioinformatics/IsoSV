"""
Interval tree implementation for candidate structuring.
Handles chromosome and SVTYPE-based organization of SV candidates.
"""

from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, field
from collections import defaultdict
import heapq


@dataclass
class Candidate:
    """Per-read candidate event structure as defined in project spec."""
    chrom: str
    pos: int
    end: int
    svtype: str  # INS, DEL, INV, BND, UNRES
    svlen: int
    read: str
    strand: str
    mapq: int
    cigar: str
    prim: bool
    mate: Optional[str] = None  # for BND
    flags: Dict[str, Any] = field(default_factory=dict)


@dataclass
class ClusterSV:
    """Clustered SV structure as defined in project spec."""
    chrom: str
    pos: int
    end: int
    svtype: str
    svlen: int
    support_sr: int  # split-read support
    support_pr: int  # paired-end support
    su: int  # total support
    cipos: Tuple[int, int]  # confidence interval for position
    ciend: Tuple[int, int]  # confidence interval for end
    quals: Dict[str, float]  # quality metrics
    annot: Dict[str, Any]  # annotations
    reads: List[str] = field(default_factory=list)  # supporting read IDs


class IntervalNode:
    """Node in the interval tree."""
    
    def __init__(self, start: int, end: int, candidate: Candidate):
        self.start = start
        self.end = end
        self.max_end = end
        self.candidates = [candidate]
        self.left = None
        self.right = None


class IntervalTree:
    """Interval tree for efficient range queries."""
    
    def __init__(self):
        self.root = None
    
    def insert(self, start: int, end: int, candidate: Candidate):
        """Insert a candidate into the interval tree."""
        if self.root is None:
            self.root = IntervalNode(start, end, candidate)
        else:
            self._insert_recursive(self.root, start, end, candidate)
    
    def _insert_recursive(self, node: IntervalNode, start: int, end: int, candidate: Candidate):
        """Recursively insert into the interval tree."""
        if start < node.start:
            if node.left is None:
                node.left = IntervalNode(start, end, candidate)
            else:
                self._insert_recursive(node.left, start, end, candidate)
        else:
            if node.right is None:
                node.right = IntervalNode(start, end, candidate)
            else:
                self._insert_recursive(node.right, start, end, candidate)
        
        # Update max_end
        node.max_end = max(node.max_end, end)
    
    def query(self, start: int, end: int) -> List[Candidate]:
        """Query for overlapping intervals."""
        results = []
        if self.root:
            self._query_recursive(self.root, start, end, results)
        return results
    
    def _query_recursive(self, node: IntervalNode, start: int, end: int, results: List[Candidate]):
        """Recursively query the interval tree."""
        if node is None:
            return
        
        # Check if current node overlaps with query
        if not (node.end < start or node.start > end):
            results.extend(node.candidates)
        
        # Check left subtree
        if node.left and node.left.max_end >= start:
            self._query_recursive(node.left, start, end, results)
        
        # Check right subtree
        if node.right and node.right.start <= end:
            self._query_recursive(node.right, start, end, results)


class CandidateStructurer:
    """Main class for structuring candidates using interval trees."""
    
    def __init__(self):
        # Separate interval trees for each chromosome and SVTYPE combination
        self.trees: Dict[Tuple[str, str], IntervalTree] = defaultdict(IntervalTree)
        self.candidates: List[Candidate] = []
    
    def add_candidate(self, candidate: Candidate):
        """Add a candidate to the appropriate interval tree."""
        self.candidates.append(candidate)
        key = (candidate.chrom, candidate.svtype)
        self.trees[key].insert(candidate.pos, candidate.end, candidate)
    
    def query_region(self, chrom: str, svtype: str, start: int, end: int) -> List[Candidate]:
        """Query candidates in a specific region."""
        key = (chrom, svtype)
        if key in self.trees:
            return self.trees[key].query(start, end)
        return []
    
    def get_all_candidates(self) -> List[Candidate]:
        """Get all candidates."""
        return self.candidates
    
    def get_candidates_by_type(self, svtype: str) -> List[Candidate]:
        """Get all candidates of a specific SV type."""
        return [c for c in self.candidates if c.svtype == svtype]
    
    def get_candidates_by_chrom(self, chrom: str) -> List[Candidate]:
        """Get all candidates on a specific chromosome."""
        return [c for c in self.candidates if c.chrom == chrom]



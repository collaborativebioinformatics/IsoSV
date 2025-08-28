"""
Interval tree implementation for candidate structuring using the intervaltree package.
Handles chromosome and SVTYPE-based organization of SV candidates.
"""

from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, field
from collections import defaultdict
from intervaltree import IntervalTree, Interval


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


class _ITreeWrapper:
    """Thin wrapper to provide a minimal compatible API over intervaltree.IntervalTree."""

    def __init__(self) -> None:
        self._tree: IntervalTree = IntervalTree()

    def insert(self, start: int, end: int, candidate: Candidate) -> None:
        # intervaltree uses half-open [begin, end); make end inclusive by adding +1
        self._tree.add(Interval(start, end + 1, candidate))

    def query(self, start: int, end: int) -> List[Candidate]:
        # Treat provided end as inclusive
        intervals = self._tree.overlap(start, end + 1)
        return [iv.data for iv in intervals]

    def at(self, point: int) -> List[Candidate]:
        intervals = self._tree.at(point)
        return [iv.data for iv in intervals]

    def clear(self) -> None:
        self._tree.clear()

    def __iter__(self):
        return iter(self._tree)

    def __len__(self) -> int:
        return len(self._tree)


class CandidateStructurer:
    """Main class for structuring candidates using interval trees."""
    
    def __init__(self):
        # Separate interval trees for each chromosome and SVTYPE combination
        self.trees: Dict[Tuple[str, str], _ITreeWrapper] = defaultdict(_ITreeWrapper)
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

    def query_point(self, chrom: str, svtype: str, point: int) -> List[Candidate]:
        """Query candidates overlapping a single genomic point."""
        key = (chrom, svtype)
        if key in self.trees:
            return self.trees[key].at(point)
        return []

    def query_overlap(self, chrom: str, svtype: str, start: int, end: int) -> List[Candidate]:
        """Alias for querying candidates overlapping [start, end] (inclusive)."""
        return self.query_region(chrom, svtype, start, end)
    
    def get_all_candidates(self) -> List[Candidate]:
        """Get all candidates."""
        return self.candidates
    
    def get_candidates_by_type(self, svtype: str) -> List[Candidate]:
        """Get all candidates of a specific SV type."""
        return [c for c in self.candidates if c.svtype == svtype]
    
    def get_candidates_by_chrom(self, chrom: str) -> List[Candidate]:
        """Get all candidates on a specific chromosome."""
        return [c for c in self.candidates if c.chrom == chrom]

    def get_tree_stats(self) -> Dict[str, Any]:
        """Return simple stats about each interval tree."""
        stats: Dict[str, Any] = {}
        for (chrom, svtype), tree in self.trees.items():
            stats[f"{chrom}_{svtype}"] = {
                "chrom": chrom,
                "svtype": svtype,
                "size": len(tree),
            }
        return stats

    def clear(self) -> None:
        """Clear all tracked candidates and interval trees."""
        self.candidates.clear()
        for tree in self.trees.values():
            tree.clear()
        self.trees.clear()

    def remove_candidate(self, candidate: Candidate) -> None:
        """Remove a specific candidate from storage and trees (linear scan)."""
        try:
            self.candidates.remove(candidate)
        except ValueError:
            pass
        key = (candidate.chrom, candidate.svtype)
        tree = self.trees.get(key)
        if not tree:
            return
        # Rebuild the tree without this candidate (intervaltree lacks remove by data)
        remaining: List[Interval] = [iv for iv in tree._tree if iv.data is not candidate]
        tree.clear()
        for iv in remaining:
            tree._tree.add(iv)

    def get_overlapping_candidates(self, candidate: Candidate, window: int = 0) -> List[Candidate]:
        """Return candidates overlapping candidate +/- window bp (excluding self)."""
        key = (candidate.chrom, candidate.svtype)
        tree = self.trees.get(key)
        if not tree:
            return []
        start = max(0, candidate.pos - window)
        end = candidate.end + window
        hits = tree.query(start, end)
        return [c for c in hits if c is not candidate]

    def merge_trees(self, other: "CandidateStructurer") -> None:
        """Merge candidates from another structurer into this one."""
        for cand in other.get_all_candidates():
            self.add_candidate(cand)

    def get_tree_visualization(self, chrom: str, svtype: str, max_intervals: int = 50) -> str:
        """Return a simple text representation of a tree for debugging."""
        key = (chrom, svtype)
        tree = self.trees.get(key)
        if not tree or len(tree) == 0:
            return f"No intervals for {chrom}_{svtype}"
        lines: List[str] = [f"Interval Tree for {chrom}_{svtype} ({len(tree)} intervals)"]
        lines.append("-" * 60)
        for idx, iv in enumerate(sorted(tree._tree, key=lambda x: (x.begin, x.end))[:max_intervals], start=1):
            cand: Candidate = iv.data
            lines.append(
                f"{idx:3d}. [{iv.begin:8d}, {iv.end - 1:8d}) {cand.svtype:4s} len={cand.svlen:6d} read={cand.read}"
            )
        if len(tree) > max_intervals:
            lines.append(f"... and {len(tree) - max_intervals} more intervals")
        return "\n".join(lines)



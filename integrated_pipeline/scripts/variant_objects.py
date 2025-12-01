from dataclasses import dataclass
from typing import Optional

@dataclass
class IndelObs:
    """Represents an InDel observation from a single read."""
    chrom: str
    pos: int
    svtype: str
    length: int
    mapq: int
    rname: str
    seq: Optional[str] = None

@dataclass
class ClipObs:
    """Represents a soft/hard clip observation from a single read."""
    chrom: str
    pos: int
    side: str
    cliplen: int
    mapq: int
    rname: str

@dataclass
class SplitObs:
    """Represents a split-read observation from an SA tag."""
    chrom1: str
    pos1: int
    chrom2: str
    pos2: int
    mapq: int
    rname: str
    same_chrom: bool
    note: str = ""
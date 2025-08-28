"""
Interval Tree Implementation

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import intervaltree as it




@dataclass
class Interval:
    def __init__(self, chrom, pos, end, svtype, read):
        self.chrom = chrom
        self.pos = pos
        self.end = end
        self.svtype = svtype
        self.read = read

class IntervalTree:
    def __init__(self):
        self.tree = it.IntervalTree()

    def insert(self, interval):
        self.tree.add(interval)
    
    def get_overlapping_intervals(self, interval):
        return self.tree.overlap(interval.pos, interval.end)
    
    def get_tree_visualization(self, chrom, svtype):
        return self.tree.visualize(chrom, svtype)
    
    def merge_overlapping_intervals(self, interval):
        return self.tree.merge_overlaps(interval)
    
    def get_interval_tree(self):
        return self.tree
    
    def get_interval_tree_visualization(self):
        return self.tree.visualize()




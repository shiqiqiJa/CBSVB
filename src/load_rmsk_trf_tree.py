import json
from intervaltree import Interval, IntervalTree

def load_interval_trees_from_json(filename):
    with open(filename, 'r') as json_file:
        tree_data = json.load(json_file)
    
    rmsk_trees = {}
    for chrom, intervals in tree_data.items():
        tree = IntervalTree()
        for start, end in intervals:
            tree.add(Interval(start, end))
        rmsk_trees[chrom] = tree
    
    return rmsk_trees

# my_tree = load_interval_trees_from_json('/home/qiqi/detect_translocation_p_copy/data/rmsk_trf_tree/chr1_10000_10000000.json')
import pandas as pd
from intervaltree import Interval, IntervalTree
import json

rmsk_df = pd.read_csv('./data/exp.rmsk.bed', sep="\t", header=None)
print(rmsk_df)

rmsk_intervals = {}

flag = 0
# 构建染色体的区间树
for _, row in rmsk_df.iterrows():
    chrom = row[0]
    start = row[1]
    end = row[2]

    if chrom not in rmsk_intervals:
        rmsk_intervals[chrom] = IntervalTree()
    rmsk_intervals[chrom].add(Interval(start, end))

tree_data = {}
for chrom, tree in rmsk_intervals.items():
    tree_data[chrom] = [(interval.begin, interval.end) for interval in tree]

with open('./data/exp.rmsk.json', 'w') as json_file:
    json.dump(tree_data, json_file)
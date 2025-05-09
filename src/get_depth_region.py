import pysam
import numpy as np
import multiprocessing
import os
import json

def get_filtered_depth_regions(bam_path, chrom, region_start, region_end):
    """
    获取 BAM 文件中深度符合正态分布的区域，舍弃深度最小的 10% 和最大的 10% 的区域。
    :param bam_file: BAM 文件路径
    :param threshold_percentile: 用于计算深度的上下限的百分位，默认为 10%（即最小的 10% 和最大的 10%）
    :return: 筛选后的深度区域，包含低深度区和高深度区
    """
    bam_file = pysam.AlignmentFile(bam_path)
    aligns = bam_file.fetch(chrom, region_start, region_end)
    
    # 获取参考基因组长度（假设参考基因组在第一个序列上）
    reference_length = region_end - region_start
    
    # 初始化深度信息
    depths = [0] * reference_length  # 存储每个位置的深度
    
    # 遍历 BAM 文件中的每个比对
    for read in aligns:
        if bam_file.getrname(read.reference_id) != chrom:
            continue
        start_pos = read.reference_start
        end_pos = read.reference_end
        
        # 统计总深度，只统计主比对，不统计次级比对
        if not read.is_secondary:
            for pos in range(start_pos, end_pos):
                if  region_start <= pos < region_end: # 这里改成在区间
                    depths[pos - region_start] += 1
        
    return chrom, region_start, region_end, depths

def get_results(bam_file):
    regions = []
    with multiprocessing.Pool(processes=16) as pool:
        results = []
        with open('./data/chr1.10M.bed') as fin:
            for region in fin:
                region_split = region.strip().split("\t")
                chrom, region_start, region_end = region_split[0], int(region_split[1]), int(region_split[2])
                result = pool.apply_async(get_filtered_depth_regions, (bam_file, chrom, region_start, region_end))
                results.append(result)

        pool.close()
        pool.join()

        for re in results:
            if re:
                result = re.get()
                chrom, region_start, region_end, depths = result[0], result[1], result[2], result[3]
                regions.append(((chrom, region_start, region_end), depths))
    
    return regions

bam_path = "/mnt/h/exp_bams/SRR28305180_chr1.bam"
regions = get_results(bam_path)
# 将深度数据转换为 NumPy 数组，便于后续统计
all_depths = []
for data in regions:
    chrom, region_start, region_end = data[0][0], data[0][1], data[0][2]
    if not os.path.exists('./data/depths'):
        os.mkdir('./data/depths')
    file = open('./data/depths/{}_{}_{}_ont.json'.format(chrom, region_start, region_end), 'w')
    json.dump(data[1], file)
    if data[0][1] % 10000 == 0 and data[0][2] % 10000 == 0:
        all_depths.extend(data[1])
depths_array = np.array(all_depths)

# 计算深度的均值和标准差
mean_depth = np.mean(depths_array)
std_depth = np.std(depths_array)

print(mean_depth, std_depth)

# 计算深度上下限，舍弃深度最小的 10%
lower_limit = mean_depth - 1.28 * std_depth  # 正态分布的下10%约为均值减去1.28倍标准差

# 找到低于下限和高于上限的区域
low_depth_regions = []

current_low_start = None

# 遍历每个位置，找到低于下限和高于上限的区域
for region in regions:
    chrom, start, end = region[0][0], int(region[0][1]), int(region[0][2])
    depths = region[1]
    for pos in range(len(depths)):
        depth = depths[pos]
        # 低于下限的区域
        if depth < lower_limit:
            if current_low_start is None:
                current_low_start = pos
        else:
            if current_low_start is not None:
                low_depth_regions.append((chrom, current_low_start+start, start+pos - 1))
                current_low_start = None
        

    # 处理最后一个区间
    if current_low_start is not None:
        low_depth_regions.append((chrom, current_low_start+start, end - 1))

with open('./data/low_depth_region_ont_chr1.bed', 'w') as f:
    for region in low_depth_regions:
        f.write("{}\t{}\t{}\n".format(region[0], region[1], region[2]))

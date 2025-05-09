from collections import Counter
import os

def filter_juncs_trfregion(include_juncs, end_info, rmsk_trf, data_type):
    if len(include_juncs) == 0:
        return False
    if data_type == 'tgs':
        judge_number = 1
    else:
        judge_number = 0
    repeat_region_num = 0
    for i in range(len(include_juncs)+1):
        if i==0:
            chrom = end_info[0].chrom
            pos1 = end_info[0].pos
            pos2 = include_juncs[0].include_bkps[0].pos
        elif i==len(include_juncs):
            chrom = end_info[1].chrom
            pos1 = include_juncs[len(include_juncs)-1].include_bkps[1].pos
            pos2 = end_info[1].pos
        elif i < len(include_juncs)-1:
            chrom = include_juncs[i].include_bkps[1].chrom
            pos1 = include_juncs[i].include_bkps[1].pos
            pos2 = include_juncs[i+1].include_bkps[0].pos
        else:
            continue
        start = min(pos1, pos2)
        end = max(pos1, pos2)
        start_flag = 0
        end_flag = 0
        for info in rmsk_trf[chrom]:
            if abs(info[0]-start)<30:
                start_flag = 1
            if abs(info[1]-end)<30:
                end_flag = 1
            if start_flag == 1 and end_flag == 1:
                repeat_region_num += 1
                break
            if info[1] - end > 300:
                break
    if repeat_region_num > judge_number:
        return False
    return True

def get_regions(options, mode):
    regions = dict()
    if mode == 'repeat_regions':
        file_path = options.repeat_region_path
    elif mode == 'rmsk_trf':
        file_path = options.rmsk_trf_path
    else:
        ValueError("Illegible mode")
    if file_path == None:
        return None
    with open(file_path) as region_file:
        for region in region_file:
            region_info = region.split('\n')[0].split('\t')
            if region_info[0] in options.all_possible_chrs:
                if region_info[0] not in regions.keys():
                    regions[region_info[0]] = []
                regions[region_info[0]].append((int(region_info[1]), int(region_info[2]), region_info[3]))
            else:
                ValueError("Illegible chromosome: {}".format(region_info[0]))
    region_file.close()
    return regions

def get_mismatch_ratio(cigarstring, reference_start, reference_end):
    op_counter = Counter(cigarstring)

    if 'X' not in op_counter:
        op_counter['X'] = 0
    
    if 'I' not in op_counter:
        op_counter['I'] = 0

    if 'D' not in op_counter:
        op_counter['D'] = 0
    
    mismatch_ratio = (op_counter['X']+op_counter['I']+op_counter['D']) / (reference_end - reference_start)
    
    return mismatch_ratio

def report_result(sample_name, juncSe_list, options):
    with open(os.path.join(options.out_path, "{}.txt".format(sample_name)), "w") as output_file:
        output_file.write('#id\tjunc_num\tsupport_reads\tdescription\n')
        for i in range(len(juncSe_list)):
            if i < len(juncSe_list)-1:
                output_file.write('Trans.{}\t{}\t{}\t{}\n'.format(i, juncSe_list[i].get_junc_num(), len(juncSe_list[i].support_reads), juncSe_list[i].to_string()))
            else:
                output_file.write('Trans.{}\t{}\t{}\t{}'.format(i, juncSe_list[i].get_junc_num(), len(juncSe_list[i].support_reads), juncSe_list[i].to_string()))
    output_file.close()
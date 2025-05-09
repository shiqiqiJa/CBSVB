import sys
import os
import datetime
import time
import shutil
import logging

import pysam

# from src.cluster_op import cluster_juncSe_list, cluster_juncSe_multijuncs, cluster_multi_juncs, judge_merge_inner
from src.cluster_op import cluster_juncSe_list
from src.reads_op_tgs import traverse_reads_tgs
from src.reads_op_ngs import traverse_reads_ngs
# from src.reads_op_ngs_unmated import get_juncSe_from_reads
from src.other_op import get_regions, report_result, get_reliable_juncSe, get_mismatch_ratio_threshold
from src.bridge_op import get_bridge_candidates
from src.combine_op import combine_results
from src.filter_multi_junc_op import get_no_and_multi_juncSe, get_no_and_multi_junc

import argparse
import multiprocessing
import pickle

def parse_arguments(arguments = sys.argv[1:]):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)


    required_params = parser.add_argument_group("Input/Output parameters")
    required_params.add_argument('-o', dest="out_path", type=os.path.abspath, required=True, help='Absolute path to output ')
    required_params.add_argument('-b', dest='bam_paths', type=str, nargs='+', required=True, help='Absolute path to bam files')
    required_params.add_argument('-seq_type', dest="seq_type", type=str, required=True, help='The type of sequence, ngs or hifi or ont')

    optional_params = parser.add_argument_group("Optional parameters")
    optional_params.add_argument('-thread_num', dest="thread_num", type=int, default=48, help='Thread numbers')
    optional_params.add_argument('-min_support', dest="min_support", type=int, default=3, help='Minimum support read number required for SV calling')
    optional_params.add_argument('-repeat_region_path', dest="repeat_region_path", type=os.path.abspath, default='./data/hg38.repeat_telomere_centromere.bed', help='Absolute path to repeat regions \
                                                                                                                                                                        (centromere, telemere) file')
    optional_params.add_argument('-accessed_region_path', dest="accessed_region_path", type=os.path.abspath, default='./data/hg38.access.10M.bed', help='Absolute path to accessed regions')
    optional_params.add_argument('-rmsk_trf_path', dest="rmsk_trf_path", type=os.path.abspath, default='./data/hg38.rmsk_trf.bed', help='Absolute path to trf file')
    optional_params.add_argument('-dist_threshold', dest="dist_threshold", type=int, default=1000, help='Minimum threshold of bkp required for SV calling')
    optional_params.add_argument('-non_local_thresh', dest="non_local_thresh", type=int, default=10000, help='Minimum threshold of non-local SV required for SV calling')
    optional_params.add_argument('-min_mapq', dest="min_mapq", type=int, default=20, help='Minimum mapq of read required for SV calling')
    optional_params.add_argument('-min_align_length', dest="min_align_length", type=int, default=100, help='Minimum length of align required for SV calling')
    optional_params.add_argument('-all_possible_chrs', dest="all_possible_chrs", type=str, nargs='+', default=["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", \
                                                                                                    "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", \
                                                                                                    "chr22", "chrX"], help='All possible chromosomes for SV calling')
    optional_params.add_argument('-mismatch_filter', dest="mismatch_filter",type=bool, default=True)
    optional_params.add_argument('-max_mismatch_ratio', dest="max_mismatch_ratio", type=float, default=0.04)
    optional_params.add_argument('-combine_result', dest="combine_result",type=bool, default=False)
    optional_params.add_argument('-combine_length_diff', dest = "combine_length_diff", type=int, default=50000)
    optional_params.add_argument('-detect_bridge', dest="detect_bridge",type=bool, default=False)
    optional_params.add_argument('-multi_junc', dest="report_multi_junc_results",type=bool, default=False)

    # options = parser.parse_args(arguments)
    options = parser.parse_args()

    return options

def transer(options, repeat_regions):
    if not os.path.exists(options.out_path):
        os.mkdir(options.out_path)

    start_time = datetime.datetime.now()
    if options.seq_type == 'ngs':
        unmated_bams = []

    for bam_path in options.bam_paths:
        print(bam_path)
        if not os.path.exists(bam_path):
            print(f"[Warning]: No bam {bam_path}")
            continue
        sample_name = os.path.basename(bam_path).split('.')[0]

        # try:
        #     bam_file = pysam.AlignmentFile(bam_path)
        # except Exception as e:
        #     print(f"[Error]: Could not open BAM file: {e}")
        #     continue

        regions = []
        final_juncSe_list = []
        chrom_juncSe_list = {}
        with open(options.accessed_region_path) as fin:
            for region in fin:
                region_split = region.strip().split("\t")
                chrom, region_start, region_end = region_split[0], int(region_split[1]), int(region_split[2])
                regions.append((chrom, region_start, region_end))
                if os.path.exists(options.out_path+'/juncSe_data'+'/juncSe_{}_{}_{}.pkl'.format(chrom, region_start, region_end)):
                    with open(options.out_path+'/juncSe_data'+'/juncSe_{}_{}_{}.pkl'.format(chrom, region_start, region_end), 'rb') as f:
                        juncSe_list = pickle.load(f)
                        print(chrom, region_start, region_end, len(juncSe_list))
                        if options.seq_type == 'ngs':
                            for juncSe in juncSe_list:
                                for junc in juncSe.include_juncs:
                                    chrom1 = min(junc.include_bkps[0].chrom, junc.include_bkps[1].chrom)
                                    chrom2 = max(junc.include_bkps[0].chrom, junc.include_bkps[1].chrom)
                                    if (chrom1, chrom2) not in chrom_juncSe_list.keys():
                                        chrom_juncSe_list[(chrom1, chrom2)] = []
                                    chrom_juncSe_list[(chrom1, chrom2)].append(juncSe)
                        else:
                            final_juncSe_list.extend(juncSe_list)
                if options.seq_type == 'ngs' and os.path.exists(options.out_path+'/unmated_data'+'/unmated_{}_{}_{}.pkl'.format(chrom, region_start, region_end)):
                    with open(options.out_path+'/unmated_data'+'/unmated_{}_{}_{}.pkl'.format(chrom, region_start, region_end), 'rb') as f_:
                        unmated_bam = pickle.load(f_)
                        unmated_bams.append(unmated_bam)

        if len(final_juncSe_list) == 0 and len(chrom_juncSe_list.keys()) == 0:
            with multiprocessing.Pool(processes=options.thread_num) as pool:
                results = []
                for region in regions:
                    if options.seq_type == 'hifi' or options.seq_type == 'ont':
                        result = pool.apply_async(traverse_reads_tgs, (bam_path, region[0], region[1], region[2], repeat_regions, options))
                    elif options.seq_type == 'ngs':
                        result = pool.apply_async(traverse_reads_ngs, (bam_path, region[0], region[1], region[2], repeat_regions, "mated", options))
                    else:
                        raise ValueError("Wrong seq type")
                    results.append(result)

                pool.close()
                pool.join()

                # 从队列获取结果
                for re in results:
                    if re:
                        result = re.get()
                        if options.seq_type == 'hifi' or options.seq_type == 'ont':
                            juncSe_list = result
                        elif options.seq_type == 'ngs':
                            juncSe_list, unmated_bam = result[0], result[1]
                            print(unmated_bam)
                            unmated_bams.append(unmated_bam)
                        if len(juncSe_list) != 0:
                            final_juncSe_list.extend(juncSe_list)

        if options.seq_type == 'ngs':
            merged_unmated_prefix = bam_path.split("/")[-1].replace(".bam", ".unmated.bam")
            merged_unmated_sort_prefix = bam_path.split("/")[-1].replace(".bam", ".sorted.unmated.bam")
            merged_unmated_bam = os.path.join(options.out_path, merged_unmated_prefix)
            merged_unmated_sort_bam = os.path.join(options.out_path, merged_unmated_sort_prefix)
            if not os.path.exists(merged_unmated_sort_bam):
                # 在进行merge之前需要先删除之前的“.bam”文件
                cmd_str = "pysam.merge(\"-@\", \"{}\", \"{}\", ".format(options.thread_num, merged_unmated_bam) 
                for i in range(len(unmated_bams)-1):
                    if os.path.exists(unmated_bams[i]):
                        cmd_str += "\"{}\", ".format(unmated_bams[i])
                if os.path.exists(unmated_bams[len(unmated_bams)-1]):
                    cmd_str += "\"{}\" ".format(unmated_bams[len(unmated_bams)-1])
                cmd_str += ")"
                exec(cmd_str)

                pysam.sort(merged_unmated_bam, "-o", merged_unmated_sort_bam)
                pysam.index("-@", "{}".format(options.thread_num), merged_unmated_sort_bam)
            else:
                logging.info('got unmated bam')

            print(merged_unmated_sort_bam)

            # # then traverse unmated bam
            if not os.path.exists(options.out_path+'/juncSe_data'+'/juncSe_None_None_None.pkl'):
                origin_juncSe_list,  _ = traverse_reads_ngs(merged_unmated_sort_bam, None, None, None, repeat_regions, "unmated", options)
            else:
                with open(options.out_path+'/juncSe_data'+'/juncSe_None_None_None.pkl', 'rb') as f:
                    origin_juncSe_list = pickle.load(f)
            
            if len(origin_juncSe_list) != 0:
                for juncSe in origin_juncSe_list:
                    for junc in juncSe.include_juncs:
                        chrom1 = min(junc.include_bkps[0].chrom, junc.include_bkps[1].chrom)
                        chrom2 = max(junc.include_bkps[0].chrom, junc.include_bkps[1].chrom)
                        if (chrom1, chrom2) not in chrom_juncSe_list.keys():
                            chrom_juncSe_list[(chrom1, chrom2)] = []
                        chrom_juncSe_list[(chrom1, chrom2)].append(juncSe)
                print('got all results')

        if options.seq_type == 'ngs':
            for chrom_pair in chrom_juncSe_list.keys():
                chrom_juncSe_list[chrom_pair] = sorted(chrom_juncSe_list[chrom_pair],key=lambda x:x.get_junc_num(),reverse=True)
            with multiprocessing.Pool(4) as pool:
                results = []
                for chrom_pair in chrom_juncSe_list.keys():
                    print(chrom_pair, len(chrom_juncSe_list[chrom_pair]))
                    result = pool.apply_async(cluster_juncSe_list, (chrom_juncSe_list[chrom_pair], 0))
                    results.append(result)

                pool.close()
                pool.join()

                # 从队列获取结果
                for re in results:
                    if re:
                        result = re.get()
                        final_juncSe_list.extend(result)

        srt_all_juncSe_list = sorted(final_juncSe_list,key=lambda x:x.get_junc_num(),reverse=True)
        clustered_juncSe_list = cluster_juncSe_list(srt_all_juncSe_list, 0)
        
        filter_juncSe = []
        for juncSe in clustered_juncSe_list:
            if juncSe.get_support_num() >= options.min_support:
                filter_juncSe.append(juncSe)

        srt_juncSe_list = sorted(filter_juncSe, key=lambda x:x.get_support_num(), reverse=True)

        if options.seq_type == 'hifi' or options.seq_type == 'ont':
            sample_out_file = os.path.join(options.out_path, "{}_tgs.txt".format(sample_name))
        elif options.seq_type == 'ngs':
            sample_out_file = os.path.join(options.out_path, "{}_ngs.txt".format(sample_name))
        pickle_out_file = os.path.join(options.out_path, "{}_juncSe_data.pickle".format(sample_name))
        with open(sample_out_file, "w") as fout:
            for juncSe in srt_juncSe_list:
                fout.write(juncSe.to_string() + "\n")
        with open(pickle_out_file, "wb") as f:
            pickle.dump(srt_juncSe_list, f)
        report_result(sample_name, srt_juncSe_list, options)

    end_time = datetime.datetime.now()
    print(f"Cost {int((end_time - start_time).total_seconds())}s")
    return srt_juncSe_list



if __name__ == '__main__':

    start_time = datetime.datetime.now()

    options = parse_arguments()

    work_dir = options.out_path
    if not os.path.exists(work_dir):
        print(f'Create the output directory {work_dir}')
        os.mkdir(work_dir)
    
    log_format = logging.Formatter("%(asctime)s [%(levelname)-7.7s]  %(message)s")
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)

    fileHandler = logging.FileHandler("transer.log", mode="w")
    fileHandler.setFormatter(log_format)

    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(log_format)
    root_logger.addHandler(fileHandler)

    logging.info('**************************** transer ****************************')
    logging.info("CMD: {0}".format(" ".join(sys.argv)))
    logging.info("WORKDIR DIR: {0}".format(os.path.abspath(work_dir)))

    sample_name = os.path.basename(options.bam_paths[0]).split('.')[0]
    repeat_regions = get_regions(options, "repeat_regions")

    if not os.path.exists(os.path.join(options.out_path, "{}_juncSe_data.pickle".format(sample_name))):
        srt_juncSe_list = transer(options, repeat_regions)
        # mismatch_ratio_threshold = get_mismatch_ratio_threshold(os.path.join(options.out_path, 'reads_mismatch_ratio_data'))
        # output_file_path = work_dir + '/SV.txt'
        # with open(output_file_path, 'w') as output_file:
        #     for juncSe in srt_juncSe_list:
        #         output_file.write("{}\n{}\n".format(juncSe.to_string(), juncSe.get_support_num()))
        # reliable_srt_juncSe_list = get_reliable_juncSe(srt_juncSe_list, mismatch_ratio_threshold)
        # re_output_file_path = work_dir + '/reliable_SV.txt'
        # with open(re_output_file_path, 'w') as re_output_file:
        #     for juncSe in reliable_srt_juncSe_list:
        #         re_output_file.write("{}\n{}\n".format(juncSe.to_string(), juncSe.get_support_num()))

    end_time1 = datetime.datetime.now()
    logging.info("WORK TIME NO COMBINE: {0}".format(end_time1-start_time))
    
    if options.detect_bridge or options.report_multi_junc_results or options.combine_result:
        pickle_file = os.path.join(options.out_path, "{}_juncSe_data.pickle".format(sample_name))
        with open(pickle_file, 'rb') as f:
            srt_juncSe_list = pickle.load(f)

        output_file_path = work_dir + '/SV.txt'
        with open(output_file_path, 'w') as output_file:
            for juncSe in srt_juncSe_list:
                output_file.write("{}\n{}\n".format(juncSe.to_string(), juncSe.get_support_num()))
        
        mismatch_ratio_threshold = get_mismatch_ratio_threshold(os.path.join(options.out_path, 'reads_mismatch_ratio_data'))
        reliable_srt_juncSe_list = get_reliable_juncSe(srt_juncSe_list, mismatch_ratio_threshold, options.seq_type)
        print(len(reliable_srt_juncSe_list))
        # filter_juncSe_list, multi_juncSe_list, multi_SV_juncSe = get_no_and_multi_juncSe(srt_juncSe_list, options.seq_type)
        filter_junc_list, multi_junc_list = get_no_and_multi_junc(reliable_srt_juncSe_list, options.seq_type)
        
        re_output_file_path = work_dir + '/reliable_SV.txt'
        with open(re_output_file_path, 'w') as re_output_file:
            for junc in filter_junc_list:
                re_output_file.write("{}\n{}\n".format(junc.to_string(), junc.updated_support_reads))

    if options.report_multi_junc_results:
        # print(work_dir+'/no_multi_juncs_junSe.txt')
        # no_multi_junc_file = open(work_dir+'/no_multi_juncs_junSe.txt','w')
        # no_multi_junc_file.write("#id\tjunc_num\tsupport_reads\tdescription\n")
        # for i in range(len(filter_juncSe_list)):
        #     juncSe = filter_juncSe_list[i]
        #     no_multi_junc_file.write("{}\t{}\t{}\t{}\n".format(i, juncSe.get_junc_num(), len(juncSe.support_reads), juncSe.to_string()))
        # if len(multi_SV_juncSe) > 0:
        #     no_multi_junc_file.write("multi_SV_one_bkp:")
        #     num = 0
        #     for multi_SV_juncSe in multi_SV_juncSe:
        #         support_num = len(multi_SV_juncSe.support_reads[0])+len(multi_SV_juncSe.support_reads[1])
        #         no_multi_junc_file.write("{}\t{}\t{}\n".format(num, support_num, multi_SV_juncSe.to_string()))
        #         num += 1
        # if len(multi_juncSe_list) > 0:
        multi_junc_file = open(work_dir+'/multi_juncs.txt','w')
        num = 0
        multi_junc_file.write("#id\tdescription\n")
        for multi_junc in multi_junc_list:
            multi_junc_file.write("{}\t{}\n".format(num, multi_junc.to_string()))
            num += 1

    # if options.combine_result:
        # combined_juncSe_list = combine_results(filter_juncSe_list, options)
        # with open(work_dir+'/combined_results.txt','w') as my_trans:
        #     for juncSe in combined_juncSe_list:
        #         my_trans.write(juncSe.to_string()+'\n')
        # end_time2 = datetime.datetime.now()
        # logging.info("WORK TIME AFTER COMBINE: {0}".format(end_time2-start_time))
    
    # if options.detect_bridge:
    #     bridges = []
    #     bridge_candidates = get_bridge_candidates(combined_juncSe_list)
    #     for bridge_candidate in bridge_candidates:
    #         bridge_candidate.judge_bridge_type()
    #         if  bridge_candidate.type != None:
    #             flag = 0
    #             for now_bridge in bridges:
    #                 if now_bridge.juncSe1_index == bridge_candidate.juncSe2_index and now_bridge.junc1_index == bridge_candidate.junc2_index and \
    #                 now_bridge.bkp1_index == bridge_candidate.bkp2_index:
    #                     flag = 1
    #                     break
    #             if flag==0:
    #                 bridges.append(bridge_candidate)
    #     bridge_file = open(work_dir+'/bridges.txt','w')
    #     for bridge in bridges:
    #         bridge_file.write(bridge.to_string()+'\n')

    end_time = datetime.datetime.now()
    logging.info("WORK TIME: {0}".format(end_time-start_time))
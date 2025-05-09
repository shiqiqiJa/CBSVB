# from classes import Junction_Series, Multipul_Junc
from src.classes import Junction_Series

# def judge_merge_inner(multi_junc):
#     if multi_junc.judge_minus():
#         merge, _, _, _, minus_multi_junc = compare_two_juncSe(multi_junc.support_junc_series[0], multi_junc.support_junc_series[1], 1)
#         if merge == True and minus_multi_junc != None:
#             multi_junc.same_bkp_info = minus_multi_junc.same_bkp_info
#             multi_junc.overlap = minus_multi_junc.overlap
#         else:
#             return False
    
#     index_1 = multi_junc.same_bkp_info[0][0]
#     start_1 = index_1 - multi_junc.overlap
#     index_2 = multi_junc.same_bkp_info[1][0]
#     start_2 = index_2 - multi_junc.overlap
#     if start_1 != 0 and start_2 != 0:
#         return False
#     index1 = start_1
#     index2 = start_2
#     while index1<index_1 and index2<index_2:
#         merge_flag, merge_num = compare_two_junc(multi_junc.support_junc_series[0].include_juncs[index1], multi_junc.support_junc_series[1].include_juncs[index2], 0)
#         if merge_flag is False:
#             return False
#         index1 += 1
#         index2 += 1
    
#     return True

# def judge_merge_multijuncs(base_multijunc, target_multijunc, dist_thresh=1000):
#     merge_juncs_info = None
#     merge_info = []
#     merge_before_juncs, max_start_1, max_start_2, max_overlap, _ = compare_two_juncSe(base_multijunc.merge_juncs, target_multijunc.merge_juncs, 0)
#     if merge_before_juncs is False:
#         return False, merge_juncs_info, merge_info
#     else:
#         merge_juncs_info = (max_start_1, max_start_2, max_overlap)

#     merge_bkp = compare_two_bkp(base_multijunc.get_same_bkp(), target_multijunc.get_same_bkp(), dist_thresh=dist_thresh)
#     if merge_bkp is False:
#         return False, merge_juncs_info, merge_info
#     else:
#         for i in range(len(target_multijunc.diff_juncs)):
#             merge_flag = 0
#             for j in range(len(base_multijunc.diff_juncs)):
#                 if compare_two_bkp(base_multijunc.diff_juncs[j].include_bkps[1], target_multijunc.diff_juncs[i].include_bkps[1]):
#                     merge_info.append((i, j))
#                     merge_flag = 1
#                     break
#             if merge_flag == 0:
#                 merge_info.append((i, -1))

#     return True, merge_juncs_info, merge_info

# def compare_juncSe_multijuncs(juncSe, multijuncs, dist_thresh=1000):
#     merge_juncs_info = None
#     for i in range(len(juncSe.include_juncs)):
#         merge = compare_two_bkp(juncSe.include_juncs[i].include_bkps[0], multijuncs.get_same_bkp(), dist_thresh=dist_thresh)
#         merge_index = -1
#         for j in range(len(multijuncs.diff_juncs)):
#             if compare_two_bkp(juncSe.include_juncs[i].include_bkps[1], multijuncs.diff_juncs[j].include_bkps[1]):
#                 merge_index = j
#                 break
#         if merge:
#             new_end = juncSe.include_juncs[i].include_bkps[0]
#             new_juncse = Junction_Series(juncSe.include_juncs[:i], juncSe.support_reads, juncSe.whole_segments_type[:i-1], 
#                                             juncSe.segments_direct[:i+1], [juncSe.end[0], new_end])
#             merge_juncs, max_start_1, max_start_2, max_overlap, _ = compare_two_juncSe(multijuncs.merge_juncs, new_juncse, 0)
#             if merge_juncs:
#                 merge_juncs_info = (max_start_1, max_start_2, max_overlap)
#                 return True, merge_juncs_info, (i, merge_index)
    
#     return False, merge_juncs_info, None


def compare_two_bkp(base_bkp, target_bkps, dist_thresh=1000):
    """
    compare two bkps to check whether they are the same
    :param base_bkp:
    :param target_bkps:
    :return:
    """

    # different chrom, return False
    if base_bkp.chrom != target_bkps.chrom:
        return False

    # different direct, return false
    if base_bkp.direct != target_bkps.direct:
        return False

    # large distance, return false
    if abs(base_bkp.pos - target_bkps.pos) > dist_thresh:
        return False

    return True


def compare_two_junc(base_junc, target_junc, dist_thresh=1000):
    """
    compare two juncs to check whether they are the same
    :param base_junc:
    :param target_junc:
    :param dist_thresh:
    :return:
    """

    # # different junc type, return false
    if base_junc.junc_type != target_junc.junc_type:
        return False, -1

    # # different junc sub type, return false
    if base_junc.junc_sub_type != target_junc.junc_sub_type:
        return False, -1

    base_junc_bkps = base_junc.include_bkps
    target_junc_bkps = target_junc.include_bkps

    # # different bkp number, return false
    if len(base_junc_bkps) != len(target_junc_bkps):
        return False, -1

    # # different bkps, return false
    bkp_same = []
    for i in range(len(base_junc_bkps)):
        if compare_two_bkp(base_junc_bkps[i], target_junc_bkps[i], dist_thresh=dist_thresh) is True:
            bkp_same.append(i)
    
    if len(bkp_same) == 1:
        return  False, bkp_same[0]
    
    if len(bkp_same) == 0:
        return False, -1

    # # final return
    return True, -1

def compare_two_juncSe(juncSe1, juncSe2, multi_flag, min_overlap=2, dist_thresh=1000):
    multi_junc = None
    """
    compare two juncSe by finding overlaped juncs between two juncSe, we take reciprocal into account
    :param juncSe1:
    :param juncSe2:
    :param min_overlap:
    :param reciprocal:
    :return:
    """
    flag=0
    max_overlap = -1
    max_start_1 = -1
    max_start_2 = -1

    if juncSe1 == None or juncSe2 == None:
        return True, -1,-1,-1, multi_junc

    if len(juncSe2.include_juncs)>len(juncSe1.include_juncs):
        juncse1=juncSe2
        juncse2=juncSe1
        flag = 1
    else:
        juncse1=juncSe1
        juncse2=juncSe2
    for i in range(len(juncse1.include_juncs)+len(juncse2.include_juncs)-1):
        index_on_1 = max(0,len(juncse1.include_juncs)-1-i)
        index_on_2 = max(0,i-len(juncse1.include_juncs)+1)
        overlap = 0
        start_on_1 = -1
        start_on_2 = -1
        while index_on_1<len(juncse1.include_juncs) and index_on_2<len(juncse2.include_juncs):
            merge_junc, merge_num = compare_two_junc(juncse1.include_juncs[index_on_1], juncse2.include_juncs[index_on_2], dist_thresh=dist_thresh)
            if merge_junc is False:
            #     if multi_flag == 1 and merge_num != -1:
            #         multi_junc = Multipul_Junc([(index_on_1, merge_num), (index_on_2, merge_num)], overlap, [juncse1, juncse2], [juncse1.support_reads, juncse2.support_reads])
                break
            else:
                overlap += 1
                if overlap == 1:
                    start_on_1 = index_on_1
                    start_on_2 = index_on_2
                index_on_1 += 1
                index_on_2 += 1
        if overlap>= min(min_overlap,juncSe1.get_junc_num(), juncSe2.get_junc_num()) and overlap>max_overlap and (index_on_1==len(juncse1.include_juncs) or index_on_2==len(juncse2.include_juncs)):
            max_overlap = overlap
            max_start_1 = start_on_1
            max_start_2 = start_on_2
    if max_overlap == -1:
        return False,-1,-1,-1, multi_junc
    if flag==1:
        start = max_start_1
        max_start_1 = max_start_2
        max_start_2 = start
    return True, max_start_1, max_start_2, max_overlap, multi_junc


def cluster_juncSe_list(origin_juncSe_list, multi_flag=0, min_overlap=2, dist_thresh=1000):
    multi_juns_list = []
    """
    cluster trouble reads to junc series
    :param reads_list:
    :param consider_minus_strand:
    :param min_overlap:
    :param dist_thresh:
    :return:
    """
    juncSe_list = []
    #####################
    for target_juncSe in origin_juncSe_list:
        merge_flag = 0
        # traverse each base juncSe to find whether can merge
        for base_juncSe in juncSe_list:
            merge, overlap_start_on_base, overlap_start_on_target, overlap_num, multi_juncs = compare_two_juncSe(base_juncSe, target_juncSe, multi_flag, min_overlap=min_overlap, dist_thresh=dist_thresh)
            if merge is True:
                base_juncSe.merge(target_juncSe, overlap_start_on_base, overlap_start_on_target, overlap_num)
                merge_flag = 1
                break
            # if multi_flag == 1 and multi_juncs != None:
            #     multi_juns_list.append(multi_juncs)
            #     juncSe_list.remove(base_juncSe)
            #     merge_flag = 1
            #     break

            # # after compare, if not merged, then we try target read's minus
            if merge_flag == 0:
                target_juncSe_minus = target_juncSe.get_include_juncs_minus()
                merge_minus, overlap_start_on_base_minus, overlap_start_on_target_minus, overlap_num_minus, multi_juncs_minus = compare_two_juncSe(base_juncSe, target_juncSe_minus, multi_flag, min_overlap=min_overlap, dist_thresh=dist_thresh)
                if merge_minus is True:
                    merge_flag = 1
                    base_juncSe.merge(target_juncSe_minus, overlap_start_on_base_minus, overlap_start_on_target_minus, overlap_num_minus)
                    break
                # if multi_flag == 1 and multi_juncs_minus != None:
                #     multi_juns_list.append(multi_juncs_minus)
                #     juncSe_list.remove(base_juncSe)
                #     merge_flag = 1
                #     break

        # after compare minus, if still cannot merge, then we append it
        if merge_flag == 0:
            juncSe_list.append(target_juncSe)
            continue

    return juncSe_list, multi_juns_list

# def merge_separate_juncSe_lists(separate_juncSe_lists, options):
#     """
#     merge separate juncSe lists to one, and sort by support
#     :param separate_juncSe_lists:
#     :return:
#     """

#     all_juncSe_list = []
#     for region in separate_juncSe_lists:
#         all_juncSe_list.extend(separate_juncSe_lists[region])
#     srt_all_juncSe_list = sorted(all_juncSe_list,key=lambda x:x.get_junc_num(),reverse=True)#这一步是因为聚类的时候要判断断点数量，如果把数量一样的放一起更方便吗？
    
#     clustered_juncSe_list, single_same_list = cluster_juncSe_list(srt_all_juncSe_list, options) #这里后面还得改，现在‘cluster_juncSe_list’不一样了
#     srt_juncSe_list = sorted(clustered_juncSe_list, key=lambda x:x.get_support_num(), reverse=True)

#     filterd_juncSe_list = []

#     for juncSe in srt_juncSe_list:
#         if juncSe.get_support_num() >= options.min_support:
#             filterd_juncSe_list.append(juncSe)
            
#     return filterd_juncSe_list

# def cluster_multi_juncs(multi_juncs_list, dist_thresh=1000):
#     clusted_multi_juncs = []
#     for target_multijunc in multi_juncs_list:
#         merge_flag = 0
#         for base_multijunc in clusted_multi_juncs:
#             # print("target_multijunc:{}".format(target_multijunc.to_string()))
#             # print("base_multijunc:{}".format(base_multijunc.to_string()))
#             merge_multijuncs, merge_juncs_info, merge_info = judge_merge_multijuncs(base_multijunc, target_multijunc, dist_thresh=dist_thresh)
#             # print("merge_multijuncs:{}, merge_juncs_info:{}, merge_info:{}".format(merge_multijuncs, merge_juncs_info, merge_info))
#             if merge_multijuncs:
#                 base_multijunc.merge(target_multijunc, merge_juncs_info, merge_info, 0)
#                 merge_flag = 1
#                 break
        
#         if merge_flag == 0:
#             clusted_multi_juncs.append(target_multijunc)
        
#     return clusted_multi_juncs

# def cluster_juncSe_multijuncs(juncSe_list, multijuncs_list, dist_thresh=1000):
#     for juncSe in juncSe_list:
#         merge_flag = 0
#         for multijuncs in multijuncs_list:
#             merge, merg_juncs_info, merge_info = compare_juncSe_multijuncs(juncSe, multijuncs, dist_thresh=dist_thresh)
#             if merge:
#                 juncSe_list.remove(juncSe)
#                 multijuncs.merge(juncSe, merg_juncs_info, merge_info, 1)
#                 merge_flag = 1
#                 break
#             if merge_flag == 0:
#                 juncSe_minus = juncSe.get_include_juncs_minus()
#                 merge_minus, merg_juncs_minus_info, merge_minus_info = compare_juncSe_multijuncs(juncSe_minus, multijuncs, dist_thresh=dist_thresh)
#                 if merge_minus:
#                     juncSe_list.remove(juncSe)
#                     multijuncs.merge(juncSe_minus, merg_juncs_minus_info, merge_minus_info, 1)
#                     break
    
#     return juncSe_list, multijuncs_list
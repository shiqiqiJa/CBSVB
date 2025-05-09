from src.classes import Bridge_Candidate

def judge_ajacent(bkp1, bkp2):
    if bkp1.chrom != bkp2.chrom:
        return False
    if abs(bkp1.pos-bkp2.pos)>1000:
        return False
    if bkp1.direct == bkp2.direct:
        return False
    return True

def get_bridge_candidates(juncSe_list):
    bridge_candidates = []
            
    filter_juncSe_list = juncSe_list

    for i in range(len(filter_juncSe_list)):
        for j in range(len(filter_juncSe_list[i].include_juncs)):
            if filter_juncSe_list[i].include_juncs[j].junc_sub_type == 'INS':
                continue
            for m in range(2):
                flag = 0
                bkp1 = filter_juncSe_list[i].include_juncs[j].include_bkps[1-m]
                if bkp1.chrom == filter_juncSe_list[i].include_juncs[j].include_bkps[m].chrom and abs(bkp1.pos-filter_juncSe_list[i].include_juncs[j].include_bkps[m].pos)<50000:
                    continue
                for k in range(i, len(filter_juncSe_list)):
                    if k==i:
                        start = j+1
                    else:
                        start = 0
                    for l in range(start, len(filter_juncSe_list[k].include_juncs)):
                        if filter_juncSe_list[k].include_juncs[l].junc_sub_type == 'INS':
                            continue
                        if k==i and abs(j-l)==1 and abs(filter_juncSe_list[i].include_juncs[min(j,l)].include_bkps[1].pos-filter_juncSe_list[i].include_juncs[max(j,l)].include_bkps[0].pos)<3000:
                            continue
                        for n in range(2):
                            bkp2 = filter_juncSe_list[k].include_juncs[l].include_bkps[1-n]
                            if judge_ajacent(filter_juncSe_list[i].include_juncs[j].include_bkps[m], filter_juncSe_list[k].include_juncs[l].include_bkps[n]) and \
                                ((bkp1.chrom != bkp2.chrom) or (bkp1.chrom==bkp2.chrom and abs(bkp1.pos-bkp2.pos)>50000)):
                                flag = 1
                                bridge_candidates.append(Bridge_Candidate(i,j,m,k,l,n,filter_juncSe_list))
                                break
                        if flag == 1:
                            break
                    if flag == 1:
                        break
    
    return bridge_candidates


# bridges = []
# bridge_candidates = get_bridge_candidates(srt_juncSe_list_tgs)
# for bridge_candidate in bridge_candidates:
#     bridge_candidate.judge_bridge_type()
#     if  bridge_candidate.type != None:
#         flag = 0
#         for now_bridge in bridges:
#             if now_bridge.juncSe1_index == bridge_candidate.juncSe2_index and now_bridge.junc1_index == bridge_candidate.junc2_index and \
#             now_bridge.bkp1_index == bridge_candidate.bkp2_index:
#                 flag = 1
#                 break
#         if flag==0:
#             bridges.append(bridge_candidate)

# inter_DEL = 0
# inter_DUP = 0
# intra_DEL = 0
# intra_DUP = 0

# # bridge_file = open('./result_out/A_bridge.txt','w')

# for bridge in bridges:
#     print(bridge.to_string())
#     # bridge_file.write(bridge.to_string()+'\n')
#     bridge_type = bridge.to_string().split('_')[1]
#     bridge_info = bridge.to_string().split('-')
#     bkp1_chr = bridge_info[0].split(':')[0]
#     bkp2_chr = bridge_info[-2].split(':')[0]
#     if bkp1_chr == bkp2_chr:
#         if bridge_type == 'DEL-bridge':
#             intra_DEL += 1
#         else:
#             intra_DUP += 1
#     else:
#         if bridge_type == 'DEL-bridge':
#             inter_DEL += 1
#         else:
#             inter_DUP += 1

# print("inter_DEL:{}, inter_DUP:{}, intra_DEL:{}, intra_DUP:{}".format(inter_DEL, inter_DUP, intra_DEL, intra_DUP))
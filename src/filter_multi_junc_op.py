from src.classes import Bkp, Single_same_Multi_junc
from src.cluster_op import compare_two_bkp

from Bio import pairwise2
from Bio.Seq import Seq

def get_seq(junc, bkp_index):
    if junc.junc_sub_type == 'INS':
        for aligns in junc.support_reads:
            for align in aligns:
                if align.chr == None:
                    return align.seq
    else:
        for aligns in junc.support_reads:
            for align in aligns:
                if align.chr == junc.include_bkps[bkp_index].chrom:
                    if junc.include_bkps[bkp_index].direct == 'L' and abs(junc.include_bkps[bkp_index].pos-align.ref_start)<1000:
                        return align.seq
                    if junc.include_bkps[bkp_index].direct == 'S' and abs(junc.include_bkps[bkp_index].pos-align.ref_end)<1000:
                        return align.seq
    return None
    

def compare_juncSe_for_multi_junc(juncSe1, juncSe2):
    for i in range(len(juncSe1.include_juncs)):
        for l in range(2):
            for j in range(len(juncSe2.include_juncs)):
                for m in range(2):
                    if compare_two_bkp(juncSe1.include_juncs[i].include_bkps[l], juncSe2.include_juncs[j].include_bkps[m]):
                        seq1 = Seq(get_seq(juncSe1.include_juncs[i], 1-l))
                        seq2 = Seq(get_seq(juncSe2.include_juncs[j], 1-m))
                        alignments = pairwise2.align.globalxx(seq1, seq2)
                        srt_alignments = sorted(alignments, key=lambda x:x.score, reverse=True)
                        ratio = srt_alignments[0].score/min(len(seq1), len(seq2))
                        seqs = [srt_alignments[0].seqA, srt_alignments[0].seqB]
                        if ratio > 0.8:
                            return 1, Bkp(juncSe1.include_juncs[i].include_bkps[l].chrom, int((juncSe1.include_juncs[i].include_bkps[l].pos+juncSe2.include_juncs[j].include_bkps[m].pos)/2), juncSe1.include_juncs[i].include_bkps[l].direct, juncSe1.include_juncs[i].include_bkps[l].type), \
                                juncSe1.include_juncs[i].include_bkps[1-l], juncSe2.include_juncs[j].include_bkps[1-m], seqs

    return 0, None, None, None, []

def get_no_and_multi_juncSe(juncSe_list):
    no_multi_juncSe = [] # for bridge search
    multi_juncSe = []

    for juncSe in juncSe_list:
        get_multi_flag = 0
        for target_juncSe in no_multi_juncSe:
            same_flag, same_bkp_info, juncSe1_info, juncSe2_info, seqs = compare_juncSe_for_multi_junc(target_juncSe, juncSe)
            if same_flag == 1:
                multi_juncSe.append(Single_same_Multi_junc(same_bkp_info, juncSe1_info, juncSe2_info, [target_juncSe, juncSe], [target_juncSe.support_reads, juncSe.support_reads], seqs))
                no_multi_juncSe.remove(target_juncSe)
                get_multi_flag = 1
                break
        if get_multi_flag == 1:
            continue
        no_multi_juncSe.append(juncSe)
    
    return no_multi_juncSe, multi_juncSe
# class Options:
#     def __init__(self, bam_paths, seq_type, thread_num, min_support, rmsk_trf_path, merge_dis_threshold, non_local_thresh, min_align_length, min_mapq):

#         self.bam_paths = bam_paths

#         self.out_path = './result_out'

#         self.seq_type = seq_type

#         self.thread_num = thread_num

#         self.min_support = min_support

#         self.rmsk_trf_path = rmsk_trf_path

#         self.merge_dis_threshold = merge_dis_threshold

#         self.non_local_thresh = non_local_thresh

#         # self.window_size = 10000000  # slip on the chrom

#         self.min_align_length = min_align_length

#         self.min_mapq = min_mapq

#         # self.min_repeat_overlap_ratio = 50

#         # self.min_align_num = 3

#         # self.output_seq = False

#         # self.rm_TE = True
#         self.all_possible_chrs = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
#                              "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21",
#                              "chr22", "chrX"]

class Align:

    def __init__(self, read_name, chr, ref_start, ref_end, read_start, read_end, length, direct, seq, align_type, filter_flag, mapq, mismatch_ratio, isize=-1):
        self.read_name = read_name
        self.chr = chr
        self.ref_start = ref_start
        self.ref_end = ref_end

        self.read_start = read_start
        self.read_end = read_end

        self.length = length
        self.direct = direct
        self.seq = seq

        self.align_type = align_type
        self.filter_flag = filter_flag
        self.mapq = mapq
        self.mismatch_ratio = mismatch_ratio

        self.index = -1
        self.isize = isize

    def set_insert_size(self, isize):
        self.isize = isize

    def set_index(self, index):
        self.index = index

    def to_string(self):
        return "{}-{}-{}-{}".format(self.chr, self.ref_start, self.ref_end, self.direct)
    
class Bkp:

    def __init__(self, chrom, pos, direct, type):
        self.chrom = chrom
        self.pos = pos
        self.direct = direct
        self.type = type

    def to_string(self):
        return "{}:{}-{}".format(self.chrom, self.pos, self.direct)

    # rp_bkp is not accurate while sp_bkp is accurate, use sp_bkp if accessed and use the mean when both sp_bkp, use the nearest when both rp_bkp
    # tgs_bkp use the mean
    def merge(self, target_bkp):
        """
        Merge two bkps, calculate avg cord
        :param target_bkp:
        :return:
        """
        if self.type=="ngs_rp" and target_bkp.type=="ngs_sp":
            self.type = "ngs_sp"
            self.pos = target_bkp.pos
        elif self.type=="ngs_sp" and target_bkp.type=="ngs_rp":
            self.pos = self.pos
        # consider the direct
        elif self.type == "ngs_rp" and target_bkp.type=="ngs_rp":
            if self.direct == "L":
                self.pos = min(self.pos, target_bkp.pos)
            elif self.direct == "S":
                self.pos = max(self.pos, target_bkp.pos)
        else: # both ngs_sp or both tgs or both end
            self.pos = int((self.pos + target_bkp.pos) / 2)

class Junction:

    def __init__(self, junc_type, junc_sub_type, support_reads, include_bkps, updated_support_reads=0):
        self.junc_type = junc_type
        self.junc_sub_type = junc_sub_type
        self.support_reads = support_reads
        self.include_bkps = include_bkps

        self.updated_support_reads = updated_support_reads

    def to_string(self):
        return "{}-{}-{}".format(self.include_bkps[0].to_string(), self.junc_sub_type, self.include_bkps[1].to_string())

    def merge(self, target_junc):
        """
        Merge two junc by merge each bkp
        :param target_junc:
        :return:
        """
        for i in range(len(self.include_bkps)):
            self.include_bkps[i].merge(target_junc.include_bkps[i])

        for sr in target_junc.support_reads:
            self.support_reads.append(sr)

class Junction_Series:

    def __init__(self, include_juncs, support_reads, end):
        self.include_juncs = include_juncs
        self.support_reads = support_reads

        self.end = end

    def to_string(self):
        string = ""
        string = "{}:{}_".format(self.end[0].chrom,self.end[0].pos)
        string += "{}:{}".format(self.include_juncs[0].include_bkps[0].chrom,self.include_juncs[0].include_bkps[0].pos)
        for i in range(self.get_junc_num()-1):
            string += "-{}:{}".format(self.include_juncs[i].include_bkps[1].chrom,self.include_juncs[i].include_bkps[1].pos)
            string += "_{}:{}".format(self.include_juncs[i+1].include_bkps[0].chrom,self.include_juncs[i+1].include_bkps[0].pos)
        string += "-{}:{}".format(self.include_juncs[self.get_junc_num()-1].include_bkps[1].chrom,self.include_juncs[self.get_junc_num()-1].include_bkps[1].pos)
        string += "_{}:{}".format(self.end[1].chrom,self.end[1].pos)
        return string

    def get_junc(self, index):
        return self.include_juncs[index]

    def get_junc_num(self):
        return len(self.include_juncs)

    def get_support_num(self):
        num = 0
        for support_read in self.support_reads:
            num += len(support_read)
        return num

    def get_include_juncs_minus(self):
        include_juncs_minus = []
        support_reads_minus = []
        for i in range(self.get_junc_num()-1,-1,-1):
            include_juncs_minus.append(Junction(self.include_juncs[i].junc_type,self.include_juncs[i].junc_sub_type,self.include_juncs[i].support_reads,[self.include_juncs[i].include_bkps[1],self.include_juncs[i].include_bkps[0]]))
            support_reads_minus.append(self.include_juncs[i].support_reads)
        end_minus = [self.end[1],self.end[0]]
        return Junction_Series(include_juncs_minus, support_reads_minus, end_minus)
    
    def merge(self, target_juncSe, overlap_start_on_self, overlap_start_on_target, overlap_num):
        """
        merge target juncSe to self
        :param target_juncSe:
        :return:
        """
        # updata each junc
        for i in range(0, overlap_num):
            self.include_juncs[overlap_start_on_self + i].merge(target_juncSe.include_juncs[overlap_start_on_target + i])

        for i in range(overlap_start_on_target - 1, -1, -1):
            self.include_juncs.insert(0, target_juncSe.include_juncs[i])

        for i in range(overlap_start_on_target + overlap_num, target_juncSe.get_junc_num()):
            self.include_juncs.append(target_juncSe.include_juncs[i])
        
        self.update_support_reads()

        if overlap_start_on_target==0 and overlap_start_on_self==0:
            if self.end[0].direct=='L':
                self.end[0] = self.end[0] if self.end[0].pos<target_juncSe.end[0].pos else target_juncSe.end[0]
            else:
                self.end[0] = self.end[0] if self.end[0].pos>target_juncSe.end[0].pos else target_juncSe.end[0]
        elif overlap_start_on_target==0:
            self.end[0] = self.end[0]
        else:
            self.end[0] = target_juncSe.end[0]
        # self.support_reads.extend(target_juncSe.support_reads) # 这个应该不对，如果是按照每个节点获取support reads的话，JuncSe的support reads不能这么直接弄
        if overlap_start_on_target + overlap_num==target_juncSe.get_junc_num() and overlap_start_on_self + overlap_num==self.get_junc_num():
            if self.end[1].direct=='L':
                self.end[1] = self.end[1] if self.end[1].pos<target_juncSe.end[1].pos else target_juncSe.end[1]
            else:
                self.end[1] = self.end[1] if self.end[1].pos>target_juncSe.end[1].pos else target_juncSe.end[1]
        elif overlap_start_on_target + overlap_num==target_juncSe.get_junc_num():
            self.end[1] = self.end[1]
        else:
            self.end[1] = target_juncSe.end[1]

    def combine(self, next_juncSe):
        for junc in next_juncSe.include_juncs:
            self.include_juncs.append(junc)
        
        for reads in next_juncSe.support_reads:
            self.support_reads.append(reads)
        
        self.end[1] = next_juncSe.end[1]

    def update_support_reads(self):
        new_support_reads = []
        for junc in self.include_juncs:
            new_support_reads.append(junc.support_reads)
        self.support_reads = new_support_reads

# class Multipul_Junc:

#     def __init__(self, same_bkp_info, overlap, support_junc_series, support_reads):
#         self.same_bkp_info = same_bkp_info
#         self.overlap = overlap
#         self.support_junc_series = support_junc_series
#         self.support_reads = support_reads

#         self.other_bkp_support = None

#         self.merge_juncs = None
#         self.diff_juncs = []

#     def get_support_reads_num(self):
#         return len(self.support_reads)

#     def judge_minus(self):
#         for i in range(len(self.same_bkp_info)):
#             if self.same_bkp_info[i][1] != 1:
#                 return False
#         support_junc_series_minus = []
#         for i in range(len(self.support_junc_series)):
#             support_junc_series_minus.append(self.support_junc_series[i].get_include_juncs_minus())
#         self.support_junc_series = support_junc_series_minus
#         return True
    
#     def merge_inner(self):
#         index_1 = self.same_bkp_info[0][0]
#         start_1 = index_1 - self.overlap
#         index_2 = self.same_bkp_info[1][0]
#         start_2 = index_2 - self.overlap
#         if self.overlap > 1:
#             new_end_1 = self.support_junc_series[0].include_juncs[index_1].include_bkps[0]
#             new_end_2 = self.support_junc_series[1].include_juncs[index_2].include_bkps[0]
#             new_junse1 = Junction_Series(self.support_junc_series[0].include_juncs[:index_1], self.support_reads[0], self.support_junc_series[0].whole_segments_type[:index_1-1], 
#                                         self.support_junc_series[0].segments_direct[:index_1+1], [self.support_junc_series[0].end[0], new_end_1])
#             new_junse2 = Junction_Series(self.support_junc_series[1].include_juncs[:index_2], self.support_reads[1], self.support_junc_series[1].whole_segments_type[:index_2-1], 
#                                         self.support_junc_series[1].segments_direct[:index_2+1], [self.support_junc_series[1].end[0], new_end_2])
#             new_junse1.merge(new_junse2, start_1, start_2, self.overlap)
#             self.merge_juncs = new_junse1
#         else:
#             if start_1 == 0 and start_2 != 0:
#                 new_end_2 = self.support_junc_series[1].include_juncs[index_2].include_bkps[0]
#                 new_junse2 = Junction_Series(self.support_junc_series[1].include_juncs[:index_2], self.support_reads[1], self.support_junc_series[1].whole_segments_type[:index_2-1], 
#                                         self.support_junc_series[1].segments_direct[:index_2+1], [self.support_junc_series[1].end[0], new_end_2])
#                 self.merge_juncs = new_junse2

#         self.diff_juncs.append(self.support_junc_series[0].include_juncs[index_1])
#         self.diff_juncs.append(self.support_junc_series[1].include_juncs[index_2])


#     def merge(self, another, merge_juncs_info, merge_info, merge_type): #merge_type: 0-two multi_juncs; 1-one juncse and one multi_juncs
#         if merge_type==0:
#             if merge_juncs_info != None:
#                 if merge_juncs_info[0] != -1:
#                     self.merge_juncs.merge(another.merge_juncs, merge_juncs_info[0], merge_juncs_info[1], merge_juncs_info[2]) #这里现在相当于是两个juncSe的merge，所以还得要相应的参数，而且要考虑两个是[]的情况
#                 else:
#                     if self.merge_juncs == None:
#                         self.merge_juncs = another.merge_juncs
#             # merge_info: (index, merge_index) 
#             for info in merge_info:
#                 if info[1] == -1:
#                     self.same_bkp_info.append(another.same_bkp_info[info[0]])
#                     self.support_junc_series.append(another.support_junc_series[info[0]])
#                     self.diff_juncs.append(another.diff_juncs[info[0]])
#                     self.support_reads.append(another.support_reads[info[0]])
#                 else:
#                     self.diff_juncs[info[1]].merge(another.diff_juncs[info[0]])
#                     self.support_reads[info[1]].extend(another.support_reads[info[0]])
#         else:
#             if merge_juncs_info != None:
#                 index = merge_info[0]
#                 new_end = another.include_juncs[index].include_bkps[1]
#                 new_juncse = Junction_Series(another.include_juncs[:index], another.support_reads, another.whole_segments_type[:index-1], 
#                                                 another.segments_direct[:index+1], [another.end[0], new_end])
#                 if merge_juncs_info[0] != -1:
#                     self.merge_juncs.merge(new_juncse, merge_juncs_info[0], merge_juncs_info[1], merge_juncs_info[2])
#                 else:
#                     if self.merge_juncs == None:
#                         self.merge_juncs = new_juncse
#             if merge_info[1] == -1:
#                 self.same_bkp_info.append((merge_info[0], 0))
#                 self.support_junc_series.append(another)
#                 self.diff_juncs.append(another.include_juncs[merge_info[0]])
#                 self.support_reads.append(another.support_reads)
#             else:
#                 self.diff_juncs[merge_info[1]].merge(another.include_juncs[merge_info[0]])
#                 self.support_reads[merge_info[1]].extend(another.support_reads)

#     def get_same_bkp(self):
#         same_bkp = None
#         if len(self.diff_juncs)>0:
#             same_bkp = self.diff_juncs[0].include_bkps[0]
#             for i in range(1, len(self.diff_juncs)):
#                 same_bkp.merge(self.diff_juncs[i].include_bkps[0])

#         return same_bkp
    
#     def get_filted_result(self):
#         if self.other_bkp_support != None:
#             return self.support_junc_series[self.other_bkp_support]
#         else:
#             max_support_reads_index = -1
#             max_support_reads_num = 0
#             for i in range(len(self.support_junc_series)):
#                 if len(self.support_reads[i])>max_support_reads_num:
#                     max_support_reads_num = len(self.support_reads[i])
#                     max_support_reads_index = i

#             return self.support_junc_series[max_support_reads_index]
    
#     def to_string(self):
#         string = ""
#         index = []
#         start = []
#         for i in range(len(self.same_bkp_info)):
#             index.append(self.same_bkp_info[i][0])
#             start.append(self.same_bkp_info[i][0]-self.overlap)
#         for i in range(len(index)):
#             string += "juncse:{}, support_reads_num:{}\n".format(self.support_junc_series[i].to_string(), len(self.support_reads[i]))
#             string += "index:{}, start:{}\n".format(index[i], start[i])
#         same_bkp = self.get_same_bkp()
#         string += "same_bkp:{}\n".format(same_bkp.to_string())
#         string += "different_bkps:\n"
#         for i in range(len(self.diff_juncs)):
#             string += self.diff_juncs[i].include_bkps[1].to_string() + '\n'
#         return string

class Single_same_Multi_junc:
    def __init__(self, same_bkp_info, juncs_info, support_reads_num, seqs):
        self.same_bkp_info = same_bkp_info
        self.juncs_info = juncs_info
        self.support_reads_num = support_reads_num
        self.seqs = seqs
    
    def get_support_num(self):
        return len(self.support_reads)

    def add(self, new_junc, new_seq, new_support_reads_num):
        self.juncs_info.append(new_junc)
        self.seqs.append(new_seq)
        self.support_reads_num = self.support_reads_num + new_support_reads_num
    
    def to_string(self):
        string = "{}-[".format(self.same_bkp_info.to_string())
        for i in range(len(self.juncs_info)-1):
            string += self.juncs_info[i].to_string()
            string += '/'
        string += '{}]'.format(self.juncs_info[len(self.juncs_info)-1].to_string())
        if len(self.seqs) != 0:
            string += '\n'
            string = string + self.seqs[0] + '\n'
            string = string + self.seqs[1]
        return string

class Bridge_Candidate:
    def __init__(self, juncSe1_index, junc1_index, bkp1_index, juncSe2_index, junc2_index, bkp2_index, juncSe_list):
        self.juncSe1_index = juncSe1_index
        self.junc1_index = junc1_index
        self.bkp1_index = bkp1_index
        self.juncSe2_index = juncSe2_index
        self.junc2_index = junc2_index
        self.bkp2_index = bkp2_index
        self.juncSe_list = juncSe_list

        self.type = None
    
    def judge_bridge_type(self):
        if self.juncSe_list[self.juncSe1_index].include_juncs[self.junc1_index].include_bkps[self.bkp1_index].direct == 'S':
            pos1 = self.juncSe_list[self.juncSe1_index].include_juncs[self.junc1_index].include_bkps[self.bkp1_index].pos
            pos2 = self.juncSe_list[self.juncSe2_index].include_juncs[self.junc2_index].include_bkps[self.bkp2_index].pos
        else:
            pos1 = self.juncSe_list[self.juncSe2_index].include_juncs[self.junc2_index].include_bkps[self.bkp2_index].pos
            pos2 = self.juncSe_list[self.juncSe1_index].include_juncs[self.junc1_index].include_bkps[self.bkp1_index].pos
        if 0 < pos1 - pos2 < 1000:
            self.type = "DUP-bridge"
        elif 0 < pos2 - pos1 < 1000:
            self.type = "DEL-bridge"
        else:
            self.type = None
    
    def to_string(self):
        junc1_string_info_init = self.juncSe_list[self.juncSe1_index].include_juncs[self.junc1_index].to_string().split('-')
        if self.bkp1_index == 0:
            junc1_string_info = [junc1_string_info_init[3], junc1_string_info_init[4], junc1_string_info_init[0], junc1_string_info_init[1]]
            junc1_string = '-'.join(junc1_string_info)
        else:
            junc1_string_info = [junc1_string_info_init[0], junc1_string_info_init[1], junc1_string_info_init[3], junc1_string_info_init[4]]
            junc1_string = '-'.join(junc1_string_info)
        junc2_string_info_init = self.juncSe_list[self.juncSe2_index].include_juncs[self.junc2_index].to_string().split('-')
        if self.bkp2_index == 1:
            junc2_string_info = [junc2_string_info_init[3], junc2_string_info_init[4], junc2_string_info_init[0], junc2_string_info_init[1]]
            junc2_string = '-'.join(junc2_string_info)
        else:
            junc2_string_info = [junc2_string_info_init[0], junc2_string_info_init[1], junc2_string_info_init[3], junc2_string_info_init[4]]
            junc2_string = '-'.join(junc2_string_info)
        string = "{}_{}_{}".format(junc1_string, self.type, junc2_string)
        return string
    
class EndsForCombine:
    def __init__(self, chrom, pos, direct, last_bkp_pos, juncSe_index, end_index):
        self.chrom = chrom
        self.pos = pos
        self.direct = direct
        self.last_bkp_pos =  last_bkp_pos
        self.juncSe_index = juncSe_index
        self.end_index = end_index
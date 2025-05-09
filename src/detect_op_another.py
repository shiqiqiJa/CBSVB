from src.classes_another import Bkp, Junction, Align

def detect_intra_inter_trans(srt_read_aligns, options, repeat_regions, pm_seq, non_local_thresh=10000):
    """
    Given read aligns (sort by read cord), apply rules to determine if ther are inter (non local)
    """

    read_include_juncs = []
    end1 = None
    end2 = None
    end = [end1, end2]

    # # detect from each two aligns
    for i in range(len(srt_read_aligns) - 1):
        flag_repeat = 0

        cur_align = srt_read_aligns[i]
        next_align = srt_read_aligns[i + 1]
        
        if cur_align.chr not in options.all_possible_chrs or next_align.chr not in options.all_possible_chrs:
            return [], end
        
        junc_sub_type = []
        bkp_type = ""
        
        if options.seq_type == "ngs":
            # decide the bkp type ([4]inside must be 1 or 2，we have format 'read1','read2','read1_sa','read2_sa')
            if cur_align.align_type[4] == next_align.align_type[4]:
                bkp_type = "ngs_sp"
            else:
                bkp_type = "ngs_rp"
        else:
            bkp_type = "tgs"

        # # Generate bkps and juncHT or HH trans by direct
        if cur_align.direct == 'T' and next_align.direct == 'T':
            bkp1 = Bkp(cur_align.chr, cur_align.ref_end, "S", bkp_type)
            bkp2 = Bkp(next_align.chr, next_align.ref_start, "L", bkp_type)

        elif cur_align.direct == 'F' and next_align.direct == 'F':
            bkp1 = Bkp(cur_align.chr, cur_align.ref_start, "L", bkp_type)
            bkp2 = Bkp(next_align.chr, next_align.ref_end, "S", bkp_type)

        elif cur_align.direct == 'T' and next_align.direct == 'F':
            bkp1 = Bkp(cur_align.chr, cur_align.ref_end, "S", bkp_type)
            bkp2 = Bkp(next_align.chr, next_align.ref_end, "S", bkp_type)

        elif cur_align.direct == 'F' and next_align.direct == 'T':
            bkp1 = Bkp(cur_align.chr, cur_align.ref_start, "L", bkp_type)
            bkp2 = Bkp(next_align.chr, next_align.ref_start, "L", bkp_type)

        else:
            bkp1 = None
            bkp2 = None

        for item in repeat_regions[bkp1.chrom]:
            if item[0] < bkp1.pos < item[1]:
                flag_repeat += 1
                break
        for item in repeat_regions[bkp2.chrom]:
            if item[0] < bkp2.pos < item[1]:
                flag_repeat += 1
                break
        if flag_repeat == 1:
            return [], end
        
        if i==0:
            pos_end1 = cur_align.ref_end if bkp1.pos==cur_align.ref_start else cur_align.ref_start
            direct_end1 = 'L' if bkp1.direct=='S' else 'S'
            end1 = Bkp(cur_align.chr,pos_end1,direct_end1,"end")
        if i==len(srt_read_aligns)-2:
            pos_end2 = next_align.ref_end if bkp2.pos==next_align.ref_start else next_align.ref_start
            direct_end2 = 'L' if bkp2.direct=='S' else 'S'
            end2 = Bkp(next_align.chr,pos_end2,direct_end2,"end")
        end = [end1, end2]

        # # different chrom, consider inter trans
        if cur_align.chr != next_align.chr:
            # print("cur_align.chr:{}, next_align.chr:{}\n".format(cur_align.chr, next_align.chr))
            junc_type = "inter"
            junc_sub_type.append("interTrans")
            # append to list
            read_include_juncs.append(Junction(junc_type, '-'.join(junc_sub_type), [[cur_align, next_align]], [bkp1, bkp2]))

        # # same chrom
        else:
            junc_type = 'intra'

            # consider large gap on read (INS)
            # if next_align.read_start - cur_align.read_end >= non_local_thresh:
            if next_align.read_start - cur_align.read_end > 10000:
                junc_sub_type.append("INS")
                align_obj = Align(cur_align.read_name, None, None, None, cur_align.read_end, next_align.read_start, next_align.read_start-cur_align.read_end, None, pm_seq[cur_align.read_end:next_align.read_start], "AP", 'PASS', 60, 0)
                read_include_juncs.append(Junction(junc_type, "-".join(junc_sub_type), [[cur_align, align_obj, next_align]], [bkp1, bkp2]))
            # consider large gap or overlap on ref (DEL or DUP)
            elif  next_align.ref_start - cur_align.ref_end >= non_local_thresh or cur_align.ref_start - next_align.ref_end >= non_local_thresh:
                junc_sub_type.append("intraTrans")
                read_include_juncs.append(Junction(junc_type, "-".join(junc_sub_type), [[cur_align, next_align]], [bkp1, bkp2]))

    return read_include_juncs, end

def judge_overlap(pos1, pos2):
    if pos2[0] - pos1[1] > 1000 or pos1[0] - pos2[1] > 1000:
        return False
    return True

def judge_noisy_aligns(srt_read_aligns, min_align_num=2):
    noise_segs_chroms = {}
    for i in range(len(srt_read_aligns) - 1):
        flag = 0
        cur_align = srt_read_aligns[i]
        chrom = cur_align.chr
        if chrom not in noise_segs_chroms.keys():
            noise_segs_chroms[chrom] = []
        start = cur_align.ref_start 
        end = cur_align.ref_end
        for j in range(len(noise_segs_chroms[chrom])):
            if judge_overlap((start, end), (noise_segs_chroms[chrom][j][0], noise_segs_chroms[chrom][j][1])):
                noise_segs_chroms[chrom][j][2] += 1
                flag = 1
                break
        if flag == 0:
            noise_segs_chroms[chrom].append([start, end, 1])
    for chrom in noise_segs_chroms.keys():
        for seg in noise_segs_chroms[chrom]:
            if seg[2] >= min_align_num:
                return True
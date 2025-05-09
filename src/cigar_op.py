import re

def cigar_to_list(cigar):

    opVals = re.findall(r'(\d+)([\w=])', cigar) #其实就是按照“数字+字母”的格式进行提取，“\d”表示数字[0-9],"\w"表示数字、字母、下划线、中文
    lengths = []
    ops = []
    for i in range(len(opVals)):
        lengths.append(int(opVals[i][0]))
        ops.append(opVals[i][1])

    return ops, lengths

def cal_align_read_start_and_length_from_cigar(cigar_ops, cigar_lens, pm_direct, cur_direct):
    """
    same as function name
    :param cigar_ops:
    :param cigar_lens:
    :param direct:
    :param read_1or2:
    :return:
    """

    def cal_align_read_start(ops, lens):

        # find first "M" in CIGAR
        first_M_index = -1

        for i in range(len(ops)):
            if ops[i] == 'M':
                first_M_index = i
                break

        # sum up lens before first_M_index
        read_start = 0
        for i in range(first_M_index):
            read_start += lens[i]

        return read_start

    def cal_align_length(ops, lens):

        length = 0
        for i in range(len(ops)):
            if ops[i] == 'M':
                length += lens[i]
        return length

    if cur_direct == pm_direct:
        # calculate sa align's read start by find the first 'M' in CIGAR
        align_read_start = cal_align_read_start(cigar_ops, cigar_lens)

        # calculate sa align's length
        align_length = cal_align_length(cigar_ops, cigar_lens)
    else:
        cigar_ops_reverse = [cigar_ops[i] for i in range(len(cigar_ops) - 1, -1, -1)]
        cigar_lens_reverse = [cigar_lens[i] for i in range(len(cigar_lens) - 1, -1, -1)]

        # calculate sa align's read start by find the first 'M' in CIGAR
        align_read_start = cal_align_read_start(cigar_ops_reverse, cigar_lens_reverse)

        # calculate sa align's length
        align_length = cal_align_length(cigar_ops_reverse, cigar_lens_reverse)
        
    return align_read_start, align_length
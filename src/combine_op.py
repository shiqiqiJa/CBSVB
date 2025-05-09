import networkx as nx
import os

def get_chrom_num(end):
    if end.chrom[3:] == 'X':
        return 23
    else:
        return int(end.chrom[3:])

def judge_combine(bkp1, bkp2, options):
    if bkp1.chrom != bkp2.chrom:
        return False
    if abs(bkp1.pos-bkp2.pos)>options.combine_length_diff:
        return False
    if bkp1.direct == bkp2.direct:
        return False
    return True

 #这个参数决定了最终能combine的结果，针对不同的模拟情况，这个参数的选择会导致不同的结果
# 因为这个数据模拟的时候reciprocal translocation过短，所以这里选择5000比较合适

def judge_in_out(juncSe1, juncSe2, options):
    # 判断juncS1和juncSe2的关系(in表示juncSe2能连到juncSe1, out表示juncSe1能连到juncSe2)
    for i in range(2):
        for j in range(2):
            if juncSe1.end[i].chrom == juncSe2.end[j].chrom and juncSe1.end[i].direct != juncSe2.end[j].direct:
                if abs(juncSe1.end[i].pos-juncSe2.end[j].pos)>options.combine_length_diff:
                    continue
                if juncSe1.end[i].direct == 'L':
                    if j==0:
                        lastbkp = juncSe2.include_juncs[0].include_bkps[0].pos
                    else:
                        lastbkp = juncSe2.include_juncs[len(juncSe2.include_juncs)-1].include_bkps[1].pos
                    if i==0:
                        another_lastbkp = juncSe1.include_juncs[0].include_bkps[0].pos
                    else:
                        another_lastbkp = juncSe1.include_juncs[len(juncSe1.include_juncs)-1].include_bkps[1].pos
                    # print(lastbkp, another_lastbkp)
                    # if abs(lastbkp-another_lastbkp)<1000:
                    #     continue
                    if juncSe1.end[i].pos - another_lastbkp < -1000 or juncSe2.end[j].pos - lastbkp > 1000:
                        continue
                    if lastbkp-1000 <= juncSe1.end[i].pos <= another_lastbkp+options.combine_length_diff:
                        return 'in'
                else:
                    if i==0:
                        lastbkp = juncSe1.include_juncs[0].include_bkps[0].pos
                    else:
                        lastbkp = juncSe1.include_juncs[len(juncSe1.include_juncs)-1].include_bkps[1].pos
                    if j==0:
                        another_lastbkp = juncSe2.include_juncs[0].include_bkps[0].pos
                    else:
                        another_lastbkp = juncSe2.include_juncs[len(juncSe2.include_juncs)-1].include_bkps[1].pos
                    if juncSe1.end[i].pos - another_lastbkp > 1000 or juncSe2.end[j].pos - lastbkp < -1000:
                        continue
                    if lastbkp-1000 <=  juncSe2.end[j].pos <= another_lastbkp+options.combine_length_diff:
                        return 'out'
    
    return None

def get_init_link(srt_juncSe_list, options):
    # 将juncSe都调整为方向是小染色体到大染色体或同一染色体上小位置到大位置
    for i in range(len(srt_juncSe_list)):
        chrom1 = get_chrom_num(srt_juncSe_list[i].end[0])
        chrom2 = get_chrom_num(srt_juncSe_list[i].end[1])
        if chrom1 > chrom2 or (chrom1 == chrom2 and srt_juncSe_list[i].end[0].pos > srt_juncSe_list[i].end[1].pos):
            srt_juncSe_list[i] = srt_juncSe_list[i].get_include_juncs_minus()
    srt_juncSe_list_again = sorted(srt_juncSe_list, key=lambda x:(x.end[0].chrom, x.end[0].pos))

    # for i in range(len(srt_juncSe_list_again)):
    #     print(i, srt_juncSe_list_again[i].to_string())
    
    G = nx.DiGraph()
    for i in range(len(srt_juncSe_list_again)):
        G.add_node(i)

    for i in range(len(srt_juncSe_list_again)-1):
        if G.degree(i) == 2:
            continue
        if 'INS' in srt_juncSe_list_again[i].to_string().split('-'):
            continue
        for j in range(i+1, len(srt_juncSe_list_again)):
            if  G.degree(j) == 2 or G.degree(i) == 2:
                continue
            if 'INS' in srt_juncSe_list_again[j].to_string().split('-'):
                continue
            if judge_in_out(srt_juncSe_list_again[j], srt_juncSe_list_again[i], options) == 'out':
                G.add_edge(j, i)
            elif judge_in_out(srt_juncSe_list_again[j], srt_juncSe_list_again[i], options) == 'in':
                G.add_edge(i, j)
    return G, srt_juncSe_list_again

def get_linked_list(srt_juncSe_list, options):
    juncSe_linked_list = []

    G, srt_juncSe_list_again = get_init_link(srt_juncSe_list, options)    
    # nx.draw_networkx(G)

    # 下面这一段整体还要再考虑，有可能一个点有两个"in"或两个"out"！
    for c in nx.weakly_connected_components(G):
        if len(c) == 1:
            juncSe_linked_list.append(list(c))
        else:
            juncSe_order = []
            c = sorted(c, key=lambda item:G.degree(item))
            now_point = c[0]
            juncSe_order.append(c[0])
            while True:
                if G.in_degree(now_point) != 0:
                    for edge in G.in_edges(now_point):
                        if edge[0] not in juncSe_order:
                            juncSe_order.append(edge[0])
                            now_point = edge[0]
                            break
                if G.out_degree(now_point) != 0:
                    for edge in G.out_edges(now_point):
                        if edge[1] not in juncSe_order:
                            juncSe_order.append(edge[1])
                            now_point = edge[1]
                            break
                if now_point == c[1]:
                    juncSe_linked_list.append(juncSe_order)
                    break
    return juncSe_linked_list, srt_juncSe_list_again

def combine_juncsSe(juncSe_list, juncSe_linked_list, options):
    combined_juncSe_list = []
    for order in juncSe_linked_list:
        current_juncSe = juncSe_list[order[0]]
        for i in range(1, len(order)):
            if not judge_combine(current_juncSe.end[1], juncSe_list[order[i]].end[0], options):
                if judge_combine(current_juncSe.end[1], juncSe_list[order[i]].end[1], options):
                    juncSe_list[order[i]] = juncSe_list[order[i]].get_include_juncs_minus()
                elif judge_combine(current_juncSe.end[0], juncSe_list[order[i]].end[0], options):
                    current_juncSe = current_juncSe.get_include_juncs_minus()
                elif judge_combine(current_juncSe.end[0], juncSe_list[order[i]].end[1], options):
                    juncSe_list[order[i]] = juncSe_list[order[i]].get_include_juncs_minus()
                    current_juncSe = current_juncSe.get_include_juncs_minus()
            current_juncSe.combine(juncSe_list[order[i]])
        combined_juncSe_list.append(current_juncSe)
        
    return combined_juncSe_list

def report_result_juncs(sample_name, juncSe_list, juncSe_linked_list, options):
    with open(os.path.join(options.out_path, "{}_juncs.txt".format(sample_name)), "w") as output_file:
        output_file.write('#id\tsupport_reads\tdescription\n')
        num = 0
        if juncSe_linked_list == None:
            for i in range(len(juncSe_list)):
                if i < len(juncSe_list)-1:
                    for j in range(juncSe_list[i].get_junc_num()):
                        output_file.write('Juncs.{}\t{}\t{}\n'.format(num, len(juncSe_list[i].support_reads), juncSe_list[i].include_juncs[j].to_string()))
                        num += 1
                else:
                    for j in range(juncSe_list[i].get_junc_num()):
                        if j < juncSe_list[i].get_junc_num() -1:
                            output_file.write('Juncs.{}\t{}\t{}\n'.format(num, len(juncSe_list[i].support_reads), juncSe_list[i].include_juncs[j].to_string()))
                            num += 1
                        else:
                            output_file.write('Juncs.{}\t{}\t{}'.format(num, len(juncSe_list[i].support_reads), juncSe_list[i].include_juncs[j].to_string()))
        else:
            for order in juncSe_linked_list:
                support_reads = sum([juncSe_list[index].get_support_num() for index in order])
                for index in order:
                    for k in range(len(juncSe_list[index].include_juncs)):
                        junc = juncSe_list[index].include_juncs[k]
                        if k==0 or k==len(juncSe_list[index].include_juncs)-1:
                            output_file.write('Juncs.{}\t{}\t{}\n'.format(num, support_reads, junc.to_string()))
                        else:
                            output_file.write('Juncs.{}\t{}\t{}\n'.format(num, len(juncSe_list[index].support_reads), junc.to_string()))
                        num += 1

    output_file.close()

def combine_results(srt_juncSe_list, options):
    juncSe_linked_list, srt_juncSe_list_again = get_linked_list(srt_juncSe_list, options)
    combined_juncSe_list = combine_juncsSe(srt_juncSe_list_again, juncSe_linked_list, options)
    report_result_juncs('combined_junc_results', srt_juncSe_list_again, juncSe_linked_list, options)
    return combined_juncSe_list


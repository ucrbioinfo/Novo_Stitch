import sys
#key method to select best hit for stitiching c1 and c2
def onepair(c1, c2, contig_len, candidate_hits, alm_len_cutoff, alm_shift_ratio):
    contig1_id = c1[1]
    contig1_name = "N"+contig1_id
    contig1_ori = c1[4]
    contig2_id = c2[1]
    contig2_name = "N"+contig2_id
    contig2_ori = c2[4]

    #check if the 2 contigs overlap
    contig1_begin = c1[2]
    contig1_end = c1[3]
    contig2_begin = c2[2]
    contig2_end = c2[3]
    n1 = contig_len[contig1_name]
    n2 = contig_len[contig2_name]
    if contig1_end <= contig2_begin:
        print "No overlap between", [contig1_name, contig2_name]			
        best_line = contig2_name+"\t"+contig1_name+"\tGAP"
        return best_line

    #store the correct hits in smaller_candi_hits
    smaller_candi_hits = []
    for line in candidate_hits:
        line = line.strip()
        cols = line.split()
        q_name = cols[0]
        r_name = cols[1]
        if contig1_name != r_name or contig2_name != q_name:
            continue
        c2_begin = q_begin = int(cols[6])
        c2_end = q_end = int(cols[7])
        c1_begin = r_begin = int(cols[8])
        c1_end = r_end = int(cols[9])

        if (contig1_ori == contig2_ori and r_begin >= r_end) or (contig1_ori != contig2_ori and r_begin <= r_end):  #check if orientation is correct
            continue
        alm_len = int(cols[3])
        if alm_len < alm_len_cutoff: #10000:
            continue
    
        if contig1_ori == "-":
            aj_c1_begin = min(n1+1-c1_begin, n1+1-c1_end)
            aj_c1_end = max(n1+1-c1_begin, n1+1-c1_end)
        else:
            aj_c1_begin = min(c1_begin, c1_end)
            aj_c1_end = max(c1_begin, c1_end)
        if contig2_ori == "-":
            aj_c2_begin = min(n2+1-c2_begin, n2+1-c2_end)
            aj_c2_end = max(n2+1-c2_begin, n2+1-c2_end)
        else:
            aj_c2_begin = min(c2_begin, c2_end)
            aj_c2_end = max(c2_begin, c2_end)
        opt_c1_begin = contig1_begin + ( aj_c1_begin - 1 )
        opt_c1_end = contig1_end - ( n1 - aj_c1_end )		
        opt_c2_begin = contig2_begin + ( aj_c2_begin - 1 )
        opt_c2_end = contig2_end - ( n2 - aj_c2_end )
        opt_c1_mid = float( opt_c1_begin + opt_c1_end )/2
        opt_c2_mid = float( opt_c2_begin + opt_c2_end )/2
        diff = opt_c1_mid - opt_c2_mid
        ratio = abs(diff)/min(n1,n2)
        if ratio > alm_shift_ratio: #0.05:
            continue
        smaller_candi_hits.append(line)
        
    #select the best hit from smaller_candi_hits
    if len(smaller_candi_hits) > 0:
        best_line = smaller_candi_hits[0]
    else:
        print "No appropriate alignment for stitiching", [contig1_name, contig2_name]			
        best_line = contig2_name+"\t"+contig1_name+"\tGAP"
        return best_line
    return best_line


def step1_best_hit(blast_file, list_file, contig_file, output_file, alm_len_cutoff, alm_shift_ratio):    
   
    #get contig sequences 
    contig_seq = {}
    contig_len = {}
    cur_contig = ""
    for line in open(contig_file):
        line = line.strip()
        if line[0] == '>':
            cols = line[1:].split()	
            contig_name = cols[0]
            cur_contig = contig_name
        else:
            seq = line 
            contig_seq[cur_contig] = seq
            contig_len[cur_contig] = len(seq)
    
    #get list_info 
    list_info = []
    for line in open(list_file):
        line = line.strip()
        cols = line.split()	
        list_info.append((cols[0], cols[1], int(cols[2]), int(cols[3]), cols[4]))
    list_info = sorted(list_info, key=lambda x:x[2])
    
    #store the top 10 hits of each type to candidate_hits to reduce the # of hits
    candidate_hits = []
    last_q_name = ""
    last_r_name = ""
    for line in open(blast_file):
        line = line.strip()
        cols = line.split()
        q_name = cols[0]
        r_name = cols[1]
        if q_name != last_q_name or r_name != last_r_name:
            candidate_hits.append(line)
            #print q_name, last_q_name
            left_times = 10
            last_q_name = q_name
            last_r_name = r_name
        else:
            if left_times > 0:
                candidate_hits.append(line)
                left_times -= 1
    
    fo = file(output_file, 'w')
    #select best hit from candidate_hits for each pair of contigs
    for i in range(0, len(list_info)-1):
        c1 = list_info[i]
        c2 = list_info[i+1]
        best_line = onepair(c1, c2, contig_len, candidate_hits, alm_len_cutoff,alm_shift_ratio)
        fo.write(best_line+"\n")
        
    
if __name__ == "__main__":
    step1_best_hit(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], len(sys.argv[5]), float(sys.argv[6]))
    
    

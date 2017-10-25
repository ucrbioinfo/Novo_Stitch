def conflict_or_not(alm_1, alm_2, qry_len, ref_len, orientation, ratio_cutoff): #return 1 means conflict, -1 no conflict 
    qry_1_begin = alm_1[6]
    qry_1_end = alm_1[7]
    ref_1_begin = alm_1[8]
    ref_1_end = alm_1[9]
    qry_2_begin = alm_2[6]
    qry_2_end = alm_2[7]
    ref_2_begin = alm_2[8]
    ref_2_end = alm_2[9]
    len_1 = alm_1[3]
    len_2 = alm_2[3]
    
    # rule of conflict between 2 alignments: the overlap of them (for either ref or qry) is larger than ratio_cutoff(e.g. 10%) * the length of shorter alignment 
    smaller_len = min(len_1, len_2)
    if orientation == "+":
        if qry_1_end - qry_2_begin < smaller_len * ratio_cutoff and ref_1_end - ref_2_begin < smaller_len * ratio_cutoff:
            return -1
        if qry_2_end - qry_1_begin < smaller_len * ratio_cutoff and ref_2_end - ref_1_begin < smaller_len * ratio_cutoff:
            return -1
    else:
        if qry_1_end - qry_2_begin < smaller_len * ratio_cutoff and ref_2_begin - ref_1_end < smaller_len * ratio_cutoff:
            return -1
        if qry_2_end - qry_1_begin < smaller_len * ratio_cutoff and ref_1_begin - ref_2_end < smaller_len * ratio_cutoff:
            return -1
    return 1

def optimal_alims(qry_name, ref_name, qry_len, ref_len,hits_list, orientation, length_cutoff, ratio_cutoff):
    # get candidate alignments stored in good_hits_list
    good_hits_list = []
    for line in hits_list:
        line = line.strip()
        cols = line.split()
        qry = cols[0]
        ref = cols[1]
        length = int(cols[3])
        r_begin = int(cols[8])
        r_end = int(cols[9])
        if qry != qry_name or ref != ref_name:
            continue
        if (orientation == "+" and r_begin > r_end) or (orientation == "-" and r_begin < r_end):
            continue
        if length < length_cutoff:
            continue
        good_hits_list.append((cols[0],cols[1],cols[2],int(cols[3]),cols[4],cols[5],int(cols[6]),int(cols[7]),int(cols[8]),int(cols[9]),cols[10],cols[11] ))
    good_hits_list = sorted(good_hits_list, key=lambda x:x[6])
    if len(good_hits_list) == 0: # no qualified alignments
        return [], -1

    # get conflicts matrix
    conflicts = [[0 for col in range(len(good_hits_list))] for row in range(len(good_hits_list))]	
    for i in range(0, len(good_hits_list)):
        for j in range(i, len(good_hits_list)):	
            conflicts[i][j] = conflict_or_not(good_hits_list[i], good_hits_list[j],qry_len, ref_len, orientation, ratio_cutoff)
    # dynamic programming algorithm for optimal subset of aligments
    S = [-1 for col in range(len(good_hits_list))] #S[i]: the score of optimal solution of subproblem 0~i if we pick alignment i
    prev = [-1 for col in range(len(good_hits_list))]
    S[0] = good_hits_list[0][3]
    prev[0] = -1
    for i in range(1, len(good_hits_list)):
        max_score = 0
        max_j = -1
        for j in range(0, i):
            if conflicts[j][i] == 1:
                continue
            if max_score < S[j]:
                max_score = S[j]
                max_j = j
        S[i] = max_score + good_hits_list[i][3]
        prev[i] = max_j
    total_score = max(S)
    last_j = S.index(total_score)
    j = last_j
    optimal_alm_list = []
    while j != -1:
        optimal_alm_list.append(good_hits_list[j])
        j = prev[j]	
    optimal_alm_list = sorted(optimal_alm_list, key=lambda x:x[6])

    return optimal_alm_list, total_score


def stitch_alims(qry_name, ref_name, qry_len, ref_len, hits_list, length_cutoff, ratio_cutoff):
    # get the optimal alignments set
    optimal_alm_list1, total_score1 = optimal_alims(qry_name, ref_name, qry_len, ref_len, hits_list, "+", length_cutoff, ratio_cutoff)
    optimal_alm_list2, total_score2 = optimal_alims(qry_name, ref_name, qry_len, ref_len, hits_list, "-", length_cutoff, ratio_cutoff)
    if total_score1 > total_score2:
        total_score = total_score1
        optimal_alm_list = optimal_alm_list1
    else:
        total_score = total_score2
        optimal_alm_list = optimal_alm_list2
    q_total_len = 0
    last_q_begin = -1
    last_q_end = -1

    # calculate confidence for optimal alignments set
    for hit in optimal_alm_list:
        q_begin = hit[6]
        q_end = hit[7]
        if q_begin > last_q_end:
            q_total_len += abs(q_end-q_begin+1)
        else:
            q_total_len += ( abs(q_end-q_begin+1) - abs(last_q_end-q_begin+1) )
        last_q_begin = q_begin
        last_q_end = q_end
    confidence = float(q_total_len)/qry_len
    return optimal_alm_list, confidence

    
def confidence(fout, mole_id, blast_file, contigs_file, merged_contigs_file, length_cutoff, ratio_cutoff):
    # some data structures
    hits_list = []
    for line in open(blast_file):
        line = line.strip()
        hits_list.append(line)
    
    qry_name_list = []
    ref_name_list = []
    for line in open(contigs_file):
        line = line.strip()
        if line[0] == '>':
            cols = line[1:].split()
            c_name = cols[0]
            c_len = int(cols[1][4:])
            qry_name_list.append((c_name, c_len))	
    for line in open(merged_contigs_file):
        line = line.strip()
        if line[0] == '>':
            cols = line[1:].split()
            c_name = cols[0]
            c_len = int(cols[1][4:])
            ref_name_list.append((c_name, c_len))	
    
    # for each contig
    for (q_name, q_len) in qry_name_list:
        max_confidence = -1
        max_r_name = ""
        # for each old contig, find the best stitched contig it belongs to 
        for (r_name, r_len) in ref_name_list:
            optimal_alm_list, confidence = stitch_alims(q_name, r_name, q_len, r_len, hits_list, length_cutoff, ratio_cutoff)	
            if confidence > max_confidence:
                max_confidence = confidence
                max_optimal_alm_list = optimal_alm_list
                max_r_name = r_name
        fout.write(mole_id + "\t" + q_name + "\t" + max_r_name + "\t" + str(max_confidence) + "\n")
        print q_name, "is mapped to", max_r_name, "with confidence =", max_confidence 
   

 

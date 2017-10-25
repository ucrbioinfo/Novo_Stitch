import sys

def revcomp(seq):
    rcseq = ""
    for i in range(0,len(seq)):
        c = seq[len(seq)-1-i]
        if c=='A':
            rcseq = rcseq + 'T'
        if c=='T':
            rcseq = rcseq + 'A'
        if c=='G':
            rcseq = rcseq + 'C'
        if c=='C':
            rcseq = rcseq + 'G'
        if c=='a':
            rcseq = rcseq + 't'
        if c=='t':
            rcseq = rcseq + 'a'
        if c=='g':
            rcseq = rcseq + 'c'
        if c=='c':
            rcseq = rcseq + 'g'
            
    return rcseq

def onepair(c1, c2, best_hits, contig_len):
    contig1_id = c1[1]
    contig1_name = "N"+contig1_id
    contig1_ori = c1[4]
    contig2_id = c2[1]
    contig2_name = "N"+contig2_id
    contig2_ori = c2[4]
    comb_name = contig1_name + "$" + contig2_name
    best_hit = best_hits[comb_name]
    cols = best_hit.split()
    if cols[2] == "GAP":
        return (-1, -1, -1, -1, -1, -1)

    c1_begin = int(cols[8])
    c1_end = int(cols[9])
    c2_begin = int(cols[6])
    c2_end = int(cols[7])
    n1 = contig_len[contig1_name]
    n2 = contig_len[contig2_name]

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

    return (aj_c1_begin, aj_c1_end, n1, aj_c2_begin, aj_c2_end, n2)
    #the begin, end positions of the alignment on contig1, contig2 (maybe on reverse complement seq) and lengths of contig1 and contig2



def step2_stitch(best_hit_file, list_file, contig_file, merged_contig_file):    
    
    #get the sequences of contigs
    contig_seqs = {}
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
            contig_seqs[cur_contig] = seq
            contig_len[cur_contig] = len(seq)
    #get the best hit for each pair of contigs
    best_hits = {}
    for line in open(best_hit_file):
        line = line.strip()
        if line == "":
            continue
        cols = line.split()
        contig1_name = cols[1]
        contig2_name = cols[0]
        comb_name = contig1_name + "$" + contig2_name
        best_hits[comb_name] = line
    #read list_file to get the coordinates on optical map for each contig
    list_info = []
    for line in open(list_file):
        line = line.strip()
        cols = line.split()	
        list_info.append((cols[0], cols[1], int(cols[2]), int(cols[3]), cols[4]))
    list_info.sort(key=lambda x:x[2]) #sort list_info by starting coordinates
    #calculate the beginning and ending positions of alignment on each pair of contigs
    pos_list = {}
    for i in range(0, len(list_info)-1):
        c1 = list_info[i]
        c2 = list_info[i+1]
        (aj_c1_begin, aj_c1_end, n1, aj_c2_begin, aj_c2_end, n2) = onepair(c1, c2, best_hits, contig_len)
        contig1_id = c1[1]
        contig1_name = "N"+contig1_id
        contig2_id = c2[1]
        contig2_name = "N"+contig2_id
        comb_name = contig1_name + "$" + contig2_name
        pos_list[comb_name] = (aj_c1_begin, aj_c1_end, n1, aj_c2_begin, aj_c2_end, n2)
    
    #merge seqs
    merged_seqs = {}
    included = set([])
    for i in range(0, len(list_info)-1):
        c1 = list_info[i]
        c2 = list_info[i+1]
        contig1_id = c1[1]
        contig1_name = "N"+contig1_id
        contig1_ori = c1[4]
        contig1_seq = contig_seqs[contig1_name]
        if contig1_ori == '-':
            contig1_seq = revcomp(contig1_seq)
        contig2_id = c2[1]
        contig2_name = "N"+contig2_id
        contig2_ori = c2[4]
        contig2_seq = contig_seqs[contig2_name]
        if contig2_ori == '-':
            contig2_seq = revcomp(contig2_seq)
        comb_name = contig1_name + "$" + contig2_name
        #merge contig1 and contig2               
        (aj_c1_begin, aj_c1_end, n1, aj_c2_begin, aj_c2_end, n2) = pos_list[comb_name]
        if (aj_c1_begin, aj_c1_end, n1, aj_c2_begin, aj_c2_end, n2) == (-1, -1, -1, -1, -1, -1):
            continue
        if contig1_id in included: # make sure each contig only appears in one stitched contig
            print [contig1_name,contig2_name], "are not stitched, because", contig1_name, "is already stitched"  
            continue
        contig1_left = contig1_seq[:aj_c1_begin-1]
        contig1_aligned = contig1_seq[aj_c1_begin-1: aj_c1_end]
        contig1_right = contig1_seq[aj_c1_end:]
        contig2_left = contig2_seq[:aj_c2_begin-1]
        contig2_aligned = contig2_seq[aj_c2_begin-1: aj_c2_end]
        contig2_right = contig2_seq[aj_c2_end:]
        merged_aligned = contig1_aligned #contig1 is from 5' and should have higher quality
        if len(contig1_left) >= len(contig2_left):
            merged_left = contig1_left
        else:
            merged_left = contig2_left
        if len(contig1_right) >= len(contig2_right):
            merged_right = contig1_right
        else:
            merged_right = contig2_right
        merged_seq = merged_left + merged_aligned + merged_right
        merged_seqs[comb_name] = merged_seq
        included.add(contig1_id)
        included.add(contig2_id)
        print [contig1_name,contig2_name], "are stitched"     
        
    # add the sequences of unstitched contigs to merged_seqs
    for i in range(0, len(list_info)):
        c = list_info[i]
        contig_id = c[1]
        contig_name = "N" + contig_id
        contig_ori = c[4]
        contig_seq = contig_seqs[contig_name]
        if contig_ori == '-':
            contig_seq = revcomp(contig_seq)
        if contig_id not in included:
            merged_seqs[contig_name] = contig_seq
    merged_seqs_list = sorted(merged_seqs.items())

    # write the fasta file
    stitch_info = []
    fo = file(merged_contig_file, 'w')
    molecule_id = list_info[0][0]
    for i in range(0, len(merged_seqs_list)):
        #write fasta file
        merged_contig_name = "M"+molecule_id+"_"+merged_seqs_list[i][0]
        title = ">"+merged_contig_name+" len="+str(len(merged_seqs_list[i][1]))
        fo.write(title+"\n")
        fo.write(merged_seqs_list[i][1]+"\n")
        #store stitch information
        cols = merged_seqs_list[i][0].split('$')
        info = molecule_id + "\t" + merged_contig_name + "\t"
        for c in cols:
            info += c 
            info += " "
        stitch_info.append(info)
    fo.close()

    return stitch_info

if __name__ == "__main__":
    step2_stitch(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
        

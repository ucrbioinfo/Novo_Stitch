import sys

def mtp_fasta(list_file, key_file, old_fasta_file, new_fasta_file): 
#obtain mtp fasta file 

    contig_title_id = {}
    for line in open(key_file):
        line = line.strip()
        if line[0] == '#' or line[0] == 'C':
            continue
        cols = line.split('\t')
        contig_id = cols[0]
        contig_title = cols[1] 
        contig_title_id[contig_title] = contig_id
    contig_ids = {}
    for line in open(list_file):
        line = line.strip()
        cols = line.split()
        contig_id = cols[1]	
        contig_ids[contig_id] = 0
    
    fo = file(new_fasta_file, 'w')
    
    cur_state = True
    for line in open(old_fasta_file):
        line = line.strip()
        if line[0] == '>':
            contig_title = line[1:]
            if contig_title not in contig_title_id:
                cur_state = False
                continue
            contig_id = contig_title_id[contig_title]
            if contig_id in contig_ids:
                cur_state = True
            else:
                cur_state = False
        if cur_state == True:
            fo.write(line+'\n')		

    fo.close() 
    
    

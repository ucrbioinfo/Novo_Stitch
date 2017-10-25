import sys

def rest_fastas(key_file, old_fasta_file, aligned_contigids_file, lowconf_contigids_file, unaligned_fasta_file, lowconf_fasta_file): 
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


    aligned_contig_ids = set([])
    for line in open(aligned_contigids_file):
        line = line.strip()
        contig_id = line	
        aligned_contig_ids.add(contig_id)


    lowconf_contig_ids = set([])
    for line in open(lowconf_contigids_file):
        line = line.strip()
        contig_id = line
        lowconf_contig_ids.add(contig_id)
    
    

    fo = file(unaligned_fasta_file, 'w')
    cur_state = True
    for line in open(old_fasta_file):
        line = line.strip()
        if line[0] == '>':
            contig_title = line[1:]
            if contig_title not in contig_title_id:
                cur_state = True
            else:
                contig_id = contig_title_id[contig_title]
                if contig_id in aligned_contig_ids:
                    cur_state = False
                else:
                    cur_state = True
        if cur_state == True:
            fo.write(line+'\n')		
    fo.close() 
 
    fo = file(lowconf_fasta_file, 'w')    
    cur_state = True
    for line in open(old_fasta_file):
        line = line.strip()
        if line[0] == '>':
            contig_title = line[1:]
            if contig_title not in contig_title_id:
                cur_state = False
                continue
            contig_id = contig_title_id[contig_title]
            if contig_id in lowconf_contig_ids:
                cur_state = True
            else:
                cur_state = False
        if cur_state == True:
            fo.write(line+'\n')		
    fo.close() 
   
if __name__ == "__main__":
    rest_fastas(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])

 

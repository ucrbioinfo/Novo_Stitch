    
import sys

def chname_fasta_iter(contigs_file, key_file, contigs_chname_file):
    
    contig_name_id = {}
    for line in open(key_file):
        line = line.strip()
        if line[0] == '#' or line[0] == 'C':
            continue
        cols = line.split('\t')
        contig_name = cols[1]
        contig_id = cols[0]
        contig_name_id[contig_name] = contig_id
    
    fo = file(contigs_chname_file, 'w')
    
    contigs_list = []
    cur_new_name = ""
    cur_old_name = ""
    cur_seq = ""
    for line in open(contigs_file):
        line = line.strip()
        if line[0] == '>':
            cur_old_name = line[1:]
            contig_id = contig_name_id[cur_old_name]
            cur_new_name = "N"+contig_id
        else:
            cur_seq = line
            contigs_list.append((cur_old_name, cur_new_name, cur_seq))
    
    for item in contigs_list:
        fo.write(">"+item[1]+" len="+str(len(item[2]))+" "+item[0]+"\n")
        fo.write(item[2]+"\n")

    fo.close() 
    

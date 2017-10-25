    
import sys

def group_contigs(contigs_file, list_file, key_file, output_dir, group_list_file):

    contig_group = {}
    groups = {}
    for line in open(list_file):
        line = line.strip()
        cols = line.split()
        group_id = cols[0]
        contig_id = cols[1]
        if contig_id not in contig_group:
            contig_group[contig_id] = [group_id]
        else:
            contig_group[contig_id].append(group_id)	
        groups[group_id] = 0
    
    grouplist = []
    for group_id in groups:
        grouplist.append(int(group_id))
    grouplist.sort()
    folist = file(group_list_file, 'w')
    for group_id in grouplist:
        folist.write(str(group_id)+"\n")
    
    
    fo_list = {}
    for line in open(contigs_file):
        line = line.strip()
        if line[0] == '>':
            cols = line[1:].split()
            contig_name = cols[0]
            contig_id = contig_name[1:]#eg.N3023
            if contig_id not in contig_group:
                continue
            group_id_list = contig_group[contig_id]
            cur_fo_list = []
            for group_id in group_id_list:
                if group_id in fo_list:	
                    cur_fo = fo_list[group_id]
                    cur_fo_list.append(cur_fo)
                else:
                    fo = file(output_dir+"/molecule_"+group_id+".fasta", 'w')
                    fo_list[group_id] = fo
                    cur_fo = fo
                    cur_fo_list.append(cur_fo)
            for cur_fo in cur_fo_list:
                cur_fo.write(line+'\n')
        else:
            for cur_fo in cur_fo_list:
                cur_fo.write(line+'\n')
                
    
    

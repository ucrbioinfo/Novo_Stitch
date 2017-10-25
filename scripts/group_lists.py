    
import sys

def group_lists(list_file, group_file, output_dir):
    
    group_contig_info = {}
    for line in open(group_file):
        group_id = line.strip()
        group_contig_info[group_id] = []
    
    for line in open(list_file):
        line = line.strip()
        cols = line.split()
        group_id = cols[0]
        group_contig_info[group_id].append(line)
    
    for group_id in group_contig_info:
        contig_info = group_contig_info[group_id]
        file_name = output_dir+"/molecule_"+group_id+".txt"
        fo = file(file_name, 'w')
        for line in contig_info:
            fo.write(line+"\n")
    
    

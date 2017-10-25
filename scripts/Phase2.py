
import os

from group_contigs import group_contigs
from group_lists import group_lists
from step1_best_hit import step1_best_hit
from step2_stitch import step2_stitch

def Phase2_stitch(BLASTN, CPP_1, mtp_fasta_file_chname, list_file, key_file, group_list_file, contigs_dir, listfiles_dir, best_dir, stitched_dir, stitched_fasta_file, stitch_info_file, num_threads, threshold_1, threshold_2):
    if not os.path.isdir(contigs_dir):
        os.makedirs(contigs_dir)
    if not os.path.isdir(listfiles_dir):
        os.makedirs(listfiles_dir)
    if not os.path.isdir(best_dir):
        os.makedirs(best_dir)
    if not os.path.isdir(stitched_dir):
        os.makedirs(stitched_dir)

    group_contigs(mtp_fasta_file_chname, list_file, key_file, contigs_dir, group_list_file) 
    group_lists(list_file, group_list_file, listfiles_dir)
    
    command = "export OMP_NUM_THREADS=" + str(num_threads)
    os.system(command)

    command = CPP_1 + " " + BLASTN + " " + group_list_file + " " + contigs_dir
    os.system(command)

    fout = file(stitched_fasta_file, 'w')
    fstitch = file(stitch_info_file, 'w')
    for line in open(group_list_file):
        print line
        line = line.strip()
        one_contig_file = contigs_dir + "/molecule_" + line + ".fasta"
        one_alm_file = contigs_dir + "/self_molecule_" + line + ".blast"
        one_list_file = listfiles_dir + "/molecule_" + line + ".txt"
        one_best_file = best_dir + "/molecule_" + line + ".blast"
        one_stitched_file = stitched_dir + "/stitched_contig_" + line + ".fasta"

        step1_best_hit(one_alm_file, one_list_file, one_contig_file, one_best_file, threshold_1, threshold_2)
        stitch_info = step2_stitch(one_best_file, one_list_file, one_contig_file, one_stitched_file)
   
        for info in stitch_info:
            fstitch.write(info+"\n")
        
        f = open(one_stitched_file, 'r')
        text = f.read()
        text = text.strip()
        fout.write(text+"\n")
        f.close()
    
    fout.close()
    fstitch.close()
    

import os

from confidence import confidence
from unstitch import unstitch
from step2_stitch import step2_stitch

def Phase3_check(BLASTN, CPP_2, contigs_dir, listfiles_dir, best_dir, stitched_dir, alms_dir, group_list_file, cancel_list_file, filenames_list_file, stitch_info_file, confid_file, adjusted_fasta_file, num_threads, threshold_3, threshold_4, threshold_5):

    if not os.path.isdir(alms_dir):
        os.makedirs(alms_dir)
    
    command = "export OMP_NUM_THREADS=" + str(num_threads)
    os.system(command)

    command = CPP_2 + " " + BLASTN + " " + group_list_file + " " + stitched_dir + " " + contigs_dir + " " + alms_dir
    os.system(command)

    fout = file(confid_file, 'w')
    for line in open(group_list_file):
        print "calculating confidence for molecule", line
        line = line.strip()
        one_contig_file = contigs_dir + "/molecule_" + line + ".fasta"
        one_stitched_file = stitched_dir + "/stitched_contig_" + line + ".fasta"
        one_alms_file = alms_dir + "/alignments_" + line + ".blast"
        confidence(fout, line, one_alms_file, one_contig_file, one_stitched_file, threshold_3, threshold_4)        
        print 
    fout.close()
           
    unstitch(stitch_info_file, confid_file, contigs_dir, stitched_dir, group_list_file, cancel_list_file, filenames_list_file, threshold_5)
    
    print "merge all fasta files"
    fout = file(adjusted_fasta_file, 'w')
    for line in open(filenames_list_file):
        line = line.strip()
        one_stitched_file = stitched_dir + "/" + line
        f = open(one_stitched_file, 'r')
        text = f.read()
        text = text.strip()
        fout.write(text+"\n")
        f.close()
    fout.close()

    
    



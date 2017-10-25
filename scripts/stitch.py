import os.path
from shutil import copyfile

from merge_inputs import merge_inputs
from chname_fasta import chname_fasta
from fasta_long_seqs import fasta_long_seqs
from Phase1 import Phase1_reduction
from Phase2 import Phase2_stitch
from Phase3 import Phase3_check
from post_process import post_process
from iterations_stop import iterations_stop

def output_paras(fasta_list_file, output_dir, optmap, optmap_type, num_threads, REFALIGNER, BLASTN, GLPSOL, threshold_1, threshold_2, threshold_3, threshold_4, threshold_5, threshold_6, threshold_7):
    fo = file(output_dir+"/parameters.log", 'w')
    fo.write("fasta list files:\t"+fasta_list_file+"\n")
    fo.write("output directory:\t"+output_dir+"\n")
    fo.write("optical map:\t"+optmap+"\n")
    fo.write("optical map type:\t"+optmap_type+"\n")
    fo.write("number of threads:\t"+str(num_threads)+"\n")
    fo.write("REFALIGNER tool:\t"+REFALIGNER+"\n")
    fo.write("BLASTN tool:\t"+BLASTN+"\n")
    fo.write("GLPSOL tool:\t"+GLPSOL+"\n")
    fo.write("threshold_1:\t"+str(threshold_1)+"\n")
    fo.write("threshold_2:\t"+str(threshold_2)+"\n")
    fo.write("threshold_3\t"+str(threshold_3)+"\n")
    fo.write("threshold_4\t"+str(threshold_4)+"\n")
    fo.write("threshold_5\t"+str(threshold_5)+"\n")
    fo.write("threshold_6\t"+str(threshold_6)+"\n")
    fo.write("threshold_7\t"+str(threshold_7)+"\n")
    fo.close()

def stitch(paras_file, fasta_list_file, output_dir, REFALIGNER, BLASTN, GLPSOL, optmap, optmap_type, num_threads, threshold_1, threshold_2, threshold_3, threshold_4, threshold_5, threshold_6, threshold_7):

    #tools

    DIGEST = "./tools/fa2cmap_multi.pl"
    CPP_1 = "./tools/Script_mapping_contigs2contigs"
    CPP_2 = "./tools/Script_mapping_contigs2stitched_contigs"

    #intermediate files unrelated to each iteration
    input_merged_file = output_dir + "/input_merged.fasta" 
    input_merged_file_chname = output_dir + "/input_merged_chname.fasta"
    input_merged_file_chname_long = output_dir + "/input_merged_chname_long.fasta"
    input_merged_file_chname_nochi = output_dir + "/input_merged_chname_nochi.fasta"
    sub_dir_prefix = output_dir + "/iteration_"
    last_stitched_file = output_dir + "/last_stitched.fasta"
    last_mtp_stitched_file = output_dir + "/last_mtp_stitched.fasta"
    final_stitched_file = output_dir + "/final_stitched.fasta"
    output_file = output_dir + "/final_stitched.nochi.fasta"
    last_refaligner_dir = output_dir + "/last_refaligner"
    last_refaligner_prefix = "stitched_contigs"
    last_silicomap = last_refaligner_dir + "/" + last_refaligner_prefix + "_" + optmap_type + ".cmap"  
    last_list_file = last_refaligner_dir + "/" + last_refaligner_prefix + "_" + optmap_type + "_BNG_VS_seq_list.txt" 
    last_key_file = last_refaligner_dir + "/" + last_refaligner_prefix + "_" + optmap_type + "_key.txt" 
    
    chimeric_start_dir = output_dir + "/chimeric_start"
    chimeric_end_dir = output_dir + "/chimeric_end"

    #check existence of output_dir
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    #pre-operations
    output_paras(fasta_list_file, output_dir, optmap, optmap_type, num_threads, REFALIGNER, BLASTN, GLPSOL, threshold_1, threshold_2, threshold_3, threshold_4, threshold_5, threshold_6, threshold_7) 
    merge_inputs(output_dir, fasta_list_file, input_merged_file)
    chname_fasta(input_merged_file, input_merged_file_chname)
    fasta_long_seqs(input_merged_file_chname, input_merged_file_chname_long) 

    #iterations
    num_iterations = 1
    old_adjusted_fasta_file = ""
    i = 0
    while True:
        #intermediate files and directories for this iteration
        i += 1
        sub_dir = sub_dir_prefix + str(i)
        refaligner_dir = sub_dir + "/refaligner"
        contigs_dir = sub_dir + "/grouped_contigs"
        listfiles_dir = sub_dir + "/grouped_listfiles"
        best_dir = sub_dir + "/best_hits"
        stitched_dir = sub_dir + "/stitched_contigs"
        alms_dir = sub_dir + "/alignments"

        input_fasta_file = sub_dir + "/input_contigs.fasta"
        mtp_fasta_file = sub_dir + "/input_contigs_mtp.fasta"
        mtp_fasta_file_chname = sub_dir + "/input_contigs_mtp_chname.fasta"
        unaligned_fasta_file = sub_dir + "/input_contigs_unaligned.fasta"
        lowconf_fasta_file = sub_dir + "/input_contigs_lowconf.fasta"
        stitched_fasta_file = sub_dir + "/stitched_contigs_total.fasta"
        stitch_info_file = sub_dir + "/stitch_info.log"
        adjusted_fasta_file = sub_dir + "/adjusted_contigs_total.fasta"

        refaligner_prefix = "input_contigs"
        silicomap = refaligner_dir + "/" + refaligner_prefix + "_" + optmap_type + ".cmap"  
        list_file = refaligner_dir + "/" + refaligner_prefix + "_" + optmap_type + "_BNG_VS_seq_list.txt" 
        key_file = refaligner_dir + "/" + refaligner_prefix + "_" + optmap_type + "_key.txt" 
        group_list_file = sub_dir + "/group_id_list.txt"
        cancel_list_file = sub_dir + "/cancel_id_list.txt"
        filenames_list_file = sub_dir + "/new_filenames.txt"
        confid_file = sub_dir + "/confidences.txt"

        
        #pre-processing
        if not os.path.isdir(sub_dir):
            os.makedirs(sub_dir)
        if i == 1:
            copyfile(input_merged_file_chname_long, input_fasta_file) 
        else:
            src = sub_dir_prefix + str(i-1) + "/adjusted_contigs_total.fasta"
            copyfile(src, input_fasta_file)
        
        #core operations
        Phase1_reduction(REFALIGNER, DIGEST, GLPSOL, sub_dir, refaligner_dir, input_fasta_file, optmap, optmap_type, silicomap, refaligner_prefix, list_file, key_file, mtp_fasta_file, mtp_fasta_file_chname, unaligned_fasta_file, lowconf_fasta_file, num_threads, threshold_6, threshold_7)
         
        Phase2_stitch(BLASTN, CPP_1, mtp_fasta_file_chname, list_file, key_file, group_list_file, contigs_dir, listfiles_dir, best_dir, stitched_dir, stitched_fasta_file, stitch_info_file, num_threads, threshold_1, threshold_2)
        
        Phase3_check(BLASTN, CPP_2, contigs_dir, listfiles_dir, best_dir, stitched_dir, alms_dir, group_list_file, cancel_list_file, filenames_list_file, stitch_info_file, confid_file, adjusted_fasta_file, num_threads, threshold_3, threshold_4, threshold_5)
         
        if num_iterations != 1:
            if iterations_stop(old_adjusted_fasta_file, adjusted_fasta_file) == True:
                break
        old_adjusted_fasta_file = adjusted_fasta_file
        num_iterations += 1

    #post-operations
    src = sub_dir_prefix + str(num_iterations) + "/adjusted_contigs_total.fasta"
    copyfile(src, last_stitched_file)
    post_process(REFALIGNER, DIGEST, GLPSOL, sub_dir, last_refaligner_dir, last_stitched_file, optmap, optmap_type, last_silicomap, last_refaligner_prefix, last_list_file, last_key_file, last_mtp_stitched_file, final_stitched_file, num_threads, threshold_6, threshold_7)


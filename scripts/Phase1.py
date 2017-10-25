import os.path
from shutil import copyfile

from mtp import mtp
from mtp_fasta import mtp_fasta
from rest_fastas import rest_fastas
from chname_fasta_iter import chname_fasta_iter

def Phase1_reduction(REFALIGNER, DIGEST, GLPSOL, sub_dir, refaligner_dir, input_fasta_file, optmap, optmap_type, silicomap, refaligner_prefix, list_file, key_file, mtp_fasta_file, mtp_fasta_file_chname, unaligned_fasta_file, lowconf_fasta_file, num_threads, threshold_6, threshold_7):
    fasta_file_in_refaligner = refaligner_dir + "/" + refaligner_prefix + ".fasta"
    out = refaligner_dir + "/" + refaligner_prefix + "_" + optmap_type + "_BNG_VS_seq"
    out2 = refaligner_dir + "/" + refaligner_prefix + "_" + optmap_type
    aligned_contigids_file = out + "_aligned.txt"
    lowconf_contigids_file = out + "_lowconf.txt"

    if not os.path.isdir(refaligner_dir):
        os.makedirs(refaligner_dir)
    copyfile(input_fasta_file, fasta_file_in_refaligner)

    command = "perl " + DIGEST + " -v -i " + fasta_file_in_refaligner + " -e " + optmap_type + " -m 0 -M 0"
    os.system(command)
    command = REFALIGNER + " -f -ref " + optmap + " -i " + silicomap + " -o " + out + " -endoutlier 1e-3 -outlier 1e-4 -extend 1 -FN 0.05 -FP 0.5 -sf 0.2 -sd 0.1 -sr 0.02 -res 2.9 -resSD 0.75 -mres 1.2 -A 5 -biaswt 0 -M 3 -Mfast 0 -maxmem 256 -maxthreads " + str(num_threads) + " -deltaX 9 -deltaY 9 -xmapchim 14 -RepeatMask 2 0.01 -RepeatRec 0.7 0.6 -T 1e-15 -stdout -stderr -xmaplen -indel" 
    os.system(command)
       
    mtp(out, out2, sub_dir, GLPSOL, threshold_7, threshold_6) 

    mtp_fasta(list_file, key_file, input_fasta_file, mtp_fasta_file)

    rest_fastas(key_file, input_fasta_file, aligned_contigids_file, lowconf_contigids_file, unaligned_fasta_file, lowconf_fasta_file)

    chname_fasta_iter(mtp_fasta_file, key_file, mtp_fasta_file_chname)    


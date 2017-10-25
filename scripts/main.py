#!/usr/bin/python

import sys, getopt
import os.path

from stitch import stitch

def paras_in_file(paras_file):
    if not os.path.isfile(paras_file):    
        print "ERROR! Parameter file " + paras_file + " doesn't exist!"
        exit()
    argv = []
    for line in open(paras_file):
        line = line.strip()
        cols = line.split()
        argv += cols
    return argv   


def main():

    #default input and output
    paras_file = "./paras.txt"
    fasta_list_file = "./input_files_list.txt"
    output_dir = "./output_dir"
    optmap = "/home/stelo/BIONANO_in_progress/vu_162_180K.cmap"
    optmap_type = "BspQI" 
   
    #default tools 
    REFALIGNER = "RefAligner" 
    BLASTN = "blastn"
    GLPSOL = "glpsol" 

    #default parameters
    num_threads = 32
    threshold_1 = 3000#1000#3000
    threshold_2 = 0.1#0.2#0.1
    threshold_3 = 10000#5000#10000
    threshold_4 = 0.5
    threshold_5 = 0.9
    threshold_6 = 25
    threshold_7 = 0.2

    #obtaining parameters
    paras_string = "f:i:o:m:t:p:x:y:g:a:b:c:d:e:h:r:"
    opts, args = getopt.getopt(sys.argv[1:], paras_string)
    for op, value in opts:
        if op == "-f":
            paras_file = value
            argv_in_file = paras_in_file(paras_file)
            opts_in_file, args_in_file = getopt.getopt(argv_in_file, paras_string) 
            for fop, fvalue in opts_in_file:
                if fop == "-i":
                    fasta_list_file = fvalue
                elif fop == "-o":
                    output_dir = fvalue
                elif fop == "-m":
                    optmap = fvalue
                elif fop == "-t":
                    optmap_type = fvalue
                elif fop == "-p":
                    num_threads = int(fvalue)
                elif fop == "-x":
                    REFALIGNER = fvalue
                elif fop == "-y":
                    BLASTN = fvalue
                elif fop == "-g":
                    GLPSOL = fvalue
                elif fop == "-a":
                    threshold_1 = int(fvalue)
                elif fop == "-b":
                    threshold_2 = float(fvalue)
                elif fop == "-c":
                    threshold_3 = int(fvalue)
                elif fop == "-d":
                    threshold_4 = float(fvalue)
                elif fop == "-e":
                    threshold_5 = float(fvalue)
                elif fop == "-h":
                    threshold_6 = float(fvalue)
                elif fop == "-r":
                    threshold_7 = float(fvalue)

        elif op == "-i":
            fasta_list_file = value
        elif op == "-o":
            output_dir = value
        elif op == "-m":
            optmap = value
        elif op == "-t":
            optmap_type = value
        elif op == "-p":
            num_threads = int(value)
        elif op == "-x":
            REFALIGNER = value
        elif op == "-g":
            GLPSOL = value
        elif op == "-y":
            BLASTN = value
        elif op == "-a":
            threshold_1 = int(value)
        elif op == "-b":
            threshold_2 = float(value)
        elif op == "-c":
            threshold_3 = int(value)
        elif op == "-d":
            threshold_4 = float(value)
        elif op == "-e":
            threshold_5 = float(value)
        elif op == "-h":
            threshold_6 = float(value)
        elif op == "-r":
            threshold_7 = float(value)
    stitch(paras_file, fasta_list_file, output_dir, REFALIGNER, BLASTN, GLPSOL, optmap, optmap_type, num_threads, threshold_1, threshold_2, threshold_3, threshold_4, threshold_5, threshold_6, threshold_7)        

if __name__ == "__main__":
    main()

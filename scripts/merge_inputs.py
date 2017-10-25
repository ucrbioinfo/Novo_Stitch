import sys
import os.path

def output_input(filename, fname_list):
    fo = file(filename, 'w')
    for fname in fname_list:
        fo.write(fname+"\n")
    fo.close()

def merge_inputs(output_dir, fasta_list_file, merged_fasta_file):
    if not os.path.isfile(fasta_list_file):    
        print "ERROR! Inputs list file " + fasta_list_file + " doesn't exist!"
        exit()
    fname_list = set([])
    for line in open(fasta_list_file):
        line = line.strip()
        if line == "":
            continue
        fname_list.add(line)
    output_input(output_dir+"/input.log", fname_list)
    fout = file(merged_fasta_file, 'w')
    for fn in fname_list:
        print "Reading file: ", fn
        f = open(fn, 'r')
        text = f.read()
        text = text.strip()
        fout.write(text+"\n")
        f.close()
    fout.close()

if __name__ == "__main__":
    merge_inputs(output_dir, sys.argv[1], sys.argv[2])




import sys

def chname_fasta(input_file, output_file):
    fout = file(output_file, 'w')
    seq_num = 1
    for line in open(input_file):
        line = line.strip()
        if line[0] == '>':
            fout.write(">SEQ_"+str(seq_num)+"\n")
            seq_num += 1
        else:
            fout.write(line+"\n")
    fout.close()


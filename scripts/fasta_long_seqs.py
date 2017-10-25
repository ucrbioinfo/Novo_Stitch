import sys

def fasta_long_seqs(fasta_file, new_fasta_file):
    # read fasta file
    seqs_list = []
    names_list = []
    for line in open(fasta_file):
        line = line.strip()
        if line[0] == '>':
            seq_name = line[1:]
            seqs_list.append([])
            names_list.append(seq_name)
        else:
            seqs_list[len(seqs_list)-1].append(line)
    

    long_seqs_list = []
    for i in range(0, len(seqs_list)):
        seq_list = seqs_list[i]
        whole_seq = ''.join(seq_list)
        long_seqs_list.append(whole_seq)

    fo = file(new_fasta_file, 'w')
    for i in range(0, len(names_list)):
        fo.write(">"+names_list[i]+"\n")
        fo.write(long_seqs_list[i]+"\n")
    fo.close()


if __name__ == "__main__":
    long_seqs(sys.argv[1], sys.argv[2])





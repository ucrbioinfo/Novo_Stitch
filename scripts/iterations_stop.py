
def iterations_stop(old_adjusted_fasta_file, adjusted_fasta_file):
    num_old = 0
    for line in open(old_adjusted_fasta_file):
        line = line.strip()
        if line[0] == '>':
            num_old += 1
    num_new = 0
    for line in open(adjusted_fasta_file):
        line = line.strip()
        if line[0] == '>':
            num_new += 1

    if num_new == num_old:
        return True
    elif num_new < num_old:
        return False
    else:
        print "ERROR! The number of contigs cannot increase!!!"
        exit()


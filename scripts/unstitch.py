    
def unstitch(stitch_info_file, confidence_file, contigs_dir, stitched_dir, group_id_file, cancel_id_file, new_filename_list, conf_cutoff): 
     
    # get stitch information
    stitch_info = {}
    for line in open(stitch_info_file):
        line = line.strip()
        cols = line.split('\t')
        molecule_id = cols[0]
        stitched_name = cols[1]
        contig_list = cols[2]
        list_cols = contig_list.split()
        stitch_info[molecule_id, stitched_name] = list_cols 

    # get cancel_set
    cancel_set = {}
    for line in open(confidence_file):
        line = line.strip()
        cols = line.split()
        mole_id = cols[0]
        contig = cols[1]
        stitched_contig = cols[2]
        conf = float(cols[3])
        if conf < conf_cutoff:#0.9:
            if mole_id not in cancel_set:
                cancel_set[mole_id] = set([])
            cancel_set[mole_id].add(stitched_contig)
            print "stitched contig", stitched_contig, "in molecule", mole_id, "is cancalled"   
    
    fo_cancel = file(cancel_id_file, 'w')
    for mole_id in cancel_set:
        fo_cancel.write(mole_id+"\n")
    fo_cancel.close()

    for mole_id in cancel_set:
        # generate a new fasta file to replace old fasta file
        unstitched_file = contigs_dir + "/molecule_" + mole_id + ".fasta"
        old_stitched_file = stitched_dir + "/stitched_contig_" + mole_id + ".fasta"
        new_stitched_file = stitched_dir + "/stitched_contig_" + mole_id + "_new.fasta" 
        old_stitched_seqs = {}
        current_name = ""
        for line in open(old_stitched_file):
            line = line.strip()
            if line[0] == '>':
                cols = line[1:].split()
                current_name = cols[0]
                old_stitched_seqs[current_name] = ""
            else:
                old_stitched_seqs[current_name] += line
        unstitched_seqs = {}
        current_name = ""
        for line in open(unstitched_file):
            line = line.strip()
            if line[0] == '>':
                cols = line[1:].split()
                current_name = cols[0]
                unstitched_seqs[current_name] = ""
            else:
                unstitched_seqs[current_name] += line
        new_stitched_seqs = {}
        for stitched_contig in old_stitched_seqs:
            if stitched_contig not in cancel_set[mole_id]:
                new_stitched_seqs[stitched_contig] = old_stitched_seqs[stitched_contig]
            else:
                contigs = stitch_info[mole_id, stitched_contig]
                for c in contigs:
                    new_stitched_seqs[c] = unstitched_seqs[c]
        fnew = file(new_stitched_file, 'w')
        for contig_name in new_stitched_seqs:
            seq = new_stitched_seqs[contig_name]
            title = ">"+contig_name+" len="+str(len(seq))
            fnew.write(title+"\n")
            fnew.write(seq+"\n")
        fnew.close()
        print "new stitched fasta file", new_stitched_file, "is generated"

    # keep name of fasta files 
    f_fn = file(new_filename_list, 'w')
    for line in open(group_id_file):
        line = line.strip()
        if line in cancel_set:
            f_fn.write("stitched_contig_" + line + "_new.fasta\n")
        else:
            f_fn.write("stitched_contig_" + line + ".fasta\n")
    
    
    
    
    
    
    

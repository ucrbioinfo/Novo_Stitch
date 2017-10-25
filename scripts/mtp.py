#!/usr/bin/python
import csv
import sys

from sets import Set

from collections import defaultdict

import collections
import math

from alignment import Alignment 
from get_mst import get_mst
from merge_DAGs import merge_DAGs
from subDAG import get_subDAGs
from unify_coords import unify_coords
from false_alms import false_alms



def copy_alms(old_alms, removed):
    current_alms = {}
    for ref in old_alms:
        current_alms[ref] = {}
        for qry in old_alms[ref]:
            x = old_alms[ref][qry]
            if removed[ref, qry] == False:
                current_alms[ref][qry] = x
    return current_alms

def print_alms(alms):
    for ref in alms:
        for qry in alms[ref]:
            x = alms[ref][qry]
            print "alignment:", (ref, qry), x.start, x.end

def count_alms(alms):
    num_alms = 0
    for ref in alms:
        num_alms += len(alms[ref])
    return num_alms

def output_alms(alms, filename):
    fo = file(filename, 'w')
    fo.write("ref_id\tqry_id\torientation\tconfidence\tstart\tend\n")
    for ref in alms:
        for qry in alms[ref]:
            x = alms[ref][qry]
            fo.write(str(ref)+"\t"+str(qry)+"\t"+str(x.orientation)+"\t"+str(x.confidence)+"\t"+str(x.start)+"\t"+str(x.end)+"\n")
    fo.close()

def different_contigs(alms_1, alms_2):
    contigs_1 = set([])
    for ref in alms_1:
        for qry in alms_1[ref]:
            contigs_1.add(qry)
    contigs_2 = set([])
    for ref in alms_2:
        for qry in alms_2[ref]:
            contigs_2.add(qry)
    diff = set([])
    for c in contigs_1:
        if c not in contigs_2:
            diff.add(c)
    return diff

def output_contigs(contigs, filename):
    fo = file(filename, 'w')
    for c in contigs:
        fo.write(str(c)+"\n")
    fo.close()

def output_forest(forest, vertex_orientations, filename):
    fo = file(filename, 'w')
    for root in forest:
        fo.write("root:"+str(root)+"\n")
        oppo_tree = {}
        for v in forest[root]:
            for u in forest[root][v]:
                oppo_tree[u] = v
        oppo_tree[root] = -1

        for v in forest[root]:
            fo.write("node:"+str(v)+"\t")
            fo.write(str(vertex_orientations[root][v])+"\t")
            fo.write("parent:"+str(oppo_tree[v])+"\t")
            fo.write("children: ")
            for u in forest[root][v]:
                fo.write(str(u)+" ")
            fo.write("\n")
        fo.write("\n")
    fo.close()

def output_DAGs(DAGs, filename):
    fo = file(filename, 'w')
    for i in range(0, len(DAGs)):
        DAG = DAGs[i]
        oppo_DAG = {}
        for v in DAG:
            oppo_DAG[v] = set([])
        for v in DAG:
            for u in DAG[v]:
                oppo_DAG[u].add(v)
        fo.write("source: ")
        for v in oppo_DAG:
            if oppo_DAG[v] == set([]):
                fo.write(str(v)+" ")
        fo.write("\n")
        fo.write("sink: ")
        for v in DAG:
            if DAG[v] == set([]):
                fo.write(str(v)+" ")
        fo.write("\n")
        
        for v in DAG:
            fo.write("node:"+str(v)+"\t")
            fo.write("incoming: ")
            for u in oppo_DAG[v]:
                fo.write(str(u)+" ")
            fo.write("\t")
            fo.write("outgoing: ")
            for u in DAG[v]:
                fo.write(str(u)+" ")
            fo.write("\n")
        fo.write("\n")
    fo.close()

def mtp(myfile, myfile2, output_dir, GLPSOL, false_alm_threshold, min_confidence):
    
    # discard alignments below min_confidence
    #min_confidence = 25
    header_lines = 10
    header = []
    # alignment overhangs above this number of bps are considered chimeric
    #minrefoverhang = 100000
    #minqryoverhang = 100000
	
    all_alms = {} # stores all the Alignments for all groups, all_groups[ref] should contain molecule ref
    qualify_alms = {} # only keep one alignment(the one with highest confidence) for each contig in one molecule
    removed = {} # removed[ref,qry] == True means alignment for (ref, qry) is already removed

    # collecting alignments and store in all_groups
    print '---------------read .xmap file-------------------'
    with open(myfile+'.xmap', 'rb') as csvfile:
        csvreader = csv.reader(csvfile, delimiter='\t')
        for i in range(header_lines): # 10 lines of header
            header.append(csvreader.next()) # save them
        # read the first non-header line
        while True:
            try:
                row = csvreader.next()
                x = Alignment(int(row[1]),int(row[2]),float(row[3]),float(row[4]),float(row[5]),
                                  float(row[6]),row[7],float(row[8]),row[9],float(row[10]),
                                  float(row[11]),int(row[12]),row[13])
                if x.ref not in all_alms:
                    all_alms[x.ref] = [x]
                else:
                    all_alms[x.ref].append(x)
            except StopIteration:
                break
    num_all_alms = 0
    for ref in all_alms:
        #print 'collected', len(all_alms[ref]), 'alignments for molecule', ref
        num_all_alms += len(all_alms[ref])
    print "In total, the number of alignments collected is ", num_all_alms
    
    
    # only keep one alignment(the one with highest confidence) for each contig in one molecule    
    for ref in all_alms:
        group = all_alms[ref]
        qry_bestx = {}
        for x in group:
            if x.qry not in qry_bestx:
                qry_bestx[x.qry] = x
            else:
                if x.confidence > qry_bestx[x.qry].confidence:
                    qry_bestx[x.qry] = x

        qualify_alms[ref] = {}
        for qry in qry_bestx:
            qualify_alms[ref][qry] = qry_bestx[qry]

    num_qualify_alms = 0
    for ref in qualify_alms:
        num_qualify_alms += len(qualify_alms[ref])
    print "In total, the number of alignments in qualify_alms is ", num_qualify_alms

    # initialize removed array
    for ref in qualify_alms:
        for qry in qualify_alms[ref]:
            removed[ref,qry] = False

    # find the reference-based coordinates for each alignments
    for ref in qualify_alms:
        for qry in qualify_alms[ref]:
            x = qualify_alms[ref][qry]
            if (x.orientation == '+'):
                x.qry_left_overlen = x.qrystartpos
                x.qry_right_overlen = x.qrylen - x.qryendpos
            else:
                x.qry_left_overlen = x.qrylen - x.qrystartpos
                x.qry_right_overlen = x.qryendpos
            x.start = x.refstartpos - x.qry_left_overlen
            x.end = x.refendpos + x.qry_right_overlen 

    current_alms = copy_alms(qualify_alms, removed)
    output_alms(current_alms, output_dir+"/alms_0_initial.log")
    print "Initially, the number of alignments is", count_alms(current_alms)    
    alms_0 = copy_alms(qualify_alms, removed)
    aligned_contigs = different_contigs(alms_0, {})
    output_contigs(aligned_contigs, myfile+'_aligned.txt')
    print '---------------END-------------------'


    # remove low confidence alignments
    print '---------------Remove low quality alignments---------------'
    for ref in qualify_alms:
        for qry in qualify_alms[ref]:
            x = qualify_alms[ref][qry]
            if x.confidence < min_confidence:
                removed[ref, qry] = True
                print 'alignment (', ref, ',', qry, ') is low quality and removed'
    num_alms = 0
    for ref in qualify_alms:
        for qry in qualify_alms[ref]:
            if removed[ref, qry] == False:
                num_alms += 1

    current_alms = copy_alms(qualify_alms, removed)
    output_alms(current_alms, output_dir+"/alms_1_removed_lowconf.log")
    print "After removing low confidence alignments, the number of alignments is", count_alms(current_alms)    
    alms_1 = copy_alms(qualify_alms, removed)
    lowconf_contigs = different_contigs(alms_0, alms_1)
    output_contigs(lowconf_contigs, myfile+'_lowconf.txt')
    print '---------------End---------------'

   
    print '---------------removing false positive alignments-------------------'
    current_alms = copy_alms(qualify_alms, removed)
    false_alms(GLPSOL, false_alm_threshold, current_alms, removed, output_dir)   
    current_alms = copy_alms(qualify_alms, removed)
    output_alms(current_alms, output_dir+"/alms_2_removed_false_alms.log")
    print "After removing false positive alignments, the number of alignments is", count_alms(current_alms) 
    print '---------------END-------------------'
       
    print '---------------removing contained contigs locally-------------------'
    for ref in qualify_alms:
        for q1 in qualify_alms[ref]:
            x = qualify_alms[ref][q1]
            for q2 in qualify_alms[ref]:
                if q2 <= q1:
                    continue
                y = qualify_alms[ref][q2]
                if (x.start >= y.start) and (x.end <= y.end):
                    removed[ref, q1] = True
                    print [ref, q1], "alignment is removed becasue it's contained in alignment", [ref, q2]
                elif (y.start >= x.start) and (y.end <= x.end):
                    removed[ref, q2] = True
                    print [ref, q2], "alignment is removed becasue it's contained in alignment", [ref, q1]
    current_alms = copy_alms(qualify_alms, removed)
    output_alms(current_alms, output_dir+"/alms_3_removed_contained_locally.log")
    print "After removing contained alignments locally, the number of alignments is", count_alms(current_alms)
    print '---------------END-------------------'


    #build the mst
    print '---------------building the mst-------------------'
    fo = file(output_dir+"/ugraph_1.log", 'w')
    current_alms = copy_alms(qualify_alms, removed) 
    forest, vertex_orientations = get_mst(current_alms, fo)
    fo.close()
    output_forest(forest, vertex_orientations, output_dir+"/forest_1.log")
    print '---------------END-------------------'
    # unify the coordinates
    print '---------------unifying the coordinates-------------------'
    current_alms = copy_alms(qualify_alms, removed) 
    unify_alms = unify_coords(output_dir, current_alms, forest, vertex_orientations)    

    removed_unify = {}
    for root in unify_alms:
        for qry in unify_alms[root]:
            removed_unify[root, qry] = False

    contigs = set([])
    for root in unify_alms:
        for qry in unify_alms[root]:
            if qry in contigs:
                print qry, "appears in more than 1 trees"
            contigs.add(qry)
    current_alms = copy_alms(unify_alms, removed_unify)
    output_alms(current_alms, output_dir+"/alms_4_unified.log")
    print "After unifying the coordinates, the number of alignments is", count_alms(current_alms)
    print '---------------END-------------------'

    print '---------------removing contained contigs globally-------------------'
    contained = set([])
    for root in unify_alms:
        for q1 in unify_alms[root]:
            x = unify_alms[root][q1]
            for q2 in unify_alms[root]:
                if q2 <= q1:
                    continue
                y = unify_alms[root][q2]
                if (q2 not in contained) and (x.start >= y.start) and (x.end <= y.end):
                    contained.add(q1)
                    removed_unify[root, q1] = True
                    print [root, q1], "alignment is removed becasue it's contained in alignment", [root, q2]
                elif (q1 not in contained) and (y.start >= x.start) and (y.end <= x.end):
                    contained.add(q2)
                    removed_unify[root, q2] = True
                    print [root, q2], "alignment is removed becasue it's contained in alignment", [root, q1]
    for root in unify_alms:
        for qry in unify_alms[root]:
            if qry in contained and removed_unify[root, qry] == False:
                removed_unify[root, qry] = True
                print [root, qry], "alignment is removed because qry is contained contig"

    current_alms = copy_alms(unify_alms, removed_unify)
    output_alms(current_alms, output_dir+"/alms_5_removed_contained_globally.log")
    print "After removing contained contigs globally, the number of alignments is", count_alms(current_alms)
    print '---------------END-------------------'

    #build new mst
    print '---------------building new mst-------------------'
    fo = file(output_dir+"/ugraph_2.log", 'w')
    current_alms = copy_alms(unify_alms, removed_unify) 
    forest_unify, vertex_orientations_unify = get_mst(current_alms, fo)
    fo.close()
    output_forest(forest_unify, vertex_orientations_unify, output_dir+"/forest_2.log")
    print '---------------END-------------------'

    print '---------------merging DAGs-------------------'
    current_alms = copy_alms(unify_alms, removed_unify)
    DAGs = merge_DAGs(current_alms, forest_unify, vertex_orientations_unify)
    output_DAGs(DAGs, output_dir+"/dags.log")
    print '---------------END-------------------'

    #DAG to mtp contig set
    print '---------------mtp-------------------'
    mtp_node_set = get_subDAGs(DAGs, output_dir)
    current_alms = copy_alms(unify_alms, removed_unify)
    mtp = []
    for ref in current_alms:
        for qry in current_alms[ref]:
            x = current_alms[ref][qry]
            if qry in mtp_node_set:
                mtp.append(x)
            else:
                removed_unify[ref, qry] = True
    mtp.sort(key=lambda x: (x.ref, x.start))
    print "In total, the number of alignments in mtp is", len(mtp)

    current_alms = copy_alms(unify_alms, removed_unify)
    output_alms(current_alms, output_dir+"/alms_6_mtp.log")    
    print '---------------END-------------------'

    print '---------------scaling-------------------'
    # calculating scaling
    qry_len = {}
    with open(myfile2+'_key.txt') as f_key:
        for i in range(0, 4): # 4 header lines
            f_key.readline()  
        for line in f_key:
            line = line.strip()
            cols = line.split('\t')
            qry_id = int(cols[0])
            seq_len = int(cols[2])
            qry_len[qry_id] = seq_len
    scaling = 0
    num = 0
    with open(myfile+'_q.cmap') as f_q:
        for i in range(0, 11): # 11 header lines
            f_q.readline()
        for line in f_q:
            line = line.strip()
            cols = line.split('\t')
            qry_id = int(cols[0])
            appr_len = float(cols[1])
            seq_len = qry_len[qry_id]
            scaling += appr_len/seq_len
            num += 1
    scaling /= num # scaling=1.02258059775


    print '---------------outputing-------------------'
    # save the MTP in a new xmap file and count the number of unitigs in each assembly
    with open(myfile+'_list.txt', 'wb') as listfile:
        with open(myfile+'_mtp.xmap', 'wb') as csvfile:
            csvwriter = csv.writer(csvfile, delimiter='\t')
            # copies the old xmap header
            for x in header: 
                csvwriter.writerow(x)
            i = 1 # progressive number
	    # for steve ->
#	    scaling = 1.02257561752017878915 # scaling fact from opt map to BP
            previous = 0 # previous qry contig, to remove dups
            for x in mtp:
                # save the contig in listfile only if it is a new one
                if (x.qry != previous):
                    #listfile.write(str(x.qry)+'\n')
                    previous = x.qry
		# for steve ->
		listfile.write(str(x.ref)+'\t'+str(x.qry)+'\t'+str(int(round(float(x.start)/scaling)))+'\t'+str(int(round(float(x.end)/scaling)))+'\t'+x.orientation+'\n')
                # dump the alignment
                csvwriter.writerow([i]+x.unpack()) 
                i += 1
    csvfile.close()
    listfile.close()


if __name__ == "__main__":
    mtp(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], float(sys.argv[5]), float(sys.argv[6]))


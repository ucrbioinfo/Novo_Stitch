#!/usr/bin/python
import csv
import sys

from sets import Set

from collections import defaultdict

import collections
import math
import copy
from alignment import Alignment



def change_coords(fo, alms, onetree_unify_alms, onetree, onetree_v_oris, real_root, root):
    if root == real_root:
    #for real root vertex, just copy coordiantes
        fo.write("root:"+str(real_root)+"\t"+"node:"+str(root)+"\t"+str(1)+"\tave_shift:"+str(0)+"\n")
        for qry in alms[real_root]:
            onetree_unify_alms[qry] = copy.deepcopy(alms[real_root][qry])
            y = alms[real_root][qry]
            fo.write("newqry:"+str(qry)+"\told: ("+str(y.orientation)+", "+str(y.start)+", "+str(y.end)+")\tnew: ("+str(y.orientation)+", "+str(y.start)+", "+str(y.end)+")\n")
        fo.write("\n")
    else:
    #for other vertices, change coordinates of contigs aligned to molecule root
        # use common contigs to get average shift
        ref_ori = onetree_v_oris[root]
        shifts = []
        for qry in alms[root]:
            if qry in onetree_unify_alms: # common contigs of vertex root and the set of visited vertices 
                x = onetree_unify_alms[qry]
                x_middle = float(x.start + x.end)/2
                y = alms[root][qry]
                y_middle = float(y.start + y.end)/2
                if ref_ori == 1:
                    shift = x_middle - y_middle
                elif ref_ori == -1:
                    shift = x_middle + y_middle
                else:
                    print "ERROR: orientation of this molecule is not 1 or -1!!!"
                    exit()
                shifts.append(shift)
                fo.write("root:"+str(real_root)+"\t"+"node:"+str(root)+"\t"+str(ref_ori)+"\t"+"qry:"+str(qry)+"\ttree: ("+str(x.start)+", "+str(x.end)+", "+str(x_middle)+")\tnode: ("+str(y.start)+", "+str(y.end)+", "+str(y_middle)+")\tshift:"+str(shift)+"\n")
        ave_shift = sum(shifts)/len(shifts)
        fo.write("root:"+str(real_root)+"\t"+"node:"+str(root)+"\t"+str(ref_ori)+"\tave_shift:"+str(ave_shift)+"\n")

        # use average shift to change coordinates for other non-common contigs
        for qry in alms[root]:
            if qry in onetree_unify_alms: # pick non-common contigs
                continue
            y = alms[root][qry]
            if y.orientation == '+':
                qry_ori = 1
            elif y.orientation == '-':
                qry_ori = -1
            else:
                print "ERROR: orientation of this contig is not avalable!!!"
                exit()
            qry_ori = qry_ori * ref_ori
            if qry_ori == 1:
                qry_ori = '+'
            elif qry_ori == -1:
                qry_ori = '-'
            else:
                print "ERROR: orientation of this contig is not avalable!!!!!!"
                exit()
            if ref_ori == 1:
                x_start = min(y.start + ave_shift, y.end + ave_shift)
                x_end = max(y.start + ave_shift, y.end + ave_shift)
            elif ref_ori == -1:
                x_start = min(ave_shift - y.start, ave_shift - y.end)
                x_end = max(ave_shift - y.start, ave_shift - y.end)
            else:
                print "ERROR: orientation of this contig is not avalable!!!!!!"
                exit()

            x = copy.deepcopy(y)
            x.ref = real_root
            x.start = x_start
            x.end = x_end
            x.orientation = qry_ori
            onetree_unify_alms[qry] = x
            fo.write("newqry:"+str(qry)+"\told: ("+str(y.orientation)+", "+str(y.start)+", "+str(y.end)+")\tnew: ("+str(x.orientation)+", "+str(x.start)+", "+str(x.end)+")\n")
    fo.write("\n")
    #search children
    children = onetree[root] 
    for c in children:
        change_coords(fo, alms, onetree_unify_alms, onetree, onetree_v_oris, real_root, c)

def unify_coords(output_dir, alms, forest, vertex_orientations):
    #change coordinates
    fo = file(output_dir+"/shifts.log", 'w')
    unify_alms = {}
    for root in forest:
        onetree = forest[root]
        onetree_v_oris = vertex_orientations[root]
        onetree_unify_alms = {}
        change_coords(fo, alms, onetree_unify_alms, onetree, onetree_v_oris, root, root) 
        unify_alms[root] = onetree_unify_alms
    fo.close()
    return unify_alms



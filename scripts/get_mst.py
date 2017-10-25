#!/usr/bin/python
import csv
import sys
from sets import Set
from collections import defaultdict
import collections
import math

from kruskal import Graph 

def edge_ornot(alms, r1, r2):
    common = set([])
    for qry in alms[r1]:
        if qry in alms[r2]:
            common.add(qry)
    if common == set([]):
        return False
    else:
        return True

def calculate_weight(alms, r1, r2):
    common = set([])
    for qry in alms[r1]:
        if qry in alms[r2]:
            common.add(qry)
    total_conf = 0
    for qry in common:
        total_conf += alms[r1][qry].confidence
        total_conf += alms[r2][qry].confidence
    return 1.0/float(total_conf)

def mst_traversal(mst_ugraph, visited, root, onetree, oneugraph):
    visited[root] = True
    oneugraph[root] = mst_ugraph[root]
    onetree[root] = set([])
    children = mst_ugraph[root]
    for c in children:
        if visited[c] == True:
            continue
        onetree[root].add(c)
        mst_traversal(mst_ugraph, visited, c, onetree, oneugraph)

def tree_traversal_orientation(edge_orientations, onetree_v_oris, root, onetree):
    children = onetree[root]
    for c in children:
        min_one = min(root, c)
        max_one = max(root, c) 
        onetree_v_oris[c] = onetree_v_oris[root] * edge_orientations[min_one, max_one]#e.g. -1 * 1 = -1
        tree_traversal_orientation(edge_orientations, onetree_v_oris, c, onetree)

def get_edge_orientation(alms, r1, r2):
    common = set([])
    for qry in alms[r1]:
        if qry in alms[r2]:
            common.add(qry)
    ori_same = True
    ori_oppo = True
    for qry in common:
        x1 = alms[r1][qry]
        x2 = alms[r2][qry]
        if x1.orientation != x2.orientation:
            ori_same = False
        if x1.orientation == x2.orientation:
            ori_oppo = False
    if ori_same == True:
        return 1
    if ori_oppo == True:
        return -1
    return 0 #the orientation of [r1, r2] cannot be decided becasue of data conflict

def get_mst(alms, fo):
    edges = []
    removed_edges = []
    edge_orientations = {}
    for r1 in alms:
        for r2 in alms:
            if r1 >= r2:
                continue
            if edge_ornot(alms, r1, r2) == False: # there's no edge between r1 and r2
                continue
            ori = get_edge_orientation(alms, r1, r2)
            if ori == 0: # the orientation between r1 and r2 cannot be decided becasue of data conflict
                removed_edges.append([r1, r2])
                print [r1, r2], "the orientation between r1 and r2 cannot be decided and thus the edge is removed"
                continue
            weight = calculate_weight(alms, r1, r2)
            edges.append([r1, r2, float(weight)])
            edge_orientations[r1, r2] = ori
            print "edge", [r1, r2, float(weight), ori], "is added to undirected graph"
    print "In total, the number of vertices in undirected graph is", len(alms)
    print "In total, the number of edges in undirected graph is", len(edges)
    print 
    fo.write("vertices:\n")
    for ref in alms:
        fo.write("node: "+str(ref)+"\n")
    fo.write("\n")
    fo.write("edges:\n")
    for i in range(0, len(edges)):
        [r1, r2, weight] = edges[i]
        ori = edge_orientations[r1, r2]
        fo.write("edge: ("+str(r1)+", "+str(r2)+", "+str(weight)+", "+str(ori)+")\n")
    fo.write("\n")
    fo.write("removed edges:\n")
    for i in range(0, len(removed_edges)):
        [r1, r2] = removed_edges[i]
        fo.write("removed_edge: ("+str(r1)+", "+str(r2)+")\n")
    fo.write("\n")
    #build undirected graph
    node_ref = []
    ref_node = {}
    for ref in alms:
        node_ref.append(ref)
        ref_node[ref] = len(node_ref) - 1

    g = Graph(len(node_ref))
    for item in edges:
        r1 = item[0]
        r2 = item[1]
        weight = item[2]
        node1 = ref_node[r1]
        node2 = ref_node[r2]
        g.addEdge(node1, node2, weight)
    t = g.KruskalMST()

    #build tree
    mst_edges = {}
    for item in t:
        node1 = item[0]
        node2 = item[1]
        weight = item[2]
        r1 = node_ref[node1]
        r2 = node_ref[node2]
        mst_edges[(r1, r2)] = float(weight)
        print [r1, r2], "is an edge in mst"
    print "In total, the number of edges in mst is", len(mst_edges)
    print 

    mst_ugraph = {}
    for ref in alms:
        mst_ugraph[ref] = set([])
    for r1, r2 in mst_edges:
        mst_ugraph[r1].add(r2)
        mst_ugraph[r2].add(r1)

   
    #divide tree into forest

    ref_confs = {}
    for ref in alms:
        ref_confs[ref] = 0
        for qry in alms[ref]:
            conf = alms[ref][qry].confidence
            ref_confs[ref] += conf
    ref_confs_sort = sorted(ref_confs.items(), key=lambda x:x[1])
    visited = {}
    for ref, conf in ref_confs_sort:
        visited[ref] = False

    forest_ugraph = {}
    forest = {}
    while True:
        found = False
        for ref, conf in ref_confs_sort:
            if visited[ref] == False:
                found = True
                root = ref
                break
        if found == True: # there's still unvisited vertix in mst
            onetree = {}
            oneugraph = {}
            mst_traversal(mst_ugraph, visited, root, onetree, oneugraph) 
            forest_ugraph[root] = oneugraph
            forest[root] = onetree
        else:
            break

    num_vertices_forest = 0
    num_edges_forest = 0
    for root in forest:
        onetree = forest[root]
        num_vertices_forest += len(onetree)
        for v in onetree:
            num_edges_forest += len(onetree[v])
    print "In total, the number of vertices in forest is", num_vertices_forest
    print "In total, the number of edges in forest is", num_edges_forest
    print "In total, the number of trees in forest is", len(forest)

    num_vertices_forest_ugraph = 0
    num_edges_forest_ugraph = 0
    for root in forest_ugraph:
        oneugraph = forest_ugraph[root]
        num_vertices_forest_ugraph += len(oneugraph)
        for v in oneugraph:
            num_edges_forest_ugraph += len(oneugraph[v])
    print "In total, the number of vertices in forest_ugraph is", num_vertices_forest_ugraph
    print "In total, the number of edges in forest_ugraph is", num_edges_forest_ugraph

    # get orientation for each vertix in forest
    vertex_orientations = {}
    for root in forest:
        onetree = forest[root]
        onetree_v_oris = {}
        onetree_v_oris[root] = 1
        tree_traversal_orientation(edge_orientations, onetree_v_oris , root, onetree)
        vertex_orientations[root] = onetree_v_oris
       
    return forest, vertex_orientations


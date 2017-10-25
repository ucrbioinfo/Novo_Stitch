#!/usr/bin/python
import csv
import sys

from sets import Set

from collections import defaultdict

import collections
import math


def DFS_DAG(ug, g, visited, start, subgraph):
    visited[start] = True
    subgraph[start] = g[start]
    preds = ug[start]
    for v in preds:
        if visited[v] == True:
            continue
        DFS_DAG(ug, g, visited, v, subgraph)

def splite_graph(g):
    ug = {}
    for v in g:
        ug[v] = set([])
    for v in g:
        for u in g[v]:
            ug[v].add(u)
            ug[u].add(v)
    visited = {}
    for v in g:
        visited[v] = False
    new_g_set = []
    while True:
        found = False
        for v in g:
            if visited[v] == False:
                found = True
                start = v
                break
        if found == True:
            subgraph = {}
            DFS_DAG(ug, g, visited, start, subgraph)
            new_g_set.append(subgraph)
        else:
            break
    return new_g_set

def splite_graphs(graph_set):
    new_graph_list = []
    for root in graph_set:
        g = graph_set[root]
        new_g_set = splite_graph(g)
        for new_g in new_g_set:
            new_graph_list.append(new_g)
    return new_graph_list

def DFS(g, v, u, visited):#search for path from u to v
    if u == v:
        return True
    visited[u] = True
    for w in g[u]:
        if visited[w] == True:
            continue
        if DFS(g, v, w, visited) == True:
            return True
    return False       

def merge_DAG(g1, g2):
    for v in g2:
        if v not in g1:
            g1[v] = set([])
    for v in g2:
        next_nodes = g2[v]
        for u in next_nodes:
            if u in g1[v]:
                continue
            #check if cycle will be created by adding (v,u) to g1
            visited = {}
            for x in g1:
                visited[x] = False
            ifcycle = DFS(g1, v, u, visited)
            if ifcycle == False:
                g1[v].add(u)
            else:
                print "give up edge", v, u

def build_DAG(alms, ref, ori):
    DAG = {}
    for qry in alms[ref]:
        DAG[qry] = set([])
    for q1 in alms[ref]:
        x = alms[ref][q1]
        for q2 in alms[ref]:
            y = alms[ref][q2]
            if q1 == q2:
                continue
            if x.start < y.start and x.end < y.end and y.start < x.end:
                if ori == 1:
                    DAG[q1].add(q2)
                elif ori == -1:
                    DAG[q2].add(q1)
    return DAG

def tree_traversal_merge_DAGs(alms, onetree, onetree_v_oris, root, DAG):
    this_ori = onetree_v_oris[root]
    DAG_one = build_DAG(alms, root, this_ori)
    merge_DAG(DAG, DAG_one)
    #search children
    children = onetree[root] 
    for c in children:
        tree_traversal_merge_DAGs(alms, onetree, onetree_v_oris, c, DAG) 

def merge_DAGs(alms, forest, vertex_orientations):
    #merge DAGs

    DAG_set = {}
    for root in forest:
        DAG = {}
        onetree = forest[root]
        onetree_v_oris = vertex_orientations[root]
        tree_traversal_merge_DAGs(alms, onetree, onetree_v_oris, root, DAG)
        DAG_set[root] = DAG
    
    DAG_list = splite_graphs(DAG_set)
    
    print "In total, the number of connected DAG is", len(DAG_list)

    return DAG_list


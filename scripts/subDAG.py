#!/usr/bin/python
import csv
import sys

from sets import Set

from collections import defaultdict

import collections
import math

from BFS import dGraph
import BFS


def output_mtp_nodes(fo, mtp_sets):
    for v in mtp_sets:
        fo.write(str(v)+" ")
    fo.write("\n")    

def get_subDAGs(dgraphs_set, output_dir):
    #build dgraphs_pos_set and dgraphs_neg_set for calculating sinks and sources 
    dgraphs_pos_set = []
    dgraphs_neg_set = []
    G_set = []
    for g in dgraphs_set:
        G = dGraph()
        g_pos = {}
        g_neg = {}
        for v in g:
            g_pos[v] = {}
            g_neg[v] = {}
            G.add_vertex(v)
        for v in g:
            for u in g[v]:
                g_pos[v][u] = 1
                g_neg[u][v] = 1
                G.add_edge(v, u, 1)
        dgraphs_pos_set.append(g_pos)
        dgraphs_neg_set.append(g_neg)
        G_set.append(G)
    
    #get mtp contigs(greedy algorithm)
    fo = file(output_dir+"/mtp_nodes.log", 'w')
    mtp_node_set = set([])
    for i in range(0, len(dgraphs_pos_set)):
        source = set([])
        sink = set([])
        g_pos = dgraphs_pos_set[i]
        g_neg = dgraphs_neg_set[i]
        G = G_set[i]
        for v in g_pos:
            if g_pos[v] == {}:
                sink.add(v)
        for v in g_neg:
            if g_neg[v] == {}:
                source.add(v)
        max_cost = -1
        max_path = []
        for v in source:
            for u in sink:
                path = BFS.shortest_path(G, v, u)
                cost = len(path)
                if max_cost < cost:
                    max_path = path
                    max_cost = cost
                    #[max_v, max_u, max_cost, max_path] = [v, u, cost, path]

        min_set = set(max_path)
        while True:
            max_cost = -1
            max_path = []
            for v in source:
                if v in min_set:
                    continue
                min_cost = float('inf')
                min_path = []
                for x in min_set:
                    path = BFS.shortest_path(G, v, u)
                    if path == []:
                        cost = float('inf')
                    else:
                        cost = len(path)
                    cost = len(path)
                    if path == []:
                        continue
                    if min_cost > cost:
                        min_path = path
                        min_cost = cost
                if min_path == []:
                    continue
                if max_cost < min_cost:
                    max_cost = min_cost
                    max_path = min_path 
            for u in sink:
                if u in min_set:
                    continue
                min_cost = float('inf')
                min_path = []
                for x in min_set:
                    path = BFS.shortest_path(G, v, u)
                    if path == []:
                        cost = float('inf')
                    else:
                        cost = len(path)
                    if path == []:
                        continue
                    if min_cost > cost:
                        min_path = path
                        min_cost = cost
                if min_path == []:
                    continue
                if max_cost < min_cost:
                    max_cost = min_cost
                    max_path = min_path
            for x in max_path:
                if x not in min_set:
                    min_set.add(x)
            ifallin = True
            for v in source:
                if v not in min_set:
                    ifallin = False
                    break
            for u in source:
                if u not in min_set:
                    ifallin = False
                    break
            if ifallin == True:
                break
        output_mtp_nodes(fo, min_set) 
        for v in min_set:
            mtp_node_set.add(v)
    fo.close()
    print "In total, the number of contigs in mtp is", len(mtp_node_set)
    return mtp_node_set       






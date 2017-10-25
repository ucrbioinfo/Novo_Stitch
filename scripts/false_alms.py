import os
import vertex_cover_hyper_exhaust 
import vertex_cover_hyper_LP

def get_common(alms, r1, r2):
    common = set([])
    for qry in alms[r1]:
        if qry in alms[r2]:
            common.add(qry)
    return common

def edge_ornot(alms, r1, r2, q1, q2, dist_threshold):
    a11 = alms[r1][q1]
    a12 = alms[r1][q2]
    a21 = alms[r2][q1]
    a22 = alms[r2][q2]
    #1.check orientations
    ori_same = False
    ori_oppo = False
    if a11.orientation == a21.orientation and a12.orientation == a22.orientation:
        ori_same = True
    elif a11.orientation != a21.orientation and a12.orientation != a22.orientation:
        ori_oppo = True
    if ori_same == False and ori_oppo == False:
        return True
    #2.check distance
    mid11 = (a11.start+a11.end)/2
    mid12 = (a12.start+a12.end)/2
    mid21 = (a21.start+a21.end)/2
    mid22 = (a22.start+a22.end)/2
    dist1 = mid12 - mid11
    dist2 = mid22 - mid21

    if ori_same == True:
        if abs(dist2-dist1) < 1000:
            return False
        if abs(dist2-dist1) / max(abs(dist1), abs(dist2)) > dist_threshold: #0.5
            return True
    if ori_oppo == True:
        if abs(dist1+dist2) < 1000:
            return False
        if abs(dist1+dist2) / max(abs(dist1), abs(dist2)) > dist_threshold: #0.5
            return True
    return False

def DFS_ugraph(ug, visited, start, subgraph):
    visited[start] = True
    subgraph[start] = ug[start]
    preds = ug[start]
    for pred in preds:
        (v1, v2, v3) = pred
        if visited[v1] == False:
            DFS_ugraph(ug, visited, v1, subgraph)
        if visited[v2] == False:
            DFS_ugraph(ug, visited, v2, subgraph)
        if visited[v3] == False:
            DFS_ugraph(ug, visited, v3, subgraph)


def false_alms(GLPSOL, dist_threshold, alms, removed, output_dir):

    # build the graph
    print '---------------build alignment graph-------------------'
    # data structure
    vid_vertex = {}
    vertex_vid = {}
    edges = set([])
    vid_weight = {}

    # vertices
    num = 0
    for ref in alms:
        for qry in alms[ref]:
            vid_vertex[num] = (ref, qry)
            vertex_vid[ref, qry] = num
            x = alms[ref][qry]
            weight = x.confidence
            vid_weight[num] = weight
            print (ref, qry, weight), "is a vertex."
            num += 1
    print "The number of vertices for alignment graph is", len(vid_vertex)

    # edges        
#    dist_threshold = 0.2
    num_edges = 0
    num_total = 0
    for r1 in alms:
        for r2 in alms:
            if r1 >= r2:
                continue
            common = get_common(alms, r1, r2)
            if len(common) < 2:
                continue
            for q1 in common:
                for q2 in common:
                    if q1 >= q2:
                        continue
                    num_total += 1
                    if edge_ornot(alms, r1, r2, q1, q2, dist_threshold) == False: #check if the 4 alignments conflict
                        continue            
                    print (r1, q1), (r1, q2), (r2, q1), (r2, q2), "is a hyper edge."
                    edge = (vertex_vid[r1,q1], vertex_vid[r1,q2], vertex_vid[r2,q1], vertex_vid[r2,q2])
                    edges.add(edge)
                    num_edges += 1
    print "The number of hyper edges for alignment graph is", num_edges, "out of", num_total, "potential hyper edges."

    # store the undirected graph in ugraph
    ugraph = {}
    for vid in vid_vertex:
        ugraph[vid] = set([])
    for e in edges:   
        (v1, v2, v3, v4) = e
        ugraph[v1].add((v2, v3, v4))
        ugraph[v2].add((v1, v3, v4))
        ugraph[v3].add((v1, v2, v4))
        ugraph[v4].add((v1, v2, v3))

    # output
    fv = file(output_dir+"/alm_graph_vertices.txt", 'w')
    fv.write("vid\t(ref, qry)\tweights\n")
    for vid in vid_vertex:
        fv.write(str(vid)+"\t")
        fv.write("( "+str(vid_vertex[vid][0])+", "+str(vid_vertex[vid][1])+" )\t") 
        fv.write(str(vid_weight[vid]))
        fv.write("\n")
    fv.close()

    fe = file(output_dir+"/alm_graph_edges.txt", 'w')
    fe.write("vid_11\tvid_12\vid_21\t_vid_22\n")
    for v1, v2, v3, v4 in edges:
        fe.write(str(v1)+"\t"+str(v2)+"\t"+str(v3)+"\t"+str(v4)+"\n")
    fe.close()
    print "In total, the number of vertices in graph is", len(vid_vertex)
    print "In total, the number of hyper edges in graph is", len(edges)

    print '---------------END-------------------'

    print '---------------vertex cover-------------------'
    # divide ugraph into connected component
    visited = {}
    vid_list = []
    for vid in vid_vertex:
        visited[vid] = False
        vid_list.append(vid)
    vid_list.sort()

    ugraphs = []
    while True:
        found = False
        for vid in vid_list:
            if visited[vid] == False:
                found = True
                start = vid
                break
        if found == True: # there's still unvisited vertix in ugraph
            oneugraph = {}
            DFS_ugraph(ugraph, visited, start, oneugraph)
            ugraphs.append(oneugraph)
        else:
            break
    # call vertex cover algorithm
    LP_dir = output_dir + "/LP"
    if not os.path.isdir(LP_dir):
        os.makedirs(LP_dir)

    whole_cover = set([])
    for i in range(len(ugraphs)):
        ug = ugraphs[i]
        sub_vid_weight = {}
        for vid in ug:
            sub_vid_weight[vid] = vid_weight[vid]
        if len(ug) > 20:
            cover = vertex_cover_hyper_LP.VC(GLPSOL, LP_dir, ug, i, sub_vid_weight)   
        else:
            cover = vertex_cover_hyper_exhaust.VC(ug, sub_vid_weight)
        for v in cover:
            whole_cover.add(v)

    fc = file(output_dir+"/cover.log", 'w')
    fc.write("vid\t(ref, qry)\tweights\n")
    for vid in whole_cover:
        fc.write(str(vid)+"\t")
        fc.write("( "+str(vid_vertex[vid][0])+", "+str(vid_vertex[vid][1])+" )\t") 
        fc.write(str(vid_weight[vid]))
        fc.write("\n")
    fc.close()
    print '---------------END-------------------'
    
    # removing false positive alignments
    for vid in whole_cover:
        (ref, qry) = vid_vertex[vid]
        removed[ref, qry] = True





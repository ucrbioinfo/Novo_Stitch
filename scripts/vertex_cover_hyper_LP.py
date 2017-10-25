import os
import copy

def write_lp_file(filename, Vertices, Edges, w):
    fo = file(filename, 'w')

    # write objective function
    fo.write("Minimize\n")
    fo.write("value: ")
    for v in Vertices:
        fo.write("+ " + str(w[v]) + " x" + str(v) + "\n")   
    fo.write("\n")

    # write restrictions
    fo.write("Subject To\n")
    for e in Edges:
        (v1, v2, v3, v4) = e
        fo.write("x" + str(v1) + " + x" + str(v2) + " + x" + str(v3) + " + x" + str(v4) + " >= 1\n") 
    fo.write("\n")

    # write bounds
    fo.write("Bounds\n")
    V_list = []
    for v in Vertices:
        fo.write("x" + str(v)+" >= 0\n")
        V_list.append(v)
    fo.write("\n")

    fo.write("End\n")
    fo.close()

    return V_list
    
def read_results(filename, V_list, num_constrains):
    LP_values = {}
    with open(filename) as fi:
        fi.readline()
        line = fi.readline()
        line = line.strip()
        cols = line.split()
        opt_value = float(cols[2])
        for i in range(num_constrains): 
            fi.readline()
        for i in range(len(V_list)):
            v = V_list[i] 
            line = fi.readline()
            line = line.strip()
            cols = line.split()
            value = float(cols[1])
            LP_values[v] = value

    return LP_values

def greedy(LP_values, Edges):
    currentEdges = copy.deepcopy(Edges) 
    LP_values_sorted = sorted(LP_values.items(), key=lambda x: x[1], reverse=True)
    C = set([])
    for i in range(len(LP_values_sorted)):
        (v, value) = LP_values_sorted[i]
        if value < 0.25:
            continue
        ifremovenew = False 
        tempEdges = copy.deepcopy(currentEdges)
        for e in tempEdges:
            (v1, v2, v3, v4) = e
            if v1==v or v2==v or v3==v or v4==v: 
                ifremovenew = True
                currentEdges.remove(e)
        if ifremovenew == True:
            C.add(v)
    print currentEdges
    return C

def LP_based_alg(GLPSOL, output_dir, g, i, Vertices, Edges, w):#LP based Agortihm
    lp_file = output_dir + "/lp_"+str(i)+".lp"
    lp_result = output_dir + "/lp_"+str(i)+".sol"
#    GLPSOL = "glpsol"
    V_list = write_lp_file(lp_file, Vertices, Edges, w)    
    command = GLPSOL + " --lp " + lp_file + " -w " + lp_result
    os.system(command)
    LP_values = read_results(lp_result, V_list, len(Edges))
    C = greedy(LP_values, Edges)

    return C

def VC(GLPSOL, output_dir, g, i, w):
    Vertices = set([])
    Edges = set([])
    d = {}
    for v in g:
        Vertices.add(v)
        u_set = g[v]
        for u in u_set:
            e = [v, u[0], u[1], u[2]]
            e.sort()
            Edges.add((e[0], e[1], e[2], e[3]))
    if len(Edges) == 0:
        return set([])

    C = LP_based_alg(GLPSOL, output_dir, g, i, Vertices, Edges, w)
    return C


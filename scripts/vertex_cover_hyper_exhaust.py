def search(V_list, Edges, w, vertex_pos, opt_solution, opt_value, solution, m, n):
    if m == n:
        #check if solution is feasible
        feasible = True
        for (v1, v2, v3, v4) in Edges:
            if solution[vertex_pos[v1]] == False and solution[vertex_pos[v2]] == False and solution[vertex_pos[v3]] == False and solution[vertex_pos[v4]] == False:
                feasible = False
                break

        if feasible == False:
            return
        #compare with current best solution
        value = 0 #number of picked vertices
        for i in range(0, len(solution)):
            if solution[i] == True:
                value += w[V_list[i]]
        if value < opt_value[0]:
            for i in range(0, len(solution)):
                opt_solution[i] = solution[i]
            opt_value[0] = value
            
    else:
        solution[m] = False
        search(V_list, Edges, w, vertex_pos, opt_solution, opt_value, solution, m+1, n)
        solution[m] = True
        search(V_list, Edges, w, vertex_pos, opt_solution, opt_value, solution, m+1, n)

def exhaust(g, Vertices, Edges, w):#exhaust Agortihm
    V_list = []
    for v in Vertices:
        V_list.append(v)
    vertex_pos = {}
    for i in range(0, len(V_list)):
        vertex_pos[V_list[i]] = i    

    solution = []
    opt_solution = []
    for i in range(0, len(V_list)):
        solution.append(False)
        opt_solution.append(False)
    opt_value = [0]
    for v in w:
        opt_value[0] += w[v]
    opt_value[0] += 1

    search(V_list, Edges, w, vertex_pos, opt_solution, opt_value, solution, 0, len(V_list))

    C = set([])
    for i in range(0, len(V_list)):
        if opt_solution[i] == True:
            C.add(V_list[i])
    return C

def VC(g, w):
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
    C = exhaust(g, Vertices, Edges, w)
    return C



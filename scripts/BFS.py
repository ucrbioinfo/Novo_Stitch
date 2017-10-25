#!/usr/bin/python
from sets import Set
from collections import defaultdict
import collections
import math

from Queue import Queue

class dGraph:
	
    def __init__(self):
        self.vertices = set()
 
      # makes the default value for all vertices an empty list
        self.edges = collections.defaultdict(list)
        self.weights = {}
 
    def add_vertex(self, value):
        self.vertices.add(value)
 
    def add_edge(self, from_vertex, to_vertex, distance):
        if from_vertex == to_vertex: pass  # no cycles allowed
        self.edges[from_vertex].append(to_vertex)
        self.weights[(from_vertex, to_vertex)] = distance
 
    def __str__(self):
        string = "Vertices: " + str(self.vertices) + "\n"
        string += "Edges: " + str(self.edges) + "\n"
        string += "Weights: " + str(self.weights)
        return string
 
def dijkstra(graph, start):
  # initializations
    S = set()
 
  # delta represents the length shortest distance paths from start -> v, for v in delta. 
  # We initialize it so that every vertex has a path of infinity (this line will break if you run python 2)
    delta = dict.fromkeys(list(graph.vertices), float('inf'))
    previous = dict.fromkeys(list(graph.vertices), None)
 
  # then we set the path length of the start vertex to 0
    delta[start] = 0
 
  # while there exists a vertex v not in S
    while S != graph.vertices:
    # let v be the closest vertex that has not been visited...it will begin at 'start'
        v = min((set(delta.keys()) - S), key=delta.get)
 
    # for each neighbor of v not in S
        for neighbor in set(graph.edges[v]) - S:
            new_path = delta[v] + graph.weights[v,neighbor]
      # is the new path from neighbor through 
            if new_path < delta[neighbor]:
        # since it's optimal, update the shortest path for neighbor
                delta[neighbor] = new_path
 
        # set the previous vertex of neighbor to v
                previous[neighbor] = v
        S.add(v)
		
    return (delta, previous)

def BFS(graph, start):
    previous = dict.fromkeys(list(graph.vertices), None)
    visited = set()
    q = Queue()
    q.put(start)
    visited.add(start)
    while not q.empty():
        v = q.get()
        for u in set(graph.edges[v]):
            if u not in visited:
                q.put(u)
                visited.add(u)
                previous[u] = v

    return previous
    

def shortest_path(graph, start, end):
	
    previous = BFS(graph, start)
  
    path = []
    vertex = end
 
    while vertex is not None:
        path.append(vertex)
        vertex = previous[vertex]
 
    path.reverse()
    if path[0] == start:
        return path
    else:
        return []



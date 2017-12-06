from collections import defaultdict
from queue import PriorityQueue
from heapq import *

class Graph:
    def __init__(self):
        self.nodes = set()
        self.edges = defaultdict(list)
        self.distances = {}
    
    def add_node(self, value):
        self.nodes.add(value)
        
    def add_edge(self, from_node, to_node, dist):
        self.edges[from_node].append(to_node)
        self.edges[to_node].append(from_node)
        self.distances[(from_node, to_node)] = dist
        self.distances[(to_node, from_node)] = dist

def dijkstra(g, source, dest):
    path = { source:0 }
    q = PriorityQueue()
    parents = { source:None}
    for vertex in g.nodes:
        weight = float("inf")
        if vertex == source:
            weight = 0
        path[vertex] = weight
        parents[vertex] = None
        
    q.put(([0, source]))
    
    while not q.empty():
        v = q.get()[1]
        
        for u in g.edges[v]:
            cand = path[v] + g.distances[(u,v)]
            if path[u] > cand:
                path[u] = cand
                parents[u] = v
                if cand < -1000:
                    raise Exception("Negative cycle")
                q.put([path[u],u])
        
    shortest_path = []
    end = dest
    while end is not None:
        shortest_path.append(end)
        end = parents[end]
        
    shortest_path.reverse()
    
    return shortest_path, path[dest] 

def get_dijkstra_path(graph, from_node, to_node):
    result = dijkstra(graph, from_node, to_node)
    print("Shortest path from", from_node, "to", to_node,"=", result[0])
    print("length = ", result[1])


if __name__ == "__main__":
    nodes = ["A","B","C","D","E","F","G"]
    edges = [
        ("A", "B", 7),
        ("A", "D", 5),
        ("B", "C", 8),
        ("B", "D", 9),
        ("B", "E", 7),
        ("C", "E", 5),
        ("D", "E", 15),
        ("D", "F", 6),
        ("E", "F", 8),
        ("E", "G", 9),
        ("F", "G", 11)
    ]

    g = Graph()
    for i in nodes:
        g.add_node(i)
        
    for i in edges:
        g.add_edge(i[0],i[1],i[2])
        
    get_dijkstra_path(g, "A", "G")

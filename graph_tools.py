# -*- coding: utf-8 -*-
"""
Assortment of graph-theoretic tools to facilitate studying the MRCN problem

This file contains:

1) graph construction functions including:
      i) erdos-renyi random graphs
     ii) random connected graphs
    iii) cycle graphs
     iv) complete graphs
      v) uniform generalized star graphs
   and functions to read adjacency-lists from a text files create lists of
   different-sized graphs from the same class

2) functions to compute basic graph properties including:
     i) adjacency
    ii) degree
   iii) neighborhood
    iv) vertex orderings (e.g., by degree)
     v) induced subgraphs
   and a function to perform DFS and return the resultant ordering

3) a plotting function that draws a convex rectilinear plane embedding of a
   graph about the unit circle, showing edge crossings; features include even-
   spacing and vertex numbering. the resultant drawing is either displayed or
   saved as a picture file (.png)
"""
from random import random, choice
from graph import Graph


######################
# graph construction #
######################


def erdos_renyi(n, p=1/2):
    """
    constructs erdos-renyi random graph on n vertices with probability p
    (default p=1/2)
    """
    G = Graph(V=range(n), E=[], normalize=False)
    for i in G.V:
        for j in G.V:
            if i < j and random() < p:
                G.add_edge((i, j))
    return G


def random_tree(n):
    """
    constructs a random tree by sequentially adding new edges to a random
    vertex already part of the tree
    """
    G = Graph(V=[0], E=[])
    for i in range(1, n):
        G.add_edge((choice(G.V), i))
        G.add_vertex(i)
    return G


def random_connected(n, p=1 / 2):
    """
    constructs a random connected graph by first forming a tree on n vertices
    and then adding subsequent edges with probability p (default p=1/2)
    """
    G = random_tree(n)
    for i in G.V:
        for j in G.V:
            if i < j and (i, j) not in G and random() < p:
                G.add_edge((i, j))
    return G


def C(n):
    """
    cycle graph on n vertices

    note: M(C(n)) = (n(n - 3))/2 if n odd else n(n - 4)/2 + 1
    """
    E = [(k, (k + 1) % n) for k in range(n)]
    return Graph(E=E, V=range(n))


def K(n):
    """
    complete graph on n vertices

    note: M(K(n)) = nC4
    """
    E, vs = [], range(1, n)
    for n in range(n - 1):
        for k in vs:
            E.append((n, k))
        vs = vs[1:]
    return Graph(E=E, V=range(n))


def uniform_generalized_star(b, e):
    """
    uniform generalized star graph on b branches of length e
    (in number of edges)

    the uniform generalized star graph is constructed by adding paths of e - 1
    edges to the b branches of a star graph
    """
    E, V, s = [], [0], 1
    for b in range(b):
        E.extend([(0, s)] + [(s + k, s + k + 1) for k in range(e - 1)])
        V.extend(range(s, s + e))
        s += 10    # or some other sufficiently large offset
    return Graph(E=E, V=V)


def class_range(T, a, b):
    """
    creates a list of class T graphs (e.g., C or K) in range of sizes a,..,b
    including the endpoints
    """
    return [T(i) for i in range(a, b + 1)]


def from_file(file):
    """
    constructs graph from adjacency-list in text file; the accepted format
    is very simple: line i of the file contains the neighbors of node i,
    separated by spaces. e.g.,

    1 2 3
    0 2
    0 1
    0

    is the graph specified by V = {0,1,2,3}, E = {{0,1},{0,2},{0,3},{1,2}}
    """
    adjacency_list = open(file).read().splitlines()
    E, size = [], 0
    for v, neighborhood in enumerate(adjacency_list):
        neighbors = map(int, neighborhood.split(" "))
        for n in neighbors:
            E.append((v, n))
        size = max(size, max(neighbors))

    return Graph(V=range(size + 1), E=E)


def from_file_reg_format(file):
    """
    constructs graph from adjacency-list in text file; the accepted format
    is the one used at the following website:

    http://www.mathe2.uni-bayreuth.de/markus/reggraphs.html#CRG

    the above website has many text file adjacency-lists for regular graphs;
    such graphs are condusive to study in the MRCN problem
    """
    import re

    lines = open(file).read().splitlines()
    graphs, E = [], []
    i = j = M = 0
    b = False

    for line in lines:
        if re.match(r'Gr', line):
            i += 1
            b = True
            if i > 1:
                graphs.append(Graph(E=E, V=range(M)))
                j = M = 0
                E[:] = []
        elif line == "":
            pass
        elif re.match(r'[GOT]', line):
            b = False
        elif b:
            for c in line[4:]:
                if c == ' ':
                    continue
                elif int(c) - 1 > j:
                    E.append((j, int(c) - 1))
                M = max(M, int(c))
            j += 1
    graphs.append(Graph(E=E, V=range(M)))  # get last graph
    return graphs


#############################################
# graph-theoretic properties and algorithms #
#############################################


def adj(e, f):
    """
    checks if edges e, f are adjacent
    """
    return e[0] == f[0] or e[0] == f[1] or e[1] == f[0] or e[1] == f[1]


def deg(v, G):
    """
    deg(v) for v ∈ V(G)
    """
    return len([e for e in G.E if v in e])


def N(v, G):
    """
    neighborhood of v, i.e., {u | (u,v) ∈ E(G)}
    """
    return [u for u in G.V if (u, v) in G.E or (v, u) in G.E]


def ordered(G, b=False):
    """
    returns list of vertices ordered by degree

    b: specifies ascending or descending, where default is ascending (False)
    """
    return [s[0] for s in sorted(map(lambda v: (v, deg(v, G)), G.V),
                                 key=lambda p: p[1], reverse=b)]


def induced(G, U):
    """
    returns the subgraph induced by U ⊆ V(G), i.e., {(u,v) | u,v ∈ U}
    """
    return [e for e in G.E if set(e) <= set(U)]


def N_induced(G, U, v):
    """
    returns the edges induced by the neighborhood of v in U ⊆ V(G), i.e.,
    {(u,v) | u ∈ U}
    """
    return [e for e in G.E if (e[0] in U and e[1] == v) or
                              (e[0] == v and e[1] in U)]


def dfs(G):
    """
    returns a DFS-induced ordering of V as well as the induced predecessor
    list as a dictionary

    note:
    the ordering of the gray list corresponds to discovery time, and that of
    the black list to finishing time (see CLRS 3rd Ed. pgs. 603-612)

    TODO dfs without recursion and bfs
    """
    black, gray, prev = [], [], {}

    def _dfs(v):
        gray.append(v)    # discovery time
        for w in N(v, G):
            if w not in gray and w not in black:
                prev[w] = v
                _dfs(w)
        black.append(v)    # finishing time

    for u in G.V:    # necessary only for disconnected graphs
        if u not in gray and u not in black:
            _dfs(u)

    return black, prev


############
# plotting #
############


def circle_plot(G, d=None, save_path='', circle=None):
    """
    plots `G` about the unit circle, with the vertices being placed according
    to the order of permutation P; fill in edges

    name: MRCN algorithm used
    path: desired path to save file
    """
    # silence annoying matplotlib warnings
    import warnings
    warnings.filterwarnings("ignore")

    from math import cos, sin, pi
    import matplotlib.pyplot as plt
    import networkx as nx

    G_nx = nx.Graph()
    G_nx.add_nodes_from(G.V)
    G_nx.add_edges_from(G.E)

    radius = 1
    pos = {}
    alpha = (2 * pi) / len(G.V)
    for j, v in enumerate(reversed(d or G.V)):
        # plot vertices clockwise, as equidistant points on circle
        pos[v] = -radius * sin((j + 1) * alpha), radius * cos((j + 1) * alpha)

    nx.draw(G_nx, pos, node_size=800)
    nx.draw_networkx_labels(G_nx, pos, font_size=16)

    if circle:
        circle = plt.Circle((0, 0), radius, color=circle, fill=False)
        plt.gcf().gca().add_artist(circle)

    plt.axis('off')
    plt.axis('equal')
    if save_path:
        plt.savefig(save_path)
    else:
        plt.show()
    plt.close('all')

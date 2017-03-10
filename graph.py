"""
This file contains the Graph class, which represents undirected graphs. in
addition to list structures to store vertices and edges, we augment the Graph
class with functionality to facilitate usage with MRCN algorithms
"""
from random import shuffle, sample


class Graph:
    """
    Representation of an undirected graph

    For a graph G = (V, E), letting n=|V|, V is represented as the list [n-1]
    and E as a list of pairs corresponding to the edges of E.

    Graphs are stored in a "normalized" format, meaning V is kept sorted, and
    E is sorted on both the first and second element of the edges (note that
    this does NOT imply directed edges, it's a convenience for processing)

    We also keep a list D which represents a convex drawing of G in the plane;
    such a drawing is completely determined by the relative order of the
    vertices, and hence, it suffices to keep an ordered list D of the vertices.
    D is mainly used in MRCN applications. Also note that D may change often
    since we process different convex drawings as a core part of our MRCN
    algorithms, as opposed to the "more static" V and E
    """

    def __init__(self, V=None, E=None, normalize=True):
        """
        accepts pair of lists (V,E) to construct a graph
        """
        self.V, self.E = sorted(V), (E if E else [])
        self._normalize(self.E) if normalize else E
        self.D = self.V[:]

    def _normalize(self, E=None):
        """
        maintains the current or argument edge-set E in sorted order by
        first and second element of an edge
        """
        if not E:  # normalize current edge-set
            for i, e in enumerate(self.E):
                self.E[i] = tuple(sorted(e))
            self.E = sorted(self.E, key=lambda x: (x[0], x[1]))
        else:  # normalize edge-set argument E
            for i, e in enumerate(E):
                E[i] = tuple(sorted(e))
            return sorted(E, key=lambda x: (x[0], x[1]))

    def __len__(self):
        """
        computes size of G, i.e., number of vertices
        """
        return len(self.V)

    # TODO
    # when adding a new edge or vertex, perform binary search to place it
    # correctly in proper order rather than sorting V or E each time

    def add_edge(self, e, normalize=True):
        """
        adds edge e to E
        """
        norm = tuple(sorted(e))
        if norm not in self:
            self.E.append(norm)
            self.add_vertices(norm)
        if normalize:
            self._normalize()

    def add_edges(self, edges):
        """
        adds list of edges to E
        """
        for e in edges:
            self.add_edge(e)

    def add_vertex(self, v):
        """
        adds vertex v to V
        """
        if v not in self:
            self.V = sorted(self.V + [v])
            self.D = self.V[:]
            # else warn

    def add_vertices(self, vertices):
        """
        adds list of vertices to V
        """
        for v in vertices:
            self.add_vertex(v)

    def permute(self, in_place=False):
        """
        returns random permutation of G.V, unless `in_place=True` in which
        case permute `self.D` in place
        """
        if in_place:
            shuffle(self.D)
        else:
            return sample(self.V, len(self))

    def reset(self):
        """
        resets permutation order on vertices to original ordering (1,...,n)
        """
        self.D = self.V[:]

    def next(self):
        """
        generates next lexicographic permutation of self.D

        pseudo-code:
        1. find largest k with self.D[k] < self.D[k + 1]
           if no such k exists, this permutation is the last
        2. find largest l with self.D[k] < self.D[l]
           since this is true for k+1, l is well-defined and k < l
        3. swap self.D[k] and self.D[l]
        4. reverse self.D[k + 1],...,self.D[n]
        """
        i = len(self.D) - 2
        while not (i < 0 or self.D[i] < self.D[i + 1]):
            i -= 1
        if i < 0:
            return False
        j = len(self.D) - 1
        while not (self.D[j] > self.D[i]):
            j -= 1
        self.D[i], self.D[j] = self.D[j], self.D[i]
        self.D[i + 1:] = reversed(self.D[i + 1:])
        return True

    def __contains__(self, obj):
        """
        membership for vertices and edges
        """
        if type(obj) == int:  # vertex
            return True if obj in self.V else False
        else:    # edge
            return True if obj in self.E or obj[::-1] in self.E else False

    def __str__(self):
        """
        displays V, E
        """
        return 'V - ' + str(self.V) + '\nE - ' + str(self.E)

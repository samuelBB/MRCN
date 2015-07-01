"""
Algorithms for the Maximum Rectilinear Crossing Number (MRCN) problem

=====

MRCN is the optimization problem that seeks to maximize the number of edge
crossings in rectilinear plane embeddings of graphs, which we simply call
drawings. In symbols, let G = (V,E) and let D(G) be the set of all possible
drawings for G. Let CR(D) for D ∈ D(G) be the maximum number of edge crossings
for D. We wish to compute the MRCN M_G, defined as

                         M_G = max_{D ∈ D(G)} CR(D)

Another important quantity is defined as follows: let D˚(G) be the set of all
possible convex drawings for G. Let CR˚(D) for D ∈ D˚(G) be the maximum number
of edge crossings for D. We are also interested in computing the convex-MRCN,
defined as

                        M˚_G = max_{D ∈ D˚(G)} CR˚(D)

This latter quantity is more readily computed and much easier to represent in
code.

This file contains:

1) a naive algorithm to compute the convex-MRCN

2) approximation algorithms for MRCN, including:
      i) naive randomized (random permutation and sequential)
     ii) blind greedy (with multiple vertex-ordering heuristics)
    iii) local search (moving vertices only)

3) functions to facilitate testing the above algorithms

Note that all of the algorithms in this file work in the setting of convex
drawings. This is a very convenient simplification since we can easily
represent such drawings as permutations on the vertex-sets, rather than as
sets of xy-coordinates as demanded by the general setting.
"""


from random import shuffle, randint
from GraphTheoreticTools import adjacent, ordered, plot, \
    neighborhoodInducedSubgraph


###################################
# naive algorithm for convex-MRCN #
###################################


def NAEPs(G=None, E=None):
    """
    constructs list of nonadjacent edge pairs (NAEPs) for graph G or edge-set E
    """
    if G:
        E = G.E
    pairs = []
    for i in range(len(E)):
        for j in range(i + 1, len(E)):
            if not adjacent(E[i], E[j]):
                pairs += [(E[i], E[j])]
    return pairs


def crosses(D, ep):
    """
    checks if edge-pair ep crosses in drawing D
    """

    # map the vertices of the edge-pair to their location in the drawing D
    a, b, c, d = list(map(D.index, [ep[0][0], ep[0][1], ep[1][0], ep[1][1]]))

    # check vertex locations for crossing
    return (a < c < b < d) or (b < c < a < d) or (a < d < b < c) or \
           (b < d < a < c) or (c < a < d < b) or (d < a < c < b) or \
           (c < b < d < a) or (d < b < c < a)


def CR(G=None, D=None, EP=None):
    """
    for graph G, computes CR˚(D) for drawing D ∈ D˚(G)
    """
    if G:
        D = G.D
    cr = 0
    for ep in EP:
        if crosses(D, ep):
            cr += 1
    return cr


def maxCR(G):
    """
    computes M˚_G for graph G

    let n = |V|. non-asymptotic speedups used (both follow by symmetry):
        1. keep the first vertex fixed, and stop permuting when the second
           vertex in the permutation is v_n, the highest-indexed vertex
        2. only consider permutations p where p(2) < p(n)

    The trivial algorithm checks every permutation in n! time. The above
    optimizations give:
        1. the while-loop executes (n-2)*(n-2)! times
        2. the if-statement executes (n-1)!/2 times
    """
    max_cr, maxDrawing, moreDrawingsToCheck, EP = 0, list(), True, NAEPs(G)
    while moreDrawingsToCheck and not G.D[1] == G.V[G.size - 1]:
        if G.D[1] < G.D[G.size - 1]:
            currentMax = CR(G=G, EP=EP)
            if max_cr < currentMax:
                max_cr, maxDrawing = currentMax, G.D[:]
        moreDrawingsToCheck = G.next()
    return max_cr, maxDrawing


def MMCR(Gs):
    """finds the graph(s) with minimum and maximum convex-MRCN over a set of
    graphs Gs

    TODO:
    -support ties as a list of optimums
    -make neater, if possible
    """
    maxIndex = minIndex = 0
    max_cr, maxDrawing = maxCR(Gs[0])
    min_cr, minDrawing = max_cr, maxDrawing[:]
    for i, G in enumerate(Gs[1:]):
        currentCR, currentDrawing = maxCR(G)
        if max_cr < currentCR:
            maxIndex, max_cr, maxDrawing = i, currentCR, currentDrawing[:]
        if min_cr > currentCR:
            minIndex, min_cr, minDrawing = i, currentCR, currentDrawing[:]
    return maxIndex, max_cr, maxDrawing, minIndex, min_cr, minDrawing


############################
# approximation algorithms #
############################


def randomized(G):
    """
    returns CR˚(D) for a random permutation D ∈ D˚(G)

    note:
    this yields OPT/3 in expectation (and can be derandomized)
    """
    G.permute()
    return CR(G=G, EP=NAEPs(G)), G.D


def randomizedStepped(G):
    """
    returns CR˚(D) for a random permutation D ∈ D˚(G) constructed by randomly
    selecting vertices from V(G) one-by-one
    """
    V = G.V[:]
    D = list()
    while V:
        D += [V.pop(randint(0, len(V) - 1))]
    return CR(D=D, EP=NAEPs(G))


def greedy(G, custom=None, order=False, direction=False, mix=False, io=print):
    """
    sequentially places the vertices about the unit circle, maximizing the
    number of edge crossings for each placement

    custom:    specify a custom ordering
    order:     order vertices by degree if true
    direction: order by ascending degree if False (default) else descending
    mix:       specify if ordering should be shuffled
    io:        specify IO function
    """
    D, E = [], []

    # sequential ordering may already be random (e.g., for random graphs)
    if mix:
        shuffle(G.V)

    for v in ordered(G, direction) if order else (custom if custom else G.V):
        max_for_vertex = max_pos = 0
        E += neighborhoodInducedSubgraph(G, D, v)

        for tentative_pos in range(len(D)):
            tentative_D = D[:tentative_pos + 1] + [v] + D[tentative_pos + 1:]
            tentative_cr = CR(D=tentative_D, EP=NAEPs(E=E))
            # '<' returns first position giving max
            # '<=' returns last position giving max
            if max_for_vertex < tentative_cr:
                max_for_vertex, max_pos = tentative_cr, tentative_pos

        D.insert(max_pos + 1, v)
        # print result of each greedy round
        io('\n%d -- %s | cr: %d' % (v, D, max_for_vertex))

    return CR(D=D, EP=NAEPs(E=E)), D


def localSearch(G, mix=False):
    """
    starting with a random convex drawing, moves vertices while there are still
    moves that strictly increase the number of edge crossings

    TODO:
    -vertex-ordering heuristics (e.g., degree), if useful
    -printing partial results, if useful
    """

    # sequential ordering may already be random
    if mix:
        shuffle(G.V)

    EP = NAEPs(G)
    curr_cr, curr_d = CR(D=G.D, EP=EP), G.D[:]
    for v in G.V:
        for tentative_pos in range(len(curr_d)):
            tentative_D = \
                curr_d[:tentative_pos + 1] + [v] + curr_d[tentative_pos + 1:]
            tentative_cr = CR(D=tentative_D, EP=EP)
            if curr_cr < tentative_cr:
                curr_cr, curr_d = tentative_cr, tentative_D[:]

    return CR(D=curr_d, EP=EP), curr_d


###########
# testing #
###########


def test(G, algorithm, *args, io=print, file=None):
    """
    runs supplied MRCN algorithm, prints results, and saves resultant drawing
    to a .png picture file
    """
    io('\n\n%s:\nGraph:\n%s' % (algorithm.__name__, G))
    cr, d = algorithm(G, *args, io=io) if args else algorithm(G, io=io)
    plot(G, d, algorithm.__name__, path=file)
    io('\nResults:\n%s %s' % (cr, d))
    return cr, d


def randomizedTest(G, N, io=print, file=None):
    """
    runs randomized MRCN algorithm N times, returning the result of the
    maximal round; prints results, and saves resultant drawing to a .png
    picture file
    """
    io('randomized:\nGraph:\n%s' % G)
    max_cr, max_d = randomized(G)
    for i in range(N):
        # this shouldn't matter since each round is randomized
        # G.reset()
        curr, curr_d = randomized(G)
        if curr > max_cr:
            max_cr, max_d = curr, curr_d
    plot(G, max_d, randomized.__name__, path=file)
    io('\nResults:\n%s %s' % (max_cr, max_d))

"""
Algorithms for the Maximum Rectilinear Crossing Number (MRCN) problem

=====

MRCN is the optimization problem that seeks to maximize the number of edge
crossings in rectilinear plane embeddings of graphs, which we simply call
drawings. In symbols, let G = (V,E) and let D(G) be the set of all possible
drawings for G. Let CR(D) for D ∈ D(G) be the maximum number of edge crossings
for D. We wish to compute the MRCN M(G), defined as

                         M(G) = max_{D ∈ D(G)} CR(D)

Another important quantity is defined as follows: let D˚(G) be the set of all
possible convex drawings for G. Let CR˚(D) for D ∈ D˚(G) be the maximum number
of edge crossings for D. We are also interested in computing the convex-MRCN,
defined as

                        M˚(G) = max_{D ∈ D˚(G)} CR˚(D)

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
from GraphTheoreticTools import adj, ordered, circle_plot, \
    neighborhoodInduced


###################################
# naive algorithm for convex-MRCN #
###################################


def NAEPs(G=None, E=None):
    """
    constructs list of nonadjacent edge pairs (NAEPs) for graph G or edge-set E
    """
    if G:
        E = G.E
    naeps = []
    for i in range(len(E)):
        for j in range(i + 1, len(E)):
            if not adj(E[i], E[j]):
                naeps += [(E[i], E[j])]
    return naeps


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


def maxConvexCR(G):
    """
    computes M˚(G) for graph G

    let n = |V|. non-asymptotic speedups used (both follow by symmetry):
        1. keep the first vertex fixed, and stop permuting when the second
           vertex in the permutation is v_n, the highest-indexed vertex
        2. only consider permutations p where p(2) < p(n)

    The trivial algorithm checks every permutation in n! time. The above
    optimizations give:
        1. the while-loop executes (n-2)*(n-2)! times
        2. the if-statement executes (n-1)!/2 times
    """
    max_cr, maxDrawing = 0, []
    moreDrawingsToCheck = True
    EP = NAEPs(G)
    while moreDrawingsToCheck and not G.D[1] == G.V[G.size - 1]:
        if G.D[1] < G.D[G.size - 1]:
            currentMax = CR(G=G, EP=EP)
            if max_cr < currentMax:
                max_cr, maxDrawing = currentMax, G.D[:]
        moreDrawingsToCheck = G.next()
    return max_cr, maxDrawing


# TODO - support ties as a list of optimums; make cleaner, if possible
def MMCR(Gs):
    """finds the graph(s) with minimum and maximum convex-MRCN over a set of
    graphs Gs
    """
    maxIndex = minIndex = 0
    max_cr, maxDrawing = maxConvexCR(Gs[0])
    min_cr, minDrawing = max_cr, maxDrawing[:]
    for i, G in enumerate(Gs[1:]):
        currentCR, currentDrawing = maxConvexCR(G)
        if max_cr < currentCR:
            maxIndex, max_cr, maxDrawing = i, currentCR, currentDrawing[:]
        if min_cr > currentCR:
            minIndex, min_cr, minDrawing = i, currentCR, currentDrawing[:]
    return maxIndex, max_cr, maxDrawing, minIndex, min_cr, minDrawing


###########
# testing #
###########


def test(algorithm):
    """
    decorator which adds IO functionality to MRCN algorithms; in addition to
    runnning algorithm, results are printed with specified IO function, and the
    resultant drawing is plotted and saved to a .png picture file
    """
    def _test(*args, io=print, path=None, **kwargs):
        """
        the decorated function to be returned
        """
        io('\n~%s~\nGraph:\n%s' % (algorithm.__name__, args[0]))
        cr, d = algorithm(*args, **kwargs)
        circle_plot(args[0], d, algorithm.__name__, path=path)
        io('\nResults:\n%s %s\n' % (cr, d))
        return cr, d
    return _test


def randomizedTest(G, N, io=print, path=None):
    """
    runs randomized MRCN algorithm N times, returning the result of the
    maximal round; prints results, and saves resultant drawing to a .png
    picture file at desired path
    """
    io('randomized:\nGraph:\n%s' % G)
    max_cr, max_d = randomized(G)
    for i in range(N):
        curr, curr_d = randomized(G)
        if curr > max_cr:
            max_cr, max_d = curr, curr_d[:]
    circle_plot(G, max_d, randomized.__name__, path=path)
    io('\nResults:\n%s %s' % (max_cr, max_d))
    return max_cr, max_d


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
    D = G.D[:]
    G.reset()
    return CR(D=D, EP=NAEPs(G)), D


def randomizedStepped(G):
    """
    returns CR˚(D) for a random permutation D ∈ D˚(G) constructed by randomly
    selecting vertices from V(G) one-by-one
    """
    D = []
    while G.D:
        D += [G.D.pop(randint(0, len(G.D) - 1))]
    G.reset()
    return CR(D=D, EP=NAEPs(G)), D


@test
def greedy(G, custom=None, order=False, direction=False, mix=False, IO=None):
    """
    sequentially places the vertices about the unit circle, maximizing the
    number of edge crossings for each placement

    G:         the graph to run greedy on
    custom:    specify a custom ordering
    order:     order vertices by degree if true
    direction: order by ascending degree if False (default) else descending
    mix:       specify if vertex-ordering should be shuffled
    IO:        specify IO function for partial results, if desired
    """
    D, E = [], []

    # sequential ordering may already be random (e.g., for random graphs)
    if mix:
        shuffle(G.D)

    maxCR = 0
    for v in ordered(G, direction) if order else (custom if custom else G.D):
        maxSpot = 0    # affects placement, e.g., if no crossings found
        E += neighborhoodInduced(G, D, v)
        for currentSpot in range(len(D)):
            currentDrawing = D[:currentSpot + 1] + [v] + D[currentSpot + 1:]
            currentCR = CR(D=currentDrawing, EP=NAEPs(E=E))
            # '<' / '<=' returns first / last position giving max
            if maxCR < currentCR:
                maxCR, maxSpot = currentCR, currentSpot
        D.insert(maxSpot + 1, v)

        # print result of each greedy round
        if IO:
            IO('\n%d -- %s | cr: %d' % (v, D, maxCR))

    G.reset()
    return maxCR, D


# TODO (if useful)
# investigate the approximation gain and time increase with cap
# vertex-ordering heuristics (e.g., degree)
# printing partial results
@test
def localSearch(G, mix=False, cap=float('inf')):
    """
    starting with a random convex drawing, moves vertices while there are still
    moves that strictly increase the number of edge crossings

    G:   the graph to run local-search on
    mix: specify if vertex-ordering should be shuffled
    cap: number of times to loop through the vertices
    """

    # sequential ordering may already be random
    if mix:
        shuffle(G.D)

    EP = NAEPs(G)
    maxCR, D = CR(D=G.D, EP=EP), G.D[:]
    previousCR = -1
    while previousCR != maxCR and cap > 0:
        previousCR = maxCR
        cap -= 1
        for v in G.V:
            for currentSpot in range(len(D)):
                currentDrawing = D[:currentSpot + 1] + [v] + D[currentSpot + 1:]
                currentCR = maxCR(D=currentDrawing, EP=EP)
                if maxCR < currentCR:
                    maxCR, D = currentCR, currentDrawing[:]

    G.reset()
    return maxCR, D

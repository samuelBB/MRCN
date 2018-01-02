# -*- coding: utf-8 -*-
"""
Algorithms for the Maximum Rectilinear Crossing Number (MRCN) problem

MRCN is the optimization problem that seeks to maximize the number of edge
crossings in rectilinear plane embeddings of graphs, which we simply call
drawings. In symbols, let G = (V,E) and let D(G) be the set of all possible
drawings for G. Let CR(D) for D ∈ D(G) be the maximum number of edge crossings
for D. We wish to compute the MRCN M(G), defined as

                         M(G) = max_{D ∈ D(G)} CR(D)

Another important quantity is defined as follows: let D˚(G) be the set of all
possible convex drawings for G. Let CR˚(D) for D ∈ D˚(G) be the maximum number
of edge crossings for D. We are also interested in computing the convex-MRCN:

                        M˚(G) = max_{D ∈ D˚(G)} CR˚(D)

This latter quantity is more readily computed and represented simply in code.

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
from __future__ import print_function
from os.path import join
from collections import defaultdict
from random import random, shuffle, randint


from graph_tools import adj, ordered, circle_plot, N_induced


###################################
# naive algorithm for convex-MRCN #
###################################


def naeps(G=None, E=None):
    """
    constructs list of nonadjacent edge pairs (NAEPs) for graph G or edge-set E
    """
    if G:
        E = G.E
    naeps = []
    for i in range(len(E)):
        for j in range(i + 1, len(E)):
            if not adj(E[i], E[j]):
                naeps.append((E[i], E[j]))
    return naeps


def crosses(D, ep):
    """
    checks if edge-pair ep crosses in drawing D
    """

    # map the vertices of the edge-pair to their location in the drawing D
    a, b, c, d = [D.index(v) for v in (ep[0][0], ep[0][1], ep[1][0], ep[1][1])]

    # check vertex locations for crossing
    return (a < c < b < d) or (b < c < a < d) or (a < d < b < c) or \
           (b < d < a < c) or (c < a < d < b) or (d < a < c < b) or \
           (c < b < d < a) or (d < b < c < a)


def CR(G=None, D=None, EP=None):
    """
    for graph G, computes CR˚(D) for drawing D ∈ D˚(G)
    """
    if G and not D:
        D = G.D
    if G and not EP:
        EP = naeps(G)
    cr = 0
    for ep in EP:
        if crosses(D, ep):
            cr += 1
    return cr


def MRCN_convex(G):
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
    EP, max_cr, max_drawing = naeps(G), 0, []
    with G.clean:
        while not G.D[1] == G.V[len(G) - 1]:
            if G.D[1] < G.D[len(G) - 1]:
                current_max = CR(G=G, EP=EP)
                if max_cr < current_max:
                    max_cr, max_drawing = current_max, G.D[:]
            G.next()
    return max_cr, max_drawing


def crossing_spectrum(G):
    """
    return dict of `count -> drawing` for all counts
    """
    spectrum = defaultdict(list)
    EP = naeps(G)
    with G.clean:
        while not G.D[1] == G.V[len(G) - 1]:
            if G.D[1] < G.D[len(G) - 1]:
                spectrum[CR(G=G, EP=EP)].append(G.D[:])
            G.next()
    return spectrum


def MMCR(Gs):
    """
    finds the graph(s) with min and max convex-MRCN over a set of graphs
    TODO support ties as a list of optimums; make cleaner, if possible
    """
    max_index = min_index = 0
    max_cr, max_drawing = MRCN_convex(Gs[0])
    min_cr, min_drawing = max_cr, max_drawing[:]
    for i, G in enumerate(Gs[1:]):
        current_CR, current_drawing = MRCN_convex(G)
        if max_cr < current_CR:
            max_index, max_cr, max_drawing = i, current_CR, current_drawing[:]
        if min_cr > current_CR:
            min_index, min_cr, min_drawing = i, current_CR, current_drawing[:]
    return max_index, max_cr, max_drawing, min_index, min_cr, min_drawing


###########
# testing #
###########


def test(algorithm):
    """
    decorator which adds IO functionality to MRCN algorithms; in addition to
    runnning algorithm, results are printed with specified IO function, and the
    resultant drawing is plotted and saved to a .png picture file
    """
    def _test(*args, **kwargs):
        """
        the decorated function to be returned
        """
        io = kwargs.get('io', print)
        path = kwargs.get('path', '')
        io('\n~%s~\nGraph:\n%s' % (algorithm.__name__, args[0]))
        cr, d = algorithm(*args, **kwargs)
        circle_plot(args[0], d, join(path, algorithm.__name__))
        io('\nResults:\n%s %s\n' % (cr, d))
        return cr, d
    return _test


def randomized_test(G, N, io=print, path=''):
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
    circle_plot(G, max_d, join(path, randomized.__name__))
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
    D = G.permute()
    return CR(G, D=D), D


def randomized_stepped(G):
    """
    returns CR˚(D) for a random permutation D ∈ D˚(G) constructed by randomly
    selecting vertices from V(G) one-by-one
    """
    D = []
    while G.D:
        D.append(G.D.pop(randint(0, len(G.D) - 1)))
    G.reset()
    return CR(G=G, D=D), D


def randomized_continuous(G):
    """
    returns CR˚(D) for a random permutation D ∈ D˚(G) constructed by placing
    each vertex at a random angle about the unit circle

    note:
    this yields OPT/3 in expectation (and can be derandomized)
    """
    angles = []
    for _ in G.V:
        angle = random()
        while angle in angles:
            angle = random()
        angles.append(angle)
    D = [angles.index(a) for a in sorted(angles)]
    return CR(G=G, D=D), D


# TODO clean and add derandomized (anything else?)


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

    max_CR = 0
    for v in ordered(G, direction) if order else (custom if custom else G.D):
        max_spot = 0    # affects placement, e.g., if no crossings found
        E.extend(N_induced(G, D, v))
        for current_spot in range(len(D)):
            current_drawing = D[:current_spot + 1] + [v] + D[current_spot + 1:]
            current_CR = CR(D=current_drawing, EP=naeps(E=E))
            # '<' / '<=' returns first / last position giving max
            if max_CR < current_CR:
                max_CR, max_spot = current_CR, current_spot
        D.insert(max_spot + 1, v)

        # print result of each greedy round
        if IO:
            IO('\n%d -- %s | cr: %d' % (v, D, max_CR))

    G.reset()
    return max_CR, D


@test
def local_search(G, mix=False, cap=float('inf')):
    """
    starting with a random convex drawing, moves vertices while there are still
    moves that strictly increase the number of edge crossings

    G:   the graph to run local-search on
    mix: specify if vertex-ordering should be shuffled
    cap: number of times to loop through the vertices
    """

    # sequential ordering may already be random
    if mix:
        G.permute(in_place=True)

    EP = naeps(G)
    max_CR, D = CR(G, EP=EP), G.D[:]
    previous_CR = -1
    while previous_CR != max_CR and cap > 0:
        previous_CR = max_CR
        cap -= 1
        for v in G.V:
            for current_spot in range(len(D)):
                current_drawing = D[:current_spot+1] + [v] + D[current_spot+1:]
                current_CR = max_CR(D=current_drawing, EP=EP)
                if max_CR < current_CR:
                    max_CR, D = current_CR, current_drawing[:]

    G.reset()
    return max_CR, D

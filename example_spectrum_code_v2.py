from contextlib import closing
import matplotlib.backends.backend_pdf as mpdf

from MRCN_algorithms import crossing_spectrum
from graph_tools import erdos_renyi, C, K, circle_plot


# cycle graph on 8 vertices
G = C(9)

# generate spectrum of crossings
spectrum = crossing_spectrum(G)

# print basic info (min/max crossings, crossings missed)
keys = set(spectrum.keys())
print('Min crossings: %s, Max crossings: %s' % (min(keys), max(keys)))
print('Crossings missed: %s' % sorted(set(range(max(keys))) - keys))


with closing(mpdf.PdfPages('spectrum.pdf')) as pdf:
    # plot with title in sorted order
    for cr in sorted(spectrum):
        title = 'crossings=%s, #drawings=%s' % (cr, len(spectrum[cr]))
        circle_plot(G, spectrum[cr][0], title=title, show=False, pdf=pdf)
from MRCN_algorithms import crossing_spectrum
from graph_tools import erdos_renyi, C, K, circle_plot

# cycle graph on 8 vertices
# substitute your graph of choice
G = C(8)

# generate spectrum of crossings
spectrum = crossing_spectrum(G)

print('#Crossings: #Drawings')
for cr, drawings in spectrum.items():
    print('{cr}: {n_cr}'.format(cr=cr, n_cr=len(drawings)))

# example: draw first drawing found with 15 crossings
idx, n_crossings = 0, 15
circle_plot(G, spectrum[n_crossings][idx])


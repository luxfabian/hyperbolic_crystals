"""
    ./py/post_process_open_boundary_projector.py

    Author: Fabian R. Lux
    Date:   2023-03-22

    Plot eigenvalues of projector with open boundary conditions
"""

import iomodule
import numpy as np
import matplotlib.pyplot as plt

prl_figure_width = 1.5*(3 + 3/8)  # inch
golden_ratio = 1.61803398875


access_point = iomodule.get_access_point()
p = access_point['p']
q = access_point['q']
N = abs(access_point['N'])

eig = np.sort(np.load("./out/projector_{}_{}_open_{}.npy".format(p, q, N)))

plt.plot(eig, 'o')

plt.tight_layout()

plt.savefig("./out/projector_open.png", dpi=300)

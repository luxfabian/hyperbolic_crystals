import numpy as np

from hyperbolic_disk import triangle_area

import matplotlib.pyplot as plt

n_samples = 1000000

def random_z():

    z = (2*np.random.rand()-1) + 1j * (2*np.random.rand()-1)

    while( not abs(z)<1):
        z = (2*np.random.rand()-1) + 1j * (2*np.random.rand()-1)

    return z


def random_triangle():

    z1 = random_z()
    z2 = random_z()
    z3 = random_z()

    return triangle_area(z1,z2,z3)

areas = np.array( [ random_triangle() for _ in range(n_samples)] )

plt.hist(areas, bins=200)
plt.savefig("./out/area_histogram.png", dpi=300)



# def random_triangle():

#     z1 = random_z()
#     z2 = random_z()
#     z3 = random_z()

#     return triangle_area(z1,z2,z3), triangle_area(z2,z3,z1), triangle_area(z3,z1,z2),triangle_area(z2,z1,z3),triangle_area(z1,z3,z2),triangle_area(z3,z2,z1)

# print(random_triangle())
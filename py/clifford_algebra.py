"""
    ./topological_abelianization/fuchsian.py

    Author: Fabian R. Lux
    Date:   6/27/2023
    Mail:   fabian.lux@yu.edu

    Definition of Clifford/Gamma matrices
"""

import numpy as np

# -- setup of Pauli matrices

pauli = np.zeros((4,2,2),dtype=complex)

pauli[0] = np.array(
    [[1,0],[0,1]],
    dtype=complex
)

pauli[1] = np.array(
    [[0,1],[1,0]],
    dtype=complex
)

pauli[2] = np.array(
    [[0,-1j],[1j,0]],
    dtype=complex
)

pauli[3] = np.array(
    [[1,0],[0,-1]],
    dtype=complex
)

# -- Definition of gamma matrices via tensor product

gamma = np.zeros( (5,4,4), dtype=complex)

gamma[0] = -np.kron(pauli[2], pauli[0])
gamma[1] = np.kron(pauli[1], pauli[1])
gamma[2] = np.kron(pauli[1], pauli[2])
gamma[3] = np.kron(pauli[1], pauli[3])
gamma[4] = np.kron(pauli[3], pauli[0])

# gamma[0] = np.kron(pauli[3], pauli[0])
# gamma[1] = np.kron(pauli[1], pauli[1])
# gamma[2] = np.kron(pauli[1], pauli[2])
# gamma[3] = np.kron(pauli[1], pauli[3])
# gamma[4] = -np.linalg.multi_dot( [gamma[0], gamma[1], gamma[2], gamma[3]])


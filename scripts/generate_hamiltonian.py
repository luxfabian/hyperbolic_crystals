from scipy.sparse import lil_matrix
from iomodule import get_access_point

access_point = get_access_point()

generators = [ "A", "iA", "B", "iB" ]
H = lil_matrix((d,d))
for g in generators:
    reg_fname = fname_prefix + "_" + g + ".reg" 
    reg_data  = np.loadtxt(reg_fname, dtype=int, delimiter=" ")

    for [i,j] in reg_data:
        H[i,j] += 1.0


H_fname = fname_prefix + ".hamiltonian"

with open(H_fname, 'w') as H_file:
  for i in range(d):
    for j in H.rows[i]:
      line = str(i) + " " + str(j) + " " + str(H[i,j]) + "\n"
      H_file.write(line)

















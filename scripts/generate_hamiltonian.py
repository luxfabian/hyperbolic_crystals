import os
import re
import numpy as np

from scipy.sparse import lil_matrix

HYPERBOLIC_BUILD = os.getenv('HYPERBOLIC_BUILD')

#######################################################################
# Read group_specs.inp                                                #
#######################################################################

group_specs_fname = HYPERBOLIC_BUILD + "/bin/group_specs.inp"
group_specs_file = open(group_specs_fname, 'r')
group_specs_lines = group_specs_file.readlines()
group_specs_lines = [ int(re.sub("[^0-9]", "", line)) for line in group_specs_lines]
p,q,N = group_specs_lines
print("Group specs: ", p, q, N) 
group_specs_file.close()

periodic_boundary = not (N<0)

fname_prefix = HYPERBOLIC_BUILD + "/bin/{" + str(p) + "," + str(q) + "}_"
if(periodic_boundary):
    fname_prefix += "modulo_"+str(N)
else:
    fname_prefix += "open_"+str(-N) 

#######################################################################
# Read spectral information                                           #
#######################################################################

spec_info_fname = fname_prefix + ".info"
spec_info_file = open(spec_info_fname, 'r')
spec_info_lines = spec_info_file.readlines()
spec_info_lines = [ re.sub("[^0-9]", "", line) for line in spec_info_lines]
d = int(spec_info_lines[0])
print("Hilber space dimension: ", d) 
spec_info_file.close()

#######################################################################
# Read regular represntation                                          #
#######################################################################

generators = [ "A", "iA", "B", "iB" ]
H = lil_matrix((d,d))
for g in generators:
    reg_fname = fname_prefix + "_" + g + ".reg" 
    reg_data  = np.loadtxt(reg_fname, dtype=int, delimiter=" ")

    for [i,j] in reg_data:
       H[i,j] += 1.0

for i in range(d):
  for j in H.rows[i]:
    print(i,j, H[i,j])






















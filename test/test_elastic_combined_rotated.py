import sys
import os
# Add the root directory to sys.path
script_dir = os.path.dirname(os.path.abspath(__file__))        # test/
root_dir = os.path.abspath(os.path.join(script_dir, '..'))     # project_root/
sys.path.append(root_dir)

from helper_apdl import APDL_writer

import pprint
import unittest
import numpy as np
import subprocess
import logging

#from mitc.mitc3 import incremental_solution

from mitc.mitc3 import *
from mitc.helper_vtk import write_vtk

model = Model()

# Add material
model.add_material("steel", E=2.0e11, nu=0.3)


angle = 10

a = angle * 3.1415926/180.0
# Add nodes
n0 = model.add_node(0, 0, 0)
n1 = model.add_node(np.cos(a), 0 , np.sin(a))
n2 = model.add_node(0, 1, 0)

model.add_element([n0, n1, n2], "steel", thickness=0.1)


print("Shear Interpolation Type: ",model.elements[0].mitc_type.value)  # "APDL_OCT"

# Create solver
solver = Solver(model)

# Apply boundary conditions
fixed_dofs = [0, 1, 2, 3, 4, 5, 12, 13, 14, 15, 16, 17]
solver.apply_boundary_conditions(fixed_dofs)

# Apply load
comp = 200000
forces = np.zeros(len(model.nodes) * 6)
forces[6*1: 6*1 + 3] = comp  # Force at node 1 in z-direction

print ("forces: ", forces)

# Solve
U = solver.solve(forces, num_increments=1)

disp = np.array(U).reshape((len(model.nodes), 6))


print("Displacements:\n", disp)

# ~ Displacements:
 # ~ [[ 0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00
   # ~ 0.00000000e+00  0.00000000e+00]
 # ~ [-7.61418701e-04  5.20000000e-05  4.43963744e-03 
 # ~  -1.60442539e-04 -8.85786260e-03  5.09129278e-05]
 # ~ [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00
   # ~ 0.00000000e+00  0.00000000e+00]]

# ~ Available responses: []
# ~ Ke []
# ~ Available responses: []
# ~ Element stiffness matrix: []
# ~ DispX -0.0008183482386830239
# ~ DispX 0.00044509276611413563
# ~ DispX 0.0059292599849960406
# ~ RotX -0.0018923139418104721
# ~ RotY -0.00973391494553094
# ~ RotZ 0.00048078886010730055
# ~ Stress result: []
# ~ Available responses: []
# ~ Process 0 Terminating


### HOURGLASS COEFF 0.05

# ~ Displacements:
 # ~ [[ 0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00
   # ~ 0.00000000e+00  0.00000000e+00]
 # ~ [-7.61418701e-04  5.19943700e-05  4.43963744e-03 
   # ~ -1.46904499e-04 -8.85786260e-03 -2.58651138e-05]
 # ~ [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00
   # ~ 0.00000000e+00  0.00000000e+00]]



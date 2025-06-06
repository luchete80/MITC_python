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

# Add nodes
n0 = model.add_node(0, 0, 0)
n1 = model.add_node(1, 0, 0)
n2 = model.add_node(0, 1, 0)

model.add_element([n0, n1, n2], "steel", thickness=0.1)

model.elements[0].mitc_type = MITCShearType.TYPE2


print("Shear Interpolation Type: ",model.elements[0].mitc_type.value)  # "APDL_OCT"

# Create solver
solver = Solver(model)

# Apply boundary conditions
fixed_dofs = [0, 1, 2, 3, 4, 5, 12, 13, 14, 15, 16, 17]
solver.apply_boundary_conditions(fixed_dofs)

# Apply load
forces = np.zeros(len(model.nodes) * 6)
forces[6*1 + 0] = 200000.0  # Force at node 1 in z-direction
forces[6*1 + 1] = 200000.0  # Force at node 1 in z-direction
forces[6*1 + 2] = 200000.0  # Force at node 1 in z-direction

# Solve
U = solver.solve(forces, num_increments=1)
disp = np.array(U).reshape((len(model.nodes), 6))

print("Displacements:\n", disp)

#OPENSEES (IM ON LINUX NOW)

# ~ COMBINED LOAD (Comp 20000)
# ~ Element stiffness matrix: []
# ~ DispX 2.1076923076923078e-05
# ~ DispY 4.246153846153846e-05
# ~ DispZ 0.0007373747460087081

# ~ RotX -0.0002194484760522496
# ~ RotY -0.0011999999999999997
# ~ RotZ 7.630769230769231e-05




#RESULTS HERE (MITC)
 # ~ [[0.00000000e+00 0.00000000e+00 0.00000000e+00 0.00000000e+00
  # ~ 0.00000000e+00 0.00000000e+00]
 # ~ [1.82000000e-05 5.20000000e-05 5.52240000e-03 
 #ROTATIONS      7.61135063e-21  1.09200000e-02 0.00000000e+00]
 # ~ [0.00000000e+00 0.00000000e+00 0.00000000e+00 0.00000000e+00
  # ~ 0.00000000e+00 0.00000000e+00]]

# ~ Displacements:
 # ~ [[ 0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00
   # ~ 0.00000000e+00  0.00000000e+00]
 # ~ [ 1.82000000e-05  5.20000000e-05  5.52240000e-03 -1.49253829e-20
  # ~ -1.09200000e-02  0.00000000e+00]
 # ~ [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00
   # ~ 0.00000000e+00  0.00000000e+00]]


### -----------------------------------------------------------------
#OCT STYLE WITH CROSS:

# ~ Displacements:
 # ~ [[ 0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00
   # ~ 0.00000000e+00  0.00000000e+00]
 # ~ [ 1.82000000e-05  5.20000000e-05  5.55304833e-03 
 # ~ -1.83889980e-04 -1.09200000e-02  0.00000000e+00]
 # ~ [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00
   # ~ 0.00000000e+00  0.00000000e+00]]
### -----------------------------------------------------------------

#OCT STYLE WITH NO CROSS:

# ~ Shear Interpolation Type OPS
# ~ Displacements:
 # ~ [[ 0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00
   # ~ 0.00000000e+00  0.00000000e+00]
 # ~ [ 1.82000000e-05  5.20000000e-05  5.52240000e-03 
 #-1.49253829e-20 -1.09200000e-02  0.00000000e+00]
 # ~ [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00
   # ~ 0.00000000e+00  0.00000000e+00]]



####################################
### NEW - OCT WITH CROSS + DRILLING WITH HG
# ~ DispX 2.1076923076923078e-05
# ~ DispY 4.246153846153846e-05
# ~ DispZ 0.0007373747460087081

# ~ RotX -0.0002194484760522496
# ~ RotY -0.0011999999999999997
# ~ RotZ 7.630769230769231e-05
# ~ Displacements:
 # ~ [[ 0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00
   # ~ 0.00000000e+00  0.00000000e+00]
 # ~ [ 1.82000000e-05  5.19943686e-05  5.55304833e-03 
  # ~ -1.83889980e-04 -1.09200000e-02 -1.87712252e-08]
 # ~ [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00
   # ~ 0.00000000e+00  0.00000000e+00]]


 
### NEW - OCT WITH CROSS + DRILLING WITH NO HG - CHECK Z ROTATION 
 # Displacements:
 # [[ 0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00
   # 0.00000000e+00  0.00000000e+00]
 # [ 1.82000000e-05  5.20000000e-05  5.55304833e-03 -1.83889980e-04
  # -1.09200000e-02 -7.80000000e-05]
 # [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00
   # 0.00000000e+00  0.00000000e+00]]

#### APDL_O
    # NODE       UX           UY           UZ           USUM  
       # 1   0.0000       0.0000       0.0000       0.0000     
       # 2  0.18200E-004 0.51982E-004 0.55224E-002 0.55227E-002
       # 3   0.0000       0.0000       0.0000       0.0000     

 # MAXIMUM ABSOLUTE VALUES
  # THE FOLLOWING DEGREE OF FREEDOM RESULTS ARE IN THE GLOBAL COORDINATE SYSTEM  
 
    # NODE       ROTX         ROTY         ROTZ         RSUM  
       # 1   0.0000       0.0000       0.0000       0.0000     
       # 2 -0.21373E-016-0.10920E-001 0.52065E-004 0.10920E-001
       # 3   0.0000       0.0000       0.0000       0.0000     

 # MAXIMUM ABSOLUTE VALUES
 # NODE          0            2    
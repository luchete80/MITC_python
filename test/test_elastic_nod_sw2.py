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

dofxnod = 6

#from mitc.mitc3 import incremental_solution

from mitc.mitc3 import *
from mitc.helper_vtk import write_vtk

model = Model()

# Add material
model.add_material("steel", E=2.0e11, nu=0.3)

# Add nodes
n0 = model.add_node( 0, 1, 0)
n1 = model.add_node( 0, 0, 0)
n2 = model.add_node( 1, 0, 0)


model.add_element([n0, n1, n2], "steel", thickness=0.1)

# Create solver
solver = Solver(model)

# Apply boundary conditions
fixed_dofs = [0,1,2,3,4,5,  6,7,8,9,10,11]
solver.apply_boundary_conditions(fixed_dofs)

# Apply load
forces = np.zeros(len(model.nodes) * 6)
forces[dofxnod*2 + 2] = -50  # Force at node 1 in z-direction

# Solve
U = solver.solve(forces, num_increments=1)

# Post-process
PostProcessor.write_vtk(model, "mitc3_results.vtk")

disp = np.array(U).reshape((len(model.nodes), 6))

print("Displacements:\n", disp)






# ~ F[dofxnod*2+2]  = 0.0  # Force at node 1 in the z-direction
# ~ F[dofxnod*2+0]  = 50.0  # Force at node 1 in the z-direction
# ~ free_dofs = np.setdiff1d(np.arange(len(nodes) * dofxnod), fixed_dofs)

# ~ U = incremental_solution(nodes, elements,(E, nu), t, F, fixed_dofs,1)

# ~ disp = np.array(U).reshape((len(nodes), 6))

# ~ print("Displacements MEMBRANE:", U)

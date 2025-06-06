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


from mitc.mitc3_old import incremental_solution

from mitc.helper_vtk import write_vtk

# Material properties
E = 2.0e11  # Young's modulus (Pa)
nu = 0.3   # Poisson's ratio
t = 0.1   # Thickness (m)

dofxnod = 6


# ~ #BAD DISP
nodes = np.array([
    [0, 0, 0],  # Node 0
    [1, 0, 0],  # Node 1
    [1, 1, 0]   # Node 2
])


elements = [
    [0, 1, 2]
]

# Define the load vector
F = np.zeros(len(nodes) * dofxnod)  # Load vector for all DOFs

fixed_dofs = [0,1,2,3,4,5,  6,7,8,9,10,11]
F[dofxnod*2+2]  = -50  # Force at node 1 in the z-direction
free_dofs = np.setdiff1d(np.arange(len(nodes) * dofxnod), fixed_dofs)

U = incremental_solution(nodes, elements,(E, nu), t, F, fixed_dofs,1)

disp = np.array(U).reshape((len(nodes), 6))

print("Displacements:", U)


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
from helper_mesh import *

model = Model()




division_each_side = 10
tri = TriangleMesh(division=division_each_side, reversed=True)
tri.create()
tri.set_fixed_vertex(('O',))
tri.set_loaded_vertex(('A', ))
print(tri.nodes)
#
# print('Loaded node:', tri.loaded_node)
# print('Clamped nodes:', tri.fixed_nodes)
# print('Loaded dofs:', tri.loaded_dofs())
# print('Clamped dofs:', tri.fixed_dofs())
# print('Normal vector:', tri.normal)
# print('In-plane vectors:', tri.inplane)
# print('Number of elements:', len(tri.elements))
#
# meshfile = 'meshfile.vtk'
# displacements = [[0, 0, 0, 0, 0, 0] for x in tri.nodes]
# write_vtk(meshfile, tri.nodes, tri.elements, displacements)



model = Model()

# Add material
model.add_material("steel", E=2.0e11, nu=0.3)

# Add nodes
for node in tri.nodes:
  print ("Coords ",node[0],node[1],node[2])
  n = model.add_node(node[0],node[1],node[2])
  

for elem in tri.elements:
  print (elem)
  model.add_element([model.nodes[elem[0]],model.nodes[elem[1]],model.nodes[elem[2]]], "steel", thickness=0.1)
  #model.elements[0].mitc_type = MITCShearType.TYPE1


print("Shear Interpolation Type: ",model.elements[0].mitc_type.value)  # "APDL_OCT"

# Create solver
solver = Solver(model)

# Apply boundary conditions
fixed_dofs = [0, 1, 2, 3, 4, 5, 12, 13, 14, 15, 16, 17]
solver.apply_boundary_conditions(fixed_dofs)

# Apply load
forces = np.zeros(len(model.nodes) * 6)
forces[6*1 + 2] = -50  # Force at node 1 in z-direction

# Solve
U = solver.solve(forces, num_increments=1)

# Post-process
PostProcessor.write_vtk(model, "mitc3_results.vtk")

disp = np.array(U).reshape((len(model.nodes), 6))


print("Displacements:\n", disp)

import sys
import os
# Add the root directory to sys.path
script_dir = os.path.dirname(os.path.abspath(__file__))        # test/
root_dir = os.path.abspath(os.path.join(script_dir, '..'))     # project_root/
sys.path.append(root_dir)

import pprint
import unittest
import numpy as np
import subprocess


from mitc.mitc3_old import incremental_solution
#from helper_apdl import APDL_writer

#logger = logging.getLogger(__name__)

#np.set_printoptions(precision=10, suppress=True, floatmode='fixed')

# Material properties
E = 2.0e11  # Young's modulus (Pa)
nu = 0.3  # Poisson's ratio
t = 0.1  # Thickness (m)
dofxnod = 6  # number of DOFs per node


class SingleElementTests(unittest.TestCase):

    def check_fixed_dof_displacement(self, disp, fixed_dofs):
        """
        Check if the displacements of the fixed DOFs are zero.
        """
        for dof in fixed_dofs:
            self.assertAlmostEqual(disp[dof], 0.0, places=5)


class TestElastic(SingleElementTests):

    def test_elastic_bending_behaviour(self):

        nodes = np.array([
            [0, 0, 0],  # Node 0
            [1, 0, 0],  # Node 1
            [0, 1, 0],  # Node 2
        ])

        elements = [
            [0, 1, 2]
        ]

        # Define the load vector
        F = np.zeros(len(nodes) * dofxnod)  # Load vector for all DOFs

        # Boundary conditions: fix DOFs for node 0, 2
        fixed_dofs = [0, 1, 2, 3, 4, 10, 11, 12, 13, 14]
        F[7] = -50  # Force at node 1 in the z-direction
        U = incremental_solution(nodes, elements, (E, nu), t, F, fixed_dofs, 1)

        # checking the displacements: fixed nodes have all their displacements 0
        self.check_fixed_dof_displacement(U, fixed_dofs)

        # running APDL and checking the displacements
        elements_APDL = [[x+1 for x in elements[0]]]  # APDL uses 1-based indexing
        # instantiating the APDL writer
        #apdl = APDL_writer(E=E, nu=nu, t=t, nodes=nodes, elements=elements_APDL, fixed_dofs=fixed_dofs, forces=F)
        # writing the input file and the runner batch
        #apdl.write_script(filename='_input.txt', contents=apdl.create_content())
        #runnerfile = apdl.write_headless_runner()

        # run the APDL script using a subprocess call and retreive the exit code
        #result = subprocess.run(runnerfile, capture_output=True, text=True)
        #exit_code = result.returncode
        
        # read the displacements from the APDL output file
        APDL_displacements = apdl.read_displacements()
        # cleanup: deleting the files created by APDL
        apdl.cleanup()
        # turning the nodenumber x 6 matrix from the APDL run into a 1D array for direct comparison
        APDL_displacements = np.array(APDL_displacements).flatten()

        # printing for comparison
        pprint.pprint(U)
        pprint.pprint(APDL_displacements)
        print(len(U), len(APDL_displacements))
        print(U[7], APDL_displacements[7])
        self.assertEqual(U[7], APDL_displacements[7])  # this will fail


if __name__ == '__main__':
    unittest.main()

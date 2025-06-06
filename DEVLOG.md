# MITC_Shells
20250414 - Added curvature terms to MITC3 with covariant base
20250416 - Joined Assembly and example to curved shell
         - Fixed derivative names
         - Fixed several syntax and argument things to adapt curvature
         - Fixed bending and shear integration
20250423 - Working bending invariant (and with curvature). Assembly from 6x6 to 9x9 (global disp)
         - Added first xara examples.
         - Added Xara MITC3 C++ formulation to be migrated diretly on python.
         - Added first APDL example (t=0.1 combined).
20250424 - Added APDL second example (t=0.01, bending).
         - Working first xara single element examples.
         - Correcting bending stifness
20250425 - Reordering things. Added first Xara results.   
         - Implemented OCT al local configuration. 
           RESULTS ARE NOW OF THE SAME ORDER OF MAGNITUDE!
         - Continue working on Rotated definitve version.
20250426 - Begining to compile xara entierly since it does not return stiffness matrix.
         - NEWS: TYING POINTS: DISPLACEMENTS AND ROTATIONS! GIVES THE SAME OF APDL 1.38e-6.
         - Added result in tex format.
20250429 - Missed Thickness in local verification. Modified.
         - OK MEMBRANE Check of LOCAL version, 1 elem patch test. 
20250501 - Adding new 18x18 assembly in own formulation, to add to openssees eventually.
20250502 - Working with global, adding all matrices in local.
         - Getting first results. Global matrix is weaken.
20250505 - Fixed determinant in global version
         - Compiling global. Added tcl & Eigen dirs
         - Verifying Membrane behavior in global shells.
         - Fixed bending (added local cartesian e1,e2,e3) and shear in global coords.
20250506 - Created Gabor dev branch
         - Addes Abaqus transverse bendin test. 
         - Revisiting shear calc. 
         - Working version of GLOBAL mitc3.
20250508 - Beginning to add stresses calculation (in both local and global)
         - Beginning to add oop modular code for large systems
         - Verified Bending Stresses in benchmark.
         - Implemented and verifyied also in global
         - Added model class.
20250509 - Add C++ opensees example (with build info). 
         - Making example of 45 deg rotated shell. 
         - Fixed local covariant derivatives (error appears in 45deg rotated shell test case).
         - Now all matrix norms are the same at rotating them. 
20250511 - Find wrong shear in global version.
         - Found problem. If not added stab termn, this gives wrong result because of drilling dof.
20250512 - Adding drilling dof. Removed shear diag stab.
         - Inverting local to global.
         - FIXED GLOBAL FORMULATION
         - FIXED shear factor in global.
20250517 - Changed cartesian derivatives to a xara style. problem is now sw2 and sw3 it no giving the same as sw1
         - ADDED ANSYS APDL results for sw3. CHECK THAT, AS EXPECTED, RESULTS ARE NOT THE SAME THAT ORIGINAL.
         - FIXED LOCAL (COVARIANT) DERIVATIVES MATRIX FOR MEMBRANE
20250520 - Adding OpenSees version of MITC Shear Stiffness on local coords
         - Verifying opensees results.
         - Added Scale distortion.
20250521 - Fixed several matrix derivatives.
         - Changed deriv calc from classic inverse jacobian to B/C
         - Working on OOP     
20250522 - Working new format    
20250523 - IMPORTANT: FIXED GLOBAL ROTATION MATRIX. 

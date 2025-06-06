# MITC Shells

MITC Shells is a Finite Element Formulation about MITC Triangular (and then quads),
Mixed Interpolation Tensor Components (Dvorkin & Bathe 1984) Shell Formulation.

Results are Compared with APDL and OpenSees.

##Own Solver - RESULTS ARE THE SAME OF ANSYS.

  ### Bending with Transverse Load
  Check own Solver MITC3 Shear locking algorithm!!!! (LOCAL Dir).
  Gives the same results than ANSYS 181 element.<br>
  Still pendant working on rotation to global. 

  [^1]:- 50N Transverse LOads
  [^2]: - Nodes (0, 0), (1,0),(0,1)

  >  Norm Km 388520209443.158  
  >  Norm Kb 32376684.120263178  
  >  Norm Ks 30732065562.12357  
  >  Displacements: [UX,UY, UZ , RX, RY ] (LOCAL)  
  >                 [ 0.0000000e+00  0.0000000e+00  0.0000000e+00  0.0000000e+00 0.0000000e+00]  
  >                 [0.0000000e+00  0.0000000e+00 -1.3805381e-06   3.0952381e-08 -2.7300000e-06 ]  
  >                 [0.0000000e+00  0.0000000e+00 0.0000000e+00  0.0000000e+00  0.0000000e+00]  

## Membrane
>  Membrane, 50N Displacements:  
>  [0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00]
 4.55e-09 0.00e+00 0.00e+00
 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00 0.00e+00]
 
## APDL Results
  ### Bending with Transverse Load
  > LOAD STEP=     1  SUBSTEP=     1                                             
  >  TIME=    1.0000      LOAD CASE=   0                                         
 
  > THE FOLLOWING DEGREE OF FREEDOM RESULTS ARE IN THE GLOBAL COORDINATE SYSTEM  
  >  NODE       UX           UY           UZ           USUM  
  >       1   0.0000       0.0000       0.0000       0.0000     
  >       2   0.0000       0.0000     -0.13806E-005 0.13806E-005
  >       3   0.0000       0.0000       0.0000       0.0000    

  >    NODE       ROTX         ROTY         ROTZ         RSUM  
  >       1   0.0000       0.0000       0.0000       0.0000     
  >       2  0.53433E-020 0.27300E-005  0.0000      0.27300E-005  
  >       3   0.0000       0.0000       0.0000       0.0000   
  

  ### Membrane
  >      NODE       UX           UY           UZ           USUM  
  >       1   0.0000       0.0000       0.0000       0.0000     
  >       2  0.45500E-008  0.0000       0.0000      0.45500E-008  
  >       3   0.0000       0.0000       0.0000       0.0000     

 
## GLOBAL MITC3 Implementation



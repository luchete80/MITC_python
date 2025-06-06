import numpy as np


import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.ERROR)

# Degree of freedom per node (u, v, w, theta_x, theta_y, theta_z)
dofxnod = 6  

def shape_functions_3d_triangle(xi, eta):
    """
    Shape function evaluation and mapping for a 3-node triangle surface element in 3D.
    Returns:
        N: Shape functions
        dN_dxi, dN_deta: Derivatives of shape functions
    """
    # Shape functions
    N = np.array([
        1.0 - xi - eta,
        xi,
        eta
    ])
    
    dN_dxi  = np.array([-1.0, 1.0, 0.0])
    dN_deta = np.array([-1.0, 0.0, 1.0])

    return N, dN_dxi, dN_deta


def compute_local_derivatives(dN_dxi, dN_deta, a1, a2,a3):
    # Create local orthonormal basis
    e1 = a1/np.linalg.norm(a1)
    e3 = a3/np.linalg.norm(a3)
    e2 = np.cross(e3, e1)  # Ensures orthonormal system
    
    #T = setup_local_coordinate_system
    # Project a1 and a2 onto tangent plane
    a1_proj = np.array([np.dot(a1,e1), np.dot(a1,e2)])
    a2_proj = np.array([np.dot(a2,e1), np.dot(a2,e2)])
    
    # Create 2x2 Jacobian and invert
    J = np.column_stack((a1_proj, a2_proj))
    invJ = np.linalg.inv(J)
    
    # Compute derivatives in local coordinates
    dN_dg1 = invJ[0,0]*dN_dxi + invJ[0,1]*dN_deta
    dN_dg2 = invJ[1,0]*dN_dxi + invJ[1,1]*dN_deta
    
    return dN_dg1, dN_dg2
   

def setup_local_coordinate_system(a1, a2, a3, tol=1e-10):
    """
    Creates an orthonormal basis where:
    - e1 ≈ direction of a1 (normalized)
    - e2 ≈ direction of a2 (orthogonalized to e1)
    - e3 = normal vector (a1 × a2, normalized)
    """
    # Normalize normal vector
    e3 = a3 / (np.linalg.norm(a3) + tol)
    
    # First local axis (e1) along a1 (normalized)
    e1 = a1 / (np.linalg.norm(a1) + tol)
    
    # Second local axis (e2) orthogonal to e1 but as close to a2 as possible
    e2 = a2 - (a2.dot(e1)) * e1  # Subtract a2's component along e1
    e2 = e2 / (np.linalg.norm(e2) + tol)
    
    # Ensure exact orthogonality
    e3 = np.cross(e1, e2)
    e3 = e3 / (np.linalg.norm(e3) + tol)
    
    #print ("e1,e2,e3",e1,e2,e3)
    return np.vstack((e1, e2, e3)).T  # 3x3 transformation matrix


def create_transformation_matrix(a1, a2, a3):
    """Create 18×18 transformation matrix from local to global coordinates"""
    # Create local basis
    T_small = setup_local_coordinate_system(a1, a2, a3)
    
    # Create block-diagonal transformation matrix (18×18)
    T = np.zeros((18, 18))
    for i in range(3):  # Loop over 3 nodes
        # Apply T_small to translational DOFs (DOFs 0-2 per node)
        T[6*i:6*i+3, 6*i:6*i+3] = T_small
        
        # Apply T_small to rotational DOFs (DOFs 3-5 per node)
        T[6*i+3:6*i+6, 6*i+3:6*i+6] = T_small
    

    #print ("T",T)

    #print ("TtxT",np.allclose(T.T @ T, np.eye(18)))
    return T


def compute_local_membrane_stiffness(dN_dx, dN_dy, nodes, material_properties, thickness, detJ, weight):
    E, nu = material_properties

    #print ("dN_dx,dN_dy",dN_dx,dN_dy)
    # Local membrane constitutive matrix
    D_m = E / (1 - nu**2) * np.array([
        [1, nu, 0],
        [nu, 1, 0],
        [0, 0, (1 - nu)/2]
    ])
    
    B_m = np.zeros((3, 6))  # 3 strains × 6 membrane DOFs
    #print ("a1, a2, a3",a1,a2,a3)
    for i in range(3):
        B_m[0, 2*i]   = dN_dx[i]  # ε_xx = ∂u/∂x
        B_m[1, 2*i+1] = dN_dy[i]  # ε_yy = ∂v/∂y
        B_m[2, 2*i]   = dN_dy[i]  # γ_xy = ∂u/∂y + ∂v/∂x
        B_m[2, 2*i+1] = dN_dx[i]
    
    #print ("MEMBRANE Bm", B_m)
    
    K_local = np.zeros((18, 18))
    active_dofs = [0, 1, 6, 7, 12, 13]  # u,v DOFs
    
    K_condensed = B_m.T @ D_m @ B_m * detJ * weight * thickness
    for i, dof_i in enumerate(active_dofs):
        for j, dof_j in enumerate(active_dofs):
            K_local[dof_i, dof_j] = K_condensed[i, j]
    
    return K_local



def compute_local_bending_stiffness(dN_dx, dN_dy, nodes, material_properties, thickness, detJ, weight):
    """Build complete 18×18 bending stiffness in LOCAL coordinates"""
    E, nu = material_properties
    
    # Bending constitutive matrix (local coordinates)
    D_b = E * thickness**3 / (12 * (1 - nu**2)) * np.array([
        [1,    nu,     0],
        [nu,    1,     0],
        [0,     0, (1-nu)/2]
    ])



    B_b = np.zeros((3, 18))

    for i in range(3):
        #print ("i, dNdx,dNdy",i, dN_dx[i],dN_dy[i])
        # Curvatures in local basis (κ_11, κ_22, κ_12)
        B_b[0, 6*i+4] = -dN_dx[i]  # κ_11 = -∂θ2/∂g1
        B_b[1, 6*i+3] = dN_dy[i]   # κ_22 = ∂θ1/∂g2
        B_b[2, 6*i+3] = dN_dx[i]   # κ_12 terms
        B_b[2, 6*i+4] = -dN_dy[i]
    
    # Compute stiffness matrix
    K_b = B_b.T @ D_b @ B_b * detJ * weight

    # Initialize B matrix (3 curvatures × 6 rotational DOFs)
    B_bc = np.zeros((3, 6))
    
    for i in range(3):
        # Compute Cartesian derivatives
        # dN_dx = invJ[0,0]*dN_dxi[i] + invJ[0,1]*dN_deta[i]
        # dN_dy = invJ[1,0]*dN_dxi[i] + invJ[1,1]*dN_deta[i]
        
        # Standard curvature definitions:
        B_bc[0, 2*i+1] =  -dN_dx[i]  # κ_xx = -∂θy/∂x
        B_bc[1, 2*i] = dN_dy[i]     # κ_yy = ∂θx/∂y
        B_bc[2, 2*i] = dN_dx[i]     # κ_xy terms
        B_bc[2, 2*i+1] = -dN_dy[i]
    
    # Compute stiffness matrix
    K_bc = B_bc.T @ D_b @ B_bc * detJ * weight
    #print ("B_bc", B_bc)
    #print ("K_bc", K_bc)
    return K_b#, K_bc


# using all nodes in the shape function (for i in range(3)) — which contributes to all 9 DOFs per node in the shear B matrix.
#In the local version, you're using only 2 nodes at a time (node_a and node_b), giving nonzero entries for only those two nodes per tying point.
#In the global version, each tying point contributes to all three nodes — and the total stiffness matrix accumulates full coupling between them.
#In the local version, each tying point only contributes to the stiffness of the associated edge, i.e., partial matrix blocks.
#Thus, to obtain the same result, the contributions from each tying point in the local version must be counted more than once to build the full 3x3 block matrix.

def compute_local_shear_stiffness_mitc3(nodes, material_properties, thickness):
    E, nu = material_properties
    G = E / (2 * (1 + nu))
    k_shear = 5/6  # Shear correction factor

    # Covariant basis vectors (a1, a2)
    # FOR LINEAR TRIANGLES IS THE SAME AS a
    a1 = nodes[1] - nodes[0] # ∂x/∂ξ
    a2 = nodes[2] - nodes[0] # ∂x/∂η
    a3 = np.cross(a1, a2)
    #print("a1,a2 FOR SHEAR",a1,a2)
    A = 0.5 * np.linalg.norm(a3)
    detJ = 2 * A  # det(J) = 2A for triangles

    # In-plane Jacobian (2x2 metric tensor)
    J = np.array([
        [np.dot(a1, a1), np.dot(a1, a2)],  # a1·a1, a1·a2
        [np.dot(a2, a1), np.dot(a2, a2)]   # a2·a1, a2·a2
    ])
    invJ = np.linalg.inv(J)  # Now works (2x2 matrix)

    # Tying points in natural coordinates (ξ, η)
    tying_points = [(0.5, 0.0), (0.5, 0.5), (0.0, 0.5)]
    edge_nodes = [(0, 1), (1, 2), (2, 0)]

    K_shear = np.zeros((18, 18))

#    for (xi, eta), (node_a, node_b) in zip(tying_points, edge_nodes):
    for (xi, eta) in (tying_points):
        # ~ # Natural derivatives at tying point
        dN_dxi = np.array([-1, 1, 0])  # ∂N/∂ξ
        dN_deta = np.array([-1, 0, 1]) # ∂N/∂η
        #print ("SHEAR dN_dxi, dN_deta",dN_dxi, dN_deta)
        # Covariant B-matrix (γ_ξz, γ_ηz)
        B_cov = np.zeros((2, 18))
        for i in range(3):
            B_cov[0, 6*i + 2] = dN_dxi[i]   # ∂w/∂ξ
            B_cov[0, 6*i + 4] = -0.5        # θ_y (MITC3 interpolation)
            B_cov[1, 6*i + 2] = dN_deta[i]  # ∂w/∂η
            B_cov[1, 6*i + 3] = 0.5         # θ_x

        # ~ N = np.array([1 - xi - eta, xi, eta])
        # ~ # Local derivatives of shape functions
        # ~ dN_dg1 = invJ[0,0]*dN_dxi + invJ[0,1]*dN_deta
        # ~ dN_dg2 = invJ[1,0]*dN_dxi + invJ[1,1]*dN_deta

        # ~ for i in range(3):
            # ~ B_cov[0, 6*i + 2] = dN_dg1[i]    # ∂w/∂g1
            # ~ B_cov[0, 6*i + 4] = -N[i]        # θ_y interpolated at tying point
            # ~ B_cov[1, 6*i + 2] = dN_dg2[i]    # ∂w/∂g2
            # ~ B_cov[1, 6*i + 3] = N[i]         # θ_x interpolated at tying point
            
        # Transform to local Cartesian (γ_xz, γ_yz) using inv(J)
        B_shear = invJ @ B_cov  # B_shear = J⁻¹ B_cov
        #print ("Bshear orthogonal", B_shear)
        #print ("J",J, "invJ",invJ)
        #B_shear = J @ B_cov  # Replace invJ with J (direct Jacobian)

        # Stiffness contribution (integrate with edge length weight)
        #edge_vec = nodes[node_b] - nodes[node_a]
        #L = np.linalg.norm(edge_vec)
        #print ("(detJ / 6)",(detJ / 6))
        #K_shear += (B_shear.T @ B_shear) * (k_shear * G * thickness) * L
        # Changed: Use area integration (detJ/6) instead of edge length L
        K_shear += (B_shear.T @ B_shear) * (k_shear * G * thickness) * (detJ / 6)
        #print("K_shear",K_shear)
        #K_shear += (B_shear.T @ B_shear) * (k_shear * G * thickness) * (L/2)
        
    #print("K_shear",K_shear)
    return K_shear

#
def compute_local_shear_stiffness_mitc3_OPS(nodes, dNdx,dNdy, material_properties, thickness):
    E, nu = material_properties
    G = E / (2 * (1 + nu))
    k_shear = 5/6  # Shear correction factor

    # Covariant basis vectors
    a1 = nodes[1] - nodes[0]
    a2 = nodes[2] - nodes[0]
    a3 = np.cross(a1, a2)
    print ("a1,a2,a3")
    A = 0.5 * np.linalg.norm(a3)
    detJ = 2 * A  # det(J) = 2A for triangles

    # Compute local Cartesian derivatives (same as bending)
    dN_dxi = np.array([-1, 1, 0])
    dN_deta = np.array([-1, 0, 1])


    # Tying points in natural coordinates (ξ, η)
    tying_points = [(0.5, 0.0), (0.5, 0.5), (0.0, 0.5)]
    K_shear   = np.zeros((18, 18))
    K_shear_c = np.zeros((9, 9)) #Condensed
    
    T = setup_local_coordinate_system(a1,a2,a3)
    
    translated_nodes = nodes - nodes[0]  # shift node 0 to origin
    
    #loc_nodes = np.dot(T, translated_nodes.T).T
    loc_nodes =  translated_nodes @ T
    #print ("T",T)
    #print ("global nodes",nodes)
    #print ("trans nodes",translated_nodes)
    print ("loc_nodes (node0 always 0,0)",loc_nodes)
    
    
    
    # tying_edges = [(0, 1), (1, 2), (2, 0)]  # Which edge each tying point belongs to

    # for (xi, eta), (i, j) in zip(tying_points, tying_edges):
        # N = np.zeros(3)
        # N[i] = 0.5
        # N[j] = 0.5
    cross = True 
    hardc = True

    # Compute edge vectors
    x21 = loc_nodes[1][0] - loc_nodes[0][0]
    y21 = loc_nodes[1][1] - loc_nodes[0][1]
    x32 = loc_nodes[2][0] - loc_nodes[1][0]
    y32 = loc_nodes[2][1] - loc_nodes[1][1]
    x13 = loc_nodes[0][0] - loc_nodes[2][0]
    y13 = loc_nodes[0][1] - loc_nodes[2][1]
    
    x31 = -x13
    y31 = -y13
    # Compute the angles
    phi1 = np.arctan2(y21, x21)
    phi2 = 0.5 * np.pi - np.arctan2(x31, y31)

    # Create and assign BsO matrix
    BsO = np.zeros((2, 2))
    BsO[0, 0] = np.sin(phi2)
    BsO[0, 1] = -np.sin(phi1)
    BsO[1, 0] = -np.cos(phi2)
    BsO[1, 1] = np.cos(phi1)

    # Create and assign BsC matrix
    BsC = np.zeros((2, 2))
    A2 = 2.0 * A
    BsC[0, 0] = np.sqrt(x13**2 + y13**2) / A2
    BsC[1, 1] = np.sqrt(x21**2 + y21**2) / A2

    #############
      
    for xi, eta in tying_points:
        # Shape functions at tying point
        N = np.array([1 - xi - eta, xi, eta])

        print ("Cross 0 xi eta",xi,eta, (y21 + y32 * xi) / 2.0, (y21 + y13 * xi) / 2.0, -(x21 * xi) / 2.0)
        print ("Cross 1 xi eta",xi,eta, (y13 + y32 * eta) / 2.0, (y13 * eta) / 2.0, (y13 + y21 * eta) / 2.0)    

        B_s = np.zeros((2, 9))
        B_s_OCT = np.zeros((2, 9))
        
        # Row 0: shear strain γ_xz
        B_s[0, 0] = -1.0  # θ_x at node 1 
        B_s[0, 1] = -(y21 + y32 * xi) / 2.0
        B_s[0, 2] = (x21 + x32 * xi) / 2.0        
        B_s[0, 3] = 1.0  # θ_x at node 2
        B_s[0, 4] = -(y21 + y13 * xi) / 2.0
        B_s[0, 5] = (x21 + x13 * xi) / 2.0
        #if (not hardc):B_s[0, 6] = dN_dx[2]
        B_s[0, 7] = -(y21 * xi) / 2.0
        B_s[0, 8] = (x21 * xi) / 2.0
        #############################
        # Row 1: shear strain γ_yz
        B_s[1, 0] = -1.0  # θ_y at node 1
        B_s[1, 1] = (y13 + y32 * eta) / 2.0
        B_s[1, 2] = -(x13 + x32 * eta) / 2.0
        B_s[1, 4] = (y13 * eta) / 2.0
        B_s[1, 5] = -(x13 * eta) / 2.0
        B_s[1, 6] = 1.0  # θ_y at node 3
        B_s[1, 7] = (y13 + y21 * eta) / 2.0
        B_s[1, 8] = -(x13 + x21 * eta) / 2.0
        
        #print ("BS",B_s)
        # Note:" You may need to adjust indexing if your DOFs per node differ.
        
        B_s_OCT = BsO @ BsC @ B_s

        # Continue with shear stiffness assembly as before
        #K_shear_c += (B_s.T @ B_s) * (k_shear * G * thickness) * (detJ / 6)
        K_shear_c += (B_s_OCT.T @ B_s_OCT) * (k_shear * G * thickness) * (detJ / 6)

        active_dofs = [2,3,4,  8,9,10, 14,15,16]
      
        
        for i, dof_i in enumerate(active_dofs):
            for j, dof_j in enumerate(active_dofs):
                K_shear[dof_i, dof_j] = K_shear_c[i, j]

        return K_shear
        
###########################
####  DRILLING 

def compute_pure_drilling_components(dN_dx,dN_dy,material_properties, thickness, A2):
    """Compute only the drilling-related stiffness components"""
    # Shape functions at centroid (constant for linear triangle)
    N = np.array([1/3, 1/3, 1/3])
    
    # Initialize drilling B-matrix (1×18)
    B_drill = np.zeros((1, 18))  # Only γ_drill strain

    E, nu = material_properties
        
    # Drilling strain-displacement relationship
    for i in range(3):
        B_drill[0, 6*i] = -0.5 * dN_dy[i]    # ½(∂v/∂x) component
        B_drill[0, 6*i+1] = 0.5 * dN_dx[i]   # -½(∂u/∂y) component
        B_drill[0, 6*i+5] = -N[i]                 # -θz component
    
    # Drilling material term (scalar)
    alpha = 0.1  # Stabilization factor
    D_drill = alpha * E * thickness**3 / 12
    
    # Pure drilling stiffness matrix
    K_drill = B_drill.T @ B_drill * D_drill * A2
    
    return K_drill

def compute_element_distortion( nodes, A2):
    """Compute metrics for adaptive hourglass control."""
    # Edge lengths
    L1 = np.linalg.norm(nodes[1].coords - nodes[0].coords)
    L2 = np.linalg.norm(nodes[2].coords - nodes[1].coords)
    L3 = np.linalg.norm(nodes[0].coords - nodes[2].coords)
    
    # Aspect ratio (max edge / min edge)
    max_edge = max(L1, L2, L3)
    min_edge = min(L1, L2, L3)
    aspect_ratio = max_edge / min_edge
    
    # Skewness (deviation from ideal equilateral triangle)
    ideal_area = (np.sqrt(3)/4) * (min_edge**2)
    skewness = abs(A2 - ideal_area) / ideal_area
    
    return aspect_ratio, skewness

def compute_hourglass_stabilization(nodes, dN_dx,dN_dy, material_properties, thickness, A2):
    """ANSYS-style adaptive hourglass control."""
    aspect_ratio, skewness = self.compute_element_distortion(nodes,A2)

    E, nu = material_properties
    G = E / (2 * (1 + nu))
    
    
    # Base stabilization factor (similar to ANSYS defaults)
    beta_base = 0.01  # Default for well-shaped elements
    
    # Increase stabilization for distorted elements
    beta = beta_base * (1 + 0.5 * (aspect_ratio - 1) + 2.0 * skewness)
    
    # Limit to reasonable range (0.01 to 0.1)
    beta = np.clip(beta, 0.01, 0.1)
    
    # Hourglass stiffness (same formulation as before)
    Bhx = np.zeros(18)
    Bhy = np.zeros(18)
    for i in range(3):
        Bhx[6*i+5] = dN_dx[i]  # ∂θz/∂x
        Bhy[6*i+5] = dN_dy[i]  # ∂θz/∂y
    
    K_hg = beta * G * thickness * A2 * (
        np.outer(Bhx, Bhx) + np.outer(Bhy, Bhy))
    
    return K_hg


# ~ def compute_drill_stiffness(self):
  # ~ """Combined physical drilling and hourglass stabilization"""
  # ~ K_drill = self.compute_pure_drilling_components()
  # ~ K_hg = self.compute_hourglass_stabilization()
  # ~ return K_drill + K_hg
    
 
def compute_element_stiffness(nodes, material_properties, thickness, disp=None):
    """Compute the complete 18×18 element stiffness matrix in global coordinates"""
    # Initialize stiffness matrix
    K_global = np.zeros((18, 18))
    
    # Compute covariant basis vectors
    dN_dxi = np.array([-1, 1, 0])
    dN_deta = np.array([-1, 0, 1])
    a1, a2, a3, norm_a3 = covariant_basis_vectors_triangle(dN_dxi, dN_deta, nodes)
    detJ = norm_a3 # Jacobian determinant for area integration
    
    # Create transformation matrix
    T = create_transformation_matrix(a1, a2, a3)
    #print ("T",T)
    
    
    # With this (1-point centroid rule):
    gauss_points = [(1/3, 1/3)]  # Centroid of triangle
    weights = [0.5]              # Weight = detJ/2 (area of unit triangle)
    
    # Compute local stiffness matrices at integration points
    K_mem_local = np.zeros((18, 18))
    K_bend_local = np.zeros((18, 18))
                                                   

    ##### SCALING
    T = setup_local_coordinate_system(a1,a2,a3)    
    translated_nodes = nodes - nodes[0]  # shift node 0 to origin
    
    #loc_nodes = np.dot(T, translated_nodes.T).T
    loc_nodes =  translated_nodes @ T

    # Compute edge vectors
    x21 = loc_nodes[1][0] - loc_nodes[0][0]
    y21 = loc_nodes[1][1] - loc_nodes[0][1]
    x32 = loc_nodes[2][0] - loc_nodes[1][0]
    y32 = loc_nodes[2][1] - loc_nodes[1][1]
    x13 = loc_nodes[0][0] - loc_nodes[2][0]
    y13 = loc_nodes[0][1] - loc_nodes[2][1]
    

    A2 = x21 * y32 - x32 * y21
    #print ("A2",A2)
        
    b = np.array([y32, y13, y21])
    c = np.array([-x32, -x13, -x21])
    dN_dx = -b / A2
    dN_dy = -c / A2

    # Compute local stiffness matrices
    K_mem_local += compute_local_membrane_stiffness(dN_dx, dN_dy, nodes, 
                                                  material_properties, thickness, 
                                                  detJ, weight)
    
    K_bend_local += compute_local_bending_stiffness(dN_dx, dN_dy, nodes,
                                                  material_properties, thickness,
                                                  detJ, weight)
   
       

    K_shear_local = compute_local_shear_stiffness_mitc3_OPS(nodes, dN_dx,dN_dy,material_properties, thickness)

    K_mem_global = T @ K_mem_local @ T.T
    K_bend_global = T @ K_bend_local @ T.T
    K_shear_global = T @ K_shear_local @ T.T
    
    k_drill_local = np.zeros((18, 18))
    

    K_drill = compute_pure_drilling_components(dN_dx,dN_dy,material_properties, thickness, A2)
    K_hg    = compute_hourglass_stabilization(nodes, dN_dx,dN_dy, material_properties, thickness, A2)
    
    k_drill_local = K_drill + K_hg

    eigvals = np.linalg.eigvalsh(K_global)
    # print("Number of zero eigenvalues:", sum(np.abs(eigvals) < 1e-6))

    K_local = K_mem_local + K_bend_local + K_shear_local +k_drill_local 
    # Sum all components
    K_global = K_mem_global + K_bend_global + K_shear_global + T @ k_drill_local @ T.T

    #print("K_mem_global", K_mem_global)
    #print("K_bend_global", K_bend_global)

    logger.debug("Norm Km: {}".format(np.linalg.norm(K_mem_local)))
    logger.debug("Norm Kb: {}".format(np.linalg.norm(K_bend_local)))
    logger.debug("Norm Ks: {}".format(np.linalg.norm(K_shear_local)))
    
    #print ("K LOCAL", K_local)    
    #print ("K GLOBAL", K_global)
    
    return K_global

def assemble_global_stiffness(elements, nodes, material_properties, thickness, disp=None):
    """
    Assemble the global stiffness matrix for multiple elements.
    
    Parameters:
        elements: List of elements, where each element is a list of node indices
        nodes: Array of node coordinates
        material_properties: Tuple of (E, nu)
        thickness: Shell thickness
        disp: Optional displacement vector for nonlinear analysis
        
    Returns:
        K_global: Global stiffness matrix
    """
    num_nodes = len(nodes)
    num_dofs = num_nodes * dofxnod
    K_global = np.zeros((num_dofs, num_dofs))
    
    for element in elements:
        element_nodes = nodes[element]
        
        # Compute element stiffness matrix (18×18)
        K_element = compute_element_stiffness(element_nodes, material_properties, thickness, disp)
        
        # Assemble into global matrix
        for i, node_i in enumerate(element):
            for j, node_j in enumerate(element):
                # Get global DOF indices for this node pair
                dof_i = node_i * dofxnod
                dof_j = node_j * dofxnod
                
                # Add element contribution to global matrix
                K_global[dof_i:dof_i+dofxnod, dof_j:dof_j+dofxnod] += \
                    K_element[i*dofxnod:(i+1)*dofxnod, j*dofxnod:(j+1)*dofxnod]
    
    return K_global

def incremental_solution(nodes, elements, material_properties, thickness, F, fixed_dofs, num_increments):
    """
    Incremental solution with Newton-Raphson iterations for nonlinear analysis
    
    Parameters:
        nodes: Array of node coordinates
        elements: List of elements
        material_properties: Tuple of (E, nu)
        thickness: Shell thickness
        F: Global load vector
        fixed_dofs: List of constrained DOFs
        num_increments: Number of load increments
        
    Returns:
        U: Displacement vector
    """
    num_dofs = len(nodes) * dofxnod
    U = np.zeros(num_dofs)
    free_dofs = np.setdiff1d(np.arange(num_dofs), fixed_dofs)
    
    for inc in range(num_increments):
        # print(f"Increment {inc+1}/{num_increments}")
        
        # Compute tangent stiffness matrix
        K_global = assemble_global_stiffness(elements, nodes, material_properties, thickness, U)


                
        # Solve for incremental displacements
        F_inc = F / num_increments
        delta_U = np.zeros(num_dofs)
        delta_U[free_dofs] = np.linalg.solve(K_global[np.ix_(free_dofs, free_dofs)], F_inc[free_dofs])
        
        # Update displacements
        U += delta_U
        
    
    for el in elements:
      u_elem = get_u_elem(U,el) 
      stress = calculate_element_stresses(nodes, material_properties, thickness, u_elem)
      # print ("stress membrane: ", stress['membrane'])
      # print ("stress bending: ", stress['bending'])
    
    return U

####################### ADDED FOR STRESSES
def get_u_elem(U, element):
    u_elem = []
    for node_id in element:
        start = node_id * dofxnod
        u_elem.extend(U[start:start + dofxnod])
    return np.array(u_elem)

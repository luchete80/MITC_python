import numpy as np
from mitc.helper_vtk import write_vtk

from enum import Enum

#TYPES OF SHEAR MITC
class MITCShearType(Enum):
    TYPE1 = "APDL_OCT"
    TYPE2 = "OPS_OCT"
    TYPE3 = "OPS_OCT_NOCROSS"

class MITCHourglassType(Enum):
    TYPE1 = "APDL" #DISTORTION DEPENDENT ENHANCED
    TYPE2 = "CONSTANT" #OpenSees

class Node:
    def __init__(self, node_id, x, y, z):
        self.id = node_id
        self.coords = np.array([x, y, z])
        self.displacements = np.zeros(6)  # u, v, w, θx, θy, θz
        self.forces = np.zeros(6)
        
    def update_position(self):
        """Update position after displacement"""
        return self.coords + self.displacements[:3]
        
    def get_dof_indices(self):
        """Return global DOF indices for this node"""
        return np.arange(self.id * 6, (self.id + 1) * 6)

class Material:
    def __init__(self, name, E, nu):
        self.name = name
        self.E = E
        self.nu = nu
        self.G = E / (2 * (1 + nu))
        
    def get_membrane_matrix(self, thickness):
        """Return membrane constitutive matrix"""
        factor = self.E / (1 - self.nu**2)
        return factor * np.array([
            [1, self.nu, 0],
            [self.nu, 1, 0],
            [0, 0, (1 - self.nu)/2]
        ])
        
    def get_bending_matrix(self, thickness):
        """Return bending constitutive matrix"""
        factor = self.E * thickness**3 / (12 * (1 - self.nu**2))
        return factor * np.array([
            [1, self.nu, 0],
            [self.nu, 1, 0],
            [0, 0, (1 - self.nu)/2]
        ])
        
    def get_shear_matrix(self, thickness):
        """Return shear constitutive matrix"""
        k_shear = 5/6  # Shear correction factor
        return k_shear * self.G * thickness * np.eye(2)


class Model:
    def __init__(self):
        self.nodes = []       # List of Node objects
        self.elements = []    # List of Element objects
        self.materials = {}   # Dict of Material objects
        self.solver = Solver(self) 

    def add_node(self, x, y, z):
        node = Node(len(self.nodes), x, y, z)
        self.nodes.append(node)
        return node

    def add_element(self, nodes, material_name, thickness):
        material = self.materials[material_name]
        print ("ELEMENT NODES", nodes)
        elem = MITC3Element(len(self.elements), nodes, material, thickness)
        self.elements.append(elem)
        return elem

    def solve(self, forces, fixed_dofs, num_increments=1):
        # Delegate to solver
        return self.solver.solve(self, forces, fixed_dofs, num_increments)

    def write_results(self, filename):
        # Delegate to VTK writer
        node_coords = np.array([n.coords + n.displacements[:3] for n in self.nodes])
        elements = [[n.id for n in elem.nodes] for elem in self.elements]
        displacements = np.array([n.displacements for n in self.nodes])
        write_vtk(filename, node_coords, elements, displacements)
        
    def add_material(self,Name,E, nu):
        #my_dict = {'my_list': []}
        self.materials[Name] = Material(Name, E, nu)

class Solver:
    def __init__(self, model):
        self.model = model
        self.K_global = None
        self.free_dofs = None
        self.fixed_dofs = None
        
    def assemble_global_stiffness(self):
        num_dofs = len(self.model.nodes) * 6
        self.K_global = np.zeros((num_dofs, num_dofs))
        
        for elem in self.model.elements:
            K_elem = elem.compute_stiffness_matrix()
            elem_dofs = elem.get_element_dof_indices()
            
            # Vectorized assembly
            rows, cols = np.meshgrid(elem_dofs, elem_dofs, indexing='ij')
            self.K_global[rows, cols] += K_elem
            
        return self.K_global
        
    def apply_boundary_conditions(self, fixed_dofs):
        self.fixed_dofs = np.array(fixed_dofs)
        self.free_dofs = np.setdiff1d(np.arange(len(self.model.nodes)*6), self.fixed_dofs)
        
    def solve(self, forces, num_increments=1):
        U = np.zeros(len(self.model.nodes) * 6)
        
        for inc in range(num_increments):
            print(f"Increment {inc+1}/{num_increments}")
            
            K_global = self.assemble_global_stiffness()
            F_inc = forces / num_increments
            
            delta_U = np.zeros(len(self.model.nodes) * 6)
            delta_U[self.free_dofs] = np.linalg.solve(
                K_global[np.ix_(self.free_dofs, self.free_dofs)],
                F_inc[self.free_dofs]
            )
            
            U += delta_U
            
            # Update node displacements
            for i, node in enumerate(self.model.nodes):
                node.displacements = U[i*6:(i+1)*6]
                
        return U

class PostProcessor:
    @staticmethod
    def write_vtk(model, filename):
        node_coords = np.array([n.coords + n.displacements[:3] for n in model.nodes])
        elements = [[n.id for n in elem.nodes] for elem in model.elements]
        displacements = np.array([n.displacements for n in model.nodes])
        write_vtk(filename, node_coords, elements, displacements)



class MITC3Element:
    def __init__(self, elem_id, nodes, material, thickness):
        self.id = elem_id
        self.nodes = nodes
        self.material = material
        self.thickness = thickness
        self._a1 = None  # Cached covariant basis
        self._a2 = None
        self._a3 = None
        self._T_small = None  # Cached local-to-global rotation
        self._transformation_matrix = None
        self.K = np.zeros((18, 18))
        self.mitc_type = MITCShearType.TYPE2
        self.hg_type = MITCHourglassType.TYPE1 #APDL
        self.B_m = np.zeros((3, 6))  # 3 strains × 6 membrane DOFs, FOR PRERFORMANCE COULD BE 3x 9 
        self.B_b = np.zeros((3, 6))  # BENDING
        self.B_s = np.zeros((3, 9))  # SHEAR
                        
        self.D_m = np.zeros((3,3))
        self.stresses = {
            'membrane': [],
            'bending': [],
            'shear': [],        # Stresses at upper face (z = +t/2)
            'upper': [],        # Stresses at upper face (z = +t/2)
            'lower': [],        # Stresses at lower face (z = -t/2)
            'max': [],          # Maximum principal stress
            'min': []          # Minimum principal stress
        }
        self.strains = {
            'membrane': [],
            'bending': [],
            'shear': []
        }
                
    @property
    def a1(self):
        """Lazy-load covariant basis vector a1."""
        if self._a1 is None:
            self._compute_basis_vectors()
        return self._a1

    # ~ @property
    # ~ def T_small(self):
        # ~ """Lazy-load local coordinate system."""
        # ~ if self._T_small is None:
            # ~ self._setup_local_coordinate_system()
        # ~ return self._T_small

    def create_transformation_matrix(self):
        """Create 18×18 transformation matrix"""
        T = np.zeros((18, 18))
        for i in range(3):  # For each node
            # Apply to translational DOFs
            T[6*i:6*i+3, 6*i:6*i+3] = self._T_small
            # Apply to rotational DOFs
            T[6*i+3:6*i+6, 6*i+3:6*i+6] = self._T_small
        self._T = T
        return T
        
    def covariant_basis_vectors_triangle(self):
        """Compute and cache covariant basis vectors."""
        dN_dxi = np.array([-1, 1, 0])
        dN_deta = np.array([-1, 0, 1])
        node_coords = np.array([n.coords for n in self.nodes])  # Shape: (3, 3)

        self._a1 = np.dot(dN_dxi, node_coords)
        self._a2 = np.dot(dN_deta, node_coords)
        self._a3 = np.cross(self._a1, self._a2)

    def setup_local_coordinate_system(self, tol=1e-10):
        """Compute and cache local orthonormal basis."""
        e1 = self._a1 / (np.linalg.norm(self._a1) + tol)
        e2 = self._a2 - (self._a2.dot(e1)) * e1
        e2 = e2 / (np.linalg.norm(e2) + tol)
        e3 = np.cross(e1, e2)
        self._T_small = np.vstack((e1, e2, e3)).T  # Columns are e1, e2, e3

    #AFTER COMPUTING TRANSFORMATION MATRIX
    def compute_internal_geom(self):
      nodes = np.array([n.coords for n in self.nodes])
      translated_nodes = nodes - nodes[0]  # shift node 0 to origin
      
      #loc_nodes = np.dot(T, translated_nodes.T).T
      loc_nodes =  translated_nodes @ self._T_small

      # Compute edge vectors
      self.x21 = loc_nodes[1][0] - loc_nodes[0][0]
      self.y21 = loc_nodes[1][1] - loc_nodes[0][1]
      self.x32 = loc_nodes[2][0] - loc_nodes[1][0]
      self.y32 = loc_nodes[2][1] - loc_nodes[1][1]
      self.x13 = loc_nodes[0][0] - loc_nodes[2][0]
      self.y13 = loc_nodes[0][1] - loc_nodes[2][1]
      

      self.A2 = self.x21 * self.y32 - self.x32 * self.y21
      #print ("A2",A2)
          
      b = np.array([self.y32, self.y13, self.y21])
      c = np.array([-self.x32, -self.x13, -self.x21])
      self.dN_dx = -b / self.A2
      self.dN_dy = -c / self.A2
      
      print ("dN_dx, dN_dy", self.dN_dx, self.dN_dy)
      
    def compute_material_matrices():
        # Material matrices
        self.D_m = E / (1 - nu**2) * np.array([  # Membrane
            [1, nu, 0],
            [nu, 1, 0],
            [0, 0, (1 - nu)/2]
        ])
        
        self.D_b = E * thickness**3 / (12 * (1 - nu**2)) * np.array([  # Bending
            [1, nu, 0],
            [nu, 1, 0],
            [0, 0, (1 - nu)/2]
        ])
        
        self.D_s = E * thickness * (5/6) / (2 * (1 + nu)) * np.eye(2)  # Shear
        
      
    def compute_stiffness_matrix(self):
        """Optimized stiffness matrix computation."""
        # Ensure basis vectors are computed
        #_ = self.a1, self.a2, self.a3  # Triggers lazy computation if needed
        """Compute element stiffness matrix in global coordinates"""
        dN_dxi = np.array([-1, 1, 0])
        dN_deta = np.array([-1, 0, 1])
        print (self.nodes)
        #a1, a2, a3, norm_a3 = covariant_basis_vectors_triangle(dN_dxi, dN_deta, self.nodes)
      
        
        print ("NODES:", self.nodes)
        self.covariant_basis_vectors_triangle()
        self.setup_local_coordinate_system()
        self.create_transformation_matrix()
        self.compute_internal_geom()
        
        # Initialize stiffness matrices
        K_mem = np.zeros((18, 18))
        K_bend = np.zeros((18, 18))
        K_shear = np.zeros((18, 18))

        # Integration points (Gauss quadrature for triangles)
        #gauss_points = [(1/6, 1/6), (2/3, 1/6), (1/6, 2/3)]
        #weights = [1/6, 1/6, 1/6]
        gauss_points = [(1/3, 1/3)]  # Centroid of triangle
        weights = [0.5]              # Weight = detJ/2 (area of unit triangle)
        
        for (xi, eta), weight in zip(gauss_points, weights):
            # Reuse precomputed T_small and basis vectors
            K_mem += self.compute_membrane_stiffness(xi, eta, weight)
            K_bend += self.compute_bending_stiffness(xi, eta, weight)
        
        # NEW COUPLING

        # ~ """Compute element stiffness with unified membrane-drilling formulation"""
        # ~ # Combined B-matrix and material matrix
        # ~ B_md = self.compute_membrane_drilling_B_matrix()
        # ~ D_md = self.compute_combined_material_matrix()
        
        # ~ # Main stiffness contribution
        # ~ K_md = B_md.T @ D_md @ B_md * self.A2
        
        # ~ # Hourglass stabilization (optional)
        # ~ K_hg = self.compute_hourglass_stabilization()
        
        # ~ # Transform to global coordinates
        # ~ K_local = K_md + K_hg + ...  # Add other terms (bending, shear)

        K_shear = self.compute_shear_stiffness_mitc3()
        K_drill = self.compute_drill_stiffness()
        
        K_local = K_mem + K_bend + K_shear + K_drill
        
        self.K = self._T @ K_local @ self._T.T
        
        #print (self.K)
        # Transform to global coordinates (lazy-loaded T_small)
        return self.K

    def compute_membrane_stiffness(self, xi, eta, weight):
        K_mem = np.zeros((18, 18))
        """Optimized membrane stiffness computation."""
        # Reuse a1, a2, a3
        # Local membrane constitutive matrix
        self.D_m = self.material.E / (1 - self.material.nu**2) * np.array([
            [1, self.material.nu, 0],
            [self.material.nu, 1, 0],
            [0, 0, (1 - self.material.nu)/2]
        ])
        
        B_m = np.zeros((3, 6))  # 3 strains × 6 membrane DOFs
        #print ("a1, a2, a3",a1,a2,a3)
        for i in range(3):
            B_m[0, 2*i]   = self.dN_dx[i]  # ε_xx = ∂u/∂x
            B_m[1, 2*i+1] = self.dN_dy[i]  # ε_yy = ∂v/∂y
            B_m[2, 2*i]   = self.dN_dy[i]  # γ_xy = ∂u/∂y + ∂v/∂x
            B_m[2, 2*i+1] = self.dN_dx[i]
        
        self.B_m = B_m 
        #print ("MEMBRANE Bm", B_m)
        
        K_local = np.zeros((18, 18))
        active_dofs = [0, 1, 6, 7, 12, 13]  # u,v DOFs
        
        #K_condensed = B_m.T @ D_m @ B_m * self.detJ * weight * thickness
        K_condensed = B_m.T @ self.D_m @ B_m * self.A2 * weight * self.thickness
        for i, dof_i in enumerate(active_dofs):
            for j, dof_j in enumerate(active_dofs):
                K_local[dof_i, dof_j] = K_condensed[i, j]
        
        return K_local
        return K_mem

    def get_element_dof_indices(self):
        """Efficient DOF index lookup using list comprehension."""
        return np.concatenate([n.get_dof_indices() for n in self.nodes])
        
    def compute_bending_stiffness(self, xi, eta, weight):
        """Build complete 18×18 bending stiffness in LOCAL coordinates"""
        E, nu = self.material.E,self.material.nu
        
        # Bending constitutive matrix (local coordinates)
        self.D_b = E * self.thickness**3 / (12 * (1 - nu**2)) * np.array([
            [1,    nu,     0],
            [nu,    1,     0],
            [0,     0, (1-nu)/2]
        ])

        B_b = np.zeros((3, 18))

        for i in range(3):
            #print ("i, dNdx,dNdy",i, dN_dx[i],dN_dy[i])
            # Curvatures in local basis (κ_11, κ_22, κ_12)
            B_b[0, 6*i+4] = -self.dN_dx[i]  # κ_11 = -∂θ2/∂g1
            B_b[1, 6*i+3] = self.dN_dy[i]   # κ_22 = ∂θ1/∂g2
            B_b[2, 6*i+3] = self.dN_dx[i]   # κ_12 terms
            B_b[2, 6*i+4] = -self.dN_dy[i]
        
        # Compute stiffness matrix
        #K_b = B_b.T @ D_b @ B_b * detJ * weight
        K_b = B_b.T @ self.D_b @ B_b * self.A2 * weight
        
        self.B_b = B_b
        # Initialize B matrix (3 curvatures × 6 rotational DOFs)
        B_bc = np.zeros((3, 6))
        
        for i in range(3):
            # Compute Cartesian derivatives
            # dN_dx = invJ[0,0]*dN_dxi[i] + invJ[0,1]*dN_deta[i]
            # dN_dy = invJ[1,0]*dN_dxi[i] + invJ[1,1]*dN_deta[i]
            
            # Standard curvature definitions:
            B_bc[0, 2*i+1] =  -self.dN_dx[i]  # κ_xx = -∂θy/∂x
            B_bc[1, 2*i] = self.dN_dy[i]     # κ_yy = ∂θx/∂y
            B_bc[2, 2*i] = self.dN_dx[i]     # κ_xy terms
            B_bc[2, 2*i+1] = -self.dN_dy[i]
        
        self.B_b = B_bc
        
        # Compute stiffness matrix
        #K_bc = B_bc.T @ D_b @ B_bc * detJ * weight
        K_bc = B_bc.T @ self.D_b @ B_bc * self.A2 * weight
        #print ("B_bc", B_bc)
        #print ("K_bc", K_bc)
        return K_b
        
    def compute_shear_stiffness_mitc3(self):
        E, nu = self.material.E,self.material.nu
        G = E / (2 * (1 + nu))
        k_shear = 5/6  # Shear correction factor


        # Tying points in natural coordinates (ξ, η)
        tying_points = [(0.5, 0.0), (0.5, 0.5), (0.0, 0.5)]
        edge_nodes = [(0, 1), (1, 2), (2, 0)]

        K_shear   = np.zeros((18, 18))
        K_shear_c = np.zeros((9, 9))

        
        x31 = -self.x13
        y31 = -self.y13
        # Compute the angles
        phi1 = np.arctan2(self.y21, self.x21)
        phi2 = 0.5 * np.pi - np.arctan2(x31, y31)

        # Create and assign BsO matrix
        BsO = np.zeros((2, 2))
        BsO[0, 0] = np.sin(phi2)
        BsO[0, 1] = -np.sin(phi1)
        BsO[1, 0] = -np.cos(phi2)
        BsO[1, 1] = np.cos(phi1)

        # Create and assign BsC matrix
        BsC = np.zeros((2, 2))
        
        BsC[0, 0] = np.sqrt(self.x13**2 + self.y13**2) / self.A2
        BsC[1, 1] = np.sqrt(self.x21**2 + self.y21**2) / self.A2
        ################################################### SCALE 
        
        self.D_s = (k_shear * G * self.thickness) * np.eye(2)  # Shear

            # ~ b = np.array([self.y32, self.y13, self.y21])
            # ~ c = np.array([-self.x32, -self.x13, -self.x21])
            # ~ self.dN_dx = b / self.A2
            # ~ self.dN_dy = c / self.A2
  
        #print ("Shear Interpolation Type OPS")
        for xi, eta in tying_points:
            cross = True
            B_s = np.zeros((2, 9))
            B_s_OCT = np.zeros((2, 9))

            ###Row 0: shear strain γ_xz [-1 1 0] SO DOF 6 is not written
            ### DIFFERENCE: CROSS DERIVATIVES ARE NOT CONSTANT 
            ### CROSS TERMS COUPLE 0,1 0,4 and 0,
            B_s[0, 0] = -1.0  # θ_x at node 1          
            if (cross): B_s[0, 1] = -(self.y21 + self.y32 * xi) / 2.0
            B_s[0, 2] = (self.x21 + self.x32 * xi) / 2.0        
            B_s[0, 3] = 1.0  # θ_x at node 2
            if (cross): B_s[0, 4] = -(self.y21 + self.y13 * xi) / 2.0
            B_s[0, 5] = (self.x21 + self.x13 * xi) / 2.0
            if (cross): B_s[0, 7] = -(self.y21 * xi) / 2.0
            B_s[0, 8] = (self.x21 * xi) / 2.0
            # #############################
            # # Row 1: shear strain γ_yz ##[-1, 0, 1] SO DOF 3 not filled
            B_s[1, 0] = -1.0  # θ_y at node 1
            B_s[1, 1] = (self.y13 + self.y32 * eta) / 2.0
            if (cross): B_s[1, 2] = -(self.x13 + self.x32 * eta) / 2.0
            B_s[1, 4] = (self.y13 * eta) / 2.0
            if (cross): B_s[1, 5] = -(self.x13 * eta) / 2.0            
            B_s[1, 6] = 1.0  # θ_y at node 3
            B_s[1, 7] = (self.y13 + self.y21 * eta) / 2.0
            if (cross): B_s[1, 8] = -(self.x13 + self.x21 * eta) / 2.0
            #print ("BS",B_s)
            B_s_OCT = BsO @ BsC @ B_s
              
            K_shear_c += (B_s_OCT.T @ B_s_OCT) * (k_shear * G * self.thickness) * (self.A2 / 6)        

        self.B_s = B_s_OCT
        active_dofs = [2,3,4,  8,9,10, 14,15,16]
        for i, dof_i in enumerate(active_dofs):
            for j, dof_j in enumerate(active_dofs):
                K_shear[dof_i, dof_j] = K_shear_c[i, j]
            
        #print("K_shear",K_shear)
        
        return K_shear
     
    # ~ def compute_drill_stiffness(self):
    
        # ~ k_drill_local = np.zeros((18, 18))
        
            # ~ # Add drilling stiffness (Gt/2)
        # ~ E, nu = self.material.E,self.material.nu
        # ~ G = E / (2 * (1 + nu))
        # ~ drill_stiff = 0.5 * self.thickness * G
        
        # ~ for i in range(3):
            # ~ k_drill_local[6*i+5, 6*i+5] += drill_stiff * self.A2
        # ~ return k_drill_local


    def compute_pure_drilling_components(self):
        """Compute only the drilling-related stiffness components"""
        # Shape functions at centroid (constant for linear triangle)
        N = np.array([1/3, 1/3, 1/3])
        
        # Initialize drilling B-matrix (1×18)
        B_drill = np.zeros((1, 18))  # Only γ_drill strain
        
        # Drilling strain-displacement relationship
        for i in range(3):
            B_drill[0, 6*i] = -0.5 * self.dN_dy[i]    # ½(∂v/∂x) component
            B_drill[0, 6*i+1] = 0.5 * self.dN_dx[i]   # -½(∂u/∂y) component
            B_drill[0, 6*i+5] = -N[i]                 # -θz component
        
        # Drilling material term (scalar)
        alpha = 0.1  # Stabilization factor
        D_drill = alpha * self.material.E * self.thickness**3 / 12
        
        # Pure drilling stiffness matrix
        K_drill = B_drill.T @ B_drill * D_drill * self.A2
        
        return K_drill

    def compute_element_distortion(self):
        """Compute metrics for adaptive hourglass control."""
        # Edge lengths
        L1 = np.linalg.norm(self.nodes[1].coords - self.nodes[0].coords)
        L2 = np.linalg.norm(self.nodes[2].coords - self.nodes[1].coords)
        L3 = np.linalg.norm(self.nodes[0].coords - self.nodes[2].coords)
        
        # Aspect ratio (max edge / min edge)
        max_edge = max(L1, L2, L3)
        min_edge = min(L1, L2, L3)
        aspect_ratio = max_edge / min_edge
        
        # Skewness (deviation from ideal equilateral triangle)
        ideal_area = (np.sqrt(3)/4) * (min_edge**2)
        skewness = abs(self.A2 - ideal_area) / ideal_area
        
        return aspect_ratio, skewness

    def compute_hourglass_stabilization(self):
        """ANSYS-style adaptive hourglass control."""
        aspect_ratio, skewness = self.compute_element_distortion()
        G = self.material.G
        
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
            Bhx[6*i+5] = self.dN_dx[i]  # ∂θz/∂x
            Bhy[6*i+5] = self.dN_dy[i]  # ∂θz/∂y
        
        K_hg = beta * G * self.thickness * self.A2 * (
            np.outer(Bhx, Bhx) + np.outer(Bhy, Bhy))
        
        return K_hg

    # CONSTANT
    # OPENSEES
    def compute_constant_hourglass_stabilization(self):
        """Compute hourglass control matrix for drilling DOFs"""
        # Initialize hourglass B-matrices
        Bhx = np.zeros(18)  # ∂θz/∂x terms
        Bhy = np.zeros(18)  # ∂θz/∂y terms
        
        for i in range(3):
            Bhx[6*i+5] = self.dN_dx[i]  # ∂θz/∂x
            Bhy[6*i+5] = self.dN_dy[i]  # ∂θz/∂y
        
        # Stabilization parameters
        alpha = 0.05
        G = self.material.E / (2 * (1 + self.material.nu))
        
        # Hourglass stiffness
        K_hg = alpha * G * self.thickness * self.A2 * (
            np.outer(Bhx, Bhx) + np.outer(Bhy, Bhy))
        
        return K_hg
    
    def compute_drill_stiffness(self):
      """Combined physical drilling and hourglass stabilization"""
      K_drill = self.compute_pure_drilling_components()
      K_hg = self.compute_hourglass_stabilization()
      return K_drill + K_hg
      
    def get_element_dof_indices(self):
        """Return global DOF indices for all nodes in element"""
        return np.concatenate([n.get_dof_indices() for n in self.nodes])
        
    def compute_stresses(self, global_displacements):
        """Compute element stresses at integration points"""
        # Transform global displacements to local coordinates
        elem_dofs = self.get_element_dof_indices()
        print ("Element DOFS: ",elem_dofs)
        U_global = global_displacements[elem_dofs]
        U_local = U_global @ self._T
        print ("U local", U_local)
        # Compute stresses using B matrices and material properties
        # Implementation similar to your calculate_element_stresses function...
        # # Initialize stress storage

        
        # # Gauss points
        gauss_points = [(1/3, 1/3)]  # Centroid of triangle
        weights = [0.5]              # Weight = detJ/2 (area of unit triangle)
    
    
        for (xi, eta), weight in zip(gauss_points, weights):
              
            active_dofs_m = [0, 1, 6, 7, 12, 13]  # u,v DOFs for membrane part
            U_membrane = U_local[active_dofs_m]  # Condensed displacement vector (6x1)
            
            membrane_strains = self.B_m @ U_membrane
            membrane_stresses = self.D_m @ membrane_strains

            
            active_dofs_b = [3, 4, 9, 10, 15, 16]  # u,v DOFs for membrane part
            U_bending = U_local[active_dofs_b]  # Condensed displacement vector (6x1)
                
            bending_curvatures = self.B_b @ U_bending  # [κ_xx, κ_yy, 2κ_xy]
            bending_moments = self.D_b @ bending_curvatures  # [M_xx, M_yy, M_xy] bending moments per unit length

            # Convert bending moments to stresses at faces
            # Bending stress formula: σ_b = (M * z) / I, where I = t³/12
            # For z = ±t/2, this becomes σ_b = ± (6M)/t²
            z = self.thickness/2
            bending_stresses = (6.0/self.thickness**2) * bending_moments[:3]  # [σ_xx_b, σ_yy_b, σ_xy_b]
 
             # Total stresses at faces (membrane + bending)
            upper_face = membrane_stresses + bending_stresses
            lower_face = membrane_stresses - bending_stresses
            
            active_dofs_s = [2, 3, 4, 8, 9, 10, 14, 15, 16]  # u,v DOFs for membrane part
            U_shear = U_local[active_dofs_s]  # Condensed displacement vector (6x1)
            shear_strains = self.B_s @ U_shear  
            #shear_strains_local = J_pseudo @ shear_strains_cov  [γ_xz, γ_yz]
            
            ##Shear stresses
            shear_forces = self.D_s @ shear_strains  # [Q_x, Q_y] hear forces per unit length, not shear stresses yet.
            shear_stresses = shear_forces / self.thickness

            # Combine into full stress vectors for each face
            # Format: [σ_xx, σ_yy, σ_xy, σ_xz, σ_yz, σ_zz]
            upper_stress = np.array([
                upper_face[0],      # σ_xx
                upper_face[1],      # σ_yy
                upper_face[2],      # σ_xy
                shear_stresses[0],  # σ_xz
                shear_stresses[1],  # σ_yz
                0.0                # σ_zz (assumed zero)
            ])
            
            lower_stress = np.array([
                lower_face[0],      # σ_xx
                lower_face[1],      # σ_yy
                lower_face[2],      # σ_xy
                shear_stresses[0],  # σ_xz
                shear_stresses[1],  # σ_yz
                0.0                # σ_zz (assumed zero)
            ])

            # Calculate upper and lower face stresses
            z = self.thickness/2  # Distance to outer surfaces

            # ~ # Compute principal stresses
            # ~ for face_stress in [upper_stress, lower_stress]:
                # ~ # Get in-plane stress components
                # ~ stress_matrix = np.array([
                    # ~ [face_stress[0], face_stress[2], 0],
                    # ~ [face_stress[2], face_stress[1], 0],
                    # ~ [0, 0, face_stress[5]]
                # ~ ])
                
                # ~ # Compute principal stresses
                # ~ principal_stresses = np.linalg.eigvalsh(stress_matrix)
                # ~ self.stresses['max'].append(np.max(principal_stresses))
                # ~ self.stresses['min'].append(np.min(principal_stresses))
                
            
            self.stresses['membrane'].append(membrane_stresses)
            self.stresses['bending'].append(bending_stresses)
            self.stresses['shear'].append(shear_stresses)
            self.stresses['upper'].append(upper_stress)
            self.stresses['lower'].append(lower_stress)
        return self.stresses

    def get_local_stress_tensor(sigma_xx, sigma_yy, sigma_xy, sigma_zz=0, sigma_yz=0, sigma_xz=0):
        """
        Compose a 3x3 stress tensor in the local coordinate system.
        Off-diagonal terms (xz, yz, zz) can be zero for shell in-plane stress.
        """
        sigma_local = np.array([
            [sigma_xx, sigma_xy, sigma_xz],
            [sigma_xy, sigma_yy, sigma_yz],
            [sigma_xz, sigma_yz, sigma_zz]
        ])
        return sigma_local

    def rotate_stress_tensor(sigma_local, T):
        """
        Rotate a 3x3 stress tensor from local to global using T.
        sigma_global = T * sigma_local * T.T
        """
        return T @ sigma_local @ T.T

class AnalysisCase:
    def __init__(self):
        self.loads = []
        self.bcs = []
        self.load_steps = 1      

# Usage Example
if __name__ == "__main__":
    # Create model
    model = Model()
    
    # Add material
    model.add_material("steel", E=2.0e11, nu=0.3)
    
    # Add nodes
    n0 = model.add_node(0, 0, 0)
    n1 = model.add_node(1, 0, 0)
    n2 = model.add_node(0, 1, 0)
    
    # Add element
    model.add_element([n0, n1, n2], "steel", thickness=0.1)
    
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
    
    print("Displacements:", U)

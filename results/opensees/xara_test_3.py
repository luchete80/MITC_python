from openseespy import opensees as ops
import xara

ops.model('basic', '-ndm', 3, '-ndf', 6)

# Define material
E = 200.0e9
nu = 0.3
h = 0.1
ops.nDMaterial('ElasticIsotropic', 1, E, nu)

# Nodes (4 nodes for a quad shell)
ops.node(1, 0.0, 0.0, 0.0)
ops.node(2, 1.0, 0.0, 0.0)
ops.node(3, 1.0, 1.0, 0.0)
ops.node(4, 0.0, 1.0, 0.0)

# ShellMITC4 element
ops.element('ShellMITC4', 1, 1, 2, 3, 4, 1)

# Fix nodes
ops.fix(1, 1,1,1,1,1,1)
ops.fix(4, 1,1,1,1,1,1)

# Load
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)
ops.load(2, 0, 0, -50.0, 0, 0, 0)

# Analysis
ops.system('UmfPack')
ops.analyze(1)

# Now try stiffness extraction
print("Available responses:", ops.eleResponse(1, 'list'))  # Should work
Ke = ops.eleResponse(1, 'stiffness')  # 24x24 matrix (4 nodes × 6 DOFs)
print("Stiffness matrix:", Ke)

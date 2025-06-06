#import xara 
from openseespy import opensees as ops
import numpy as np
#ops = xara.Model()
#model = ops.Model(ndm=2, ndf=3)

#ops.model(ndm=3, print(xara.list_available_elements())  # If XARA provides thisndf=3)

ops.model('basic', '-ndm', 3, '-ndf', 6)
E = 200.0e9
h = 0.1
ops.section('ElasticMembranePlateSection', 1, E, 0.0, h)

angle = 10.0

a = angle *3.1415926/180.0

ops.node(1, 0.0, 0.0,0)
ops.node(2, np.cos(a), 0.0,np.sin(a))
ops.node(3, 0.0, 1.0,0)
conn = (1,2,3)
ops.element('ASDShellT3',1,  1,2,3,  1, '-linear')
#ops.element('ASDShellT3',1,  1,2,3,  1, '--linear')



ops.fix(1, 1,1,1,1,1,1)
ops.fix(3, 1,1,1,1,1,1)

#ops.eval("print -json")

ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)

ops.load(2, 200000,200000,200000, 0,0,0)




# analysis
duration = 1.0
nsteps = 1
dt = duration/nsteps
dt_record = 0.2
ops.constraints('Transformation')
ops.numberer('RCM')
ops.system('UmfPack')
ops.test('NormDispIncr', 1.0e-5, 100, 0)
ops.algorithm('Newton')
ops.integrator('LoadControl', dt)
ops.analysis('Static')

# Compute the local stiffness matrix for the element


print("Available responses:", ops.eleResponse(1, 'list'))

#ops.printModel('ele', 1, '-flag', 2)


ops.analyze(1)

Ke = ops.eleResponse(1, 'tangentStiffness')
print("Ke",Ke)
print("Available responses:", ops.eleResponse(1, 'list'))

# Get the element stiffness matrix
Ke = ops.eleResponse(1, 'stiffness')
print("Element stiffness matrix:", Ke)

print ("DispX", ops.nodeDisp(2, 1))
print ("DispX", ops.nodeDisp(2, 2))
print ("DispX", ops.nodeDisp(2, 3))
print ("RotX", ops.nodeDisp(2, 4))
print ("RotY", ops.nodeDisp(2, 5))
print ("RotZ", ops.nodeDisp(2, 6))
stress = ops.eleResponse(1, 'fiberStress')
print("Stress result:", stress)
print("Available responses:", ops.eleResponse(1, 'list'))

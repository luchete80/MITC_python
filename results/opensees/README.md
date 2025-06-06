Currently, the OpenSeesPy Python wrapper (installed via pip) does not natively support extracting the stiffness matrix of shell elements (including the 3-node ASDShellT3). However, you can implement a workaround or extend the wrapper yourself. Here’s how:

Option 1: Use eleResponse (If Available)
Some shell elements allow extracting stiffness via eleResponse. Try:

python
import openseespy.opensees as ops

ops.eleResponse(eleTag, "stiffness")  # Replace eleTag with your shell element tag
Limitations:

Not all shell elements support this (e.g., ASDShellT3 may not expose it).

Output format varies by element type.

Option 2: Modify OpenSeesPy to Expose Stiffness
To add stiffness extraction for ASDShellT3 (or other shells):

Step 1: Clone and Edit the Source
Clone the OpenSeesPy repository:

bash
git clone https://github.com/zhuminjie/OpenSeesPy.git
cd OpenSeesPy
Edit ASDShellT3.cpp (source file):
Add a method to return the stiffness matrix (e.g., getStiffness).

Expose the Method in Python:
Edit the OpenSeesPy wrapper (openseespy/opensees.py) to call this new C++ method.

Step 2: Recompile OpenSeesPy
bash
python setup.py install
Option 3: Global Stiffness Matrix Workaround
If modifying C++ isn’t feasible, assemble the global stiffness matrix column-by-column:

python
import numpy as np

def getGlobalStiffness():
    ops.wipeAnalysis()
    ops.system("FullGeneral")  # Required for accurate stiffness
    ops.analysis("Static")
    
    nDOF = ops.systemSize()  # Total DOFs in the model
    K = np.zeros((nDOF, nDOF))
    
    for j in range(nDOF):
        u = np.zeros(nDOF)
        u[j] = 1.0  # Perturb one DOF
        ops.loadConst("-time", 0.0)
        ops.integrator("DisplacementControl", 1, 1, u[j])
        ops.analyze(1)
        K[:, j] = ops.reactions()  # Stiffness column = reactions
    
    return K
Limitations:

Computationally expensive for large models.

Only works for static analyses.

Option 4: Use OpenSees Tcl + Python Subprocess
If Python access is critical, call OpenSees via Tcl and parse the output:

python
import subprocess

tcl_script = """
model BasicBuilder -ndm 3 -ndf 6
node 1 0 0 0
node 2 1 0 0
node 3 0 1 0
element ASDShellT3 1 1 2 3 1.0 0.0
print ele 1
"""

result = subprocess.run(["opensees"], input=tcl_script, text=True, capture_output=True)
print(result.stdout)
Key Notes
For ASDShellT3: Check if print ele <tag> in Tcl outputs stiffness. If not, the element may need C++ modifications.

OpenSeesPy Limitations: The Python wrapper is a thin layer over Tcl; advanced features like element-level stiffness often require custom C++ edits.

Alternative: Consider using SOFiSTiK or Abaqus if shell stiffness access is critical.

Final Recommendation
Quick Solution: Try eleResponse or global matrix assembly.

Long-Term Solution: Modify ASDShellT3.cpp and recompile OpenSeesPy (requires C++/Python hybrid coding).

Let me know if you’d like help with specific code edits!

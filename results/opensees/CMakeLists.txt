cmake_minimum_required(VERSION 3.10)
project(ShellExample)

set(CMAKE_CXX_STANDARD 11)

# Path to your OpenSees installation
set(OPENSEES_SRC_DIR /home/weldform-pc/Numerico/OpenSees/SRC)
include_directories(
    ${OPENSEES_SRC_DIR}
    ${OPENSEES_SRC_DIR}/actor
    ${OPENSEES_SRC_DIR}/actor/actor
    ${OPENSEES_SRC_DIR}/actor/address
    ${OPENSEES_SRC_DIR}/actor/channel
    ${OPENSEES_SRC_DIR}/actor/machineBroker
    ${OPENSEES_SRC_DIR}/actor/message
    ${OPENSEES_SRC_DIR}/actor/objectBroker
    ${OPENSEES_SRC_DIR}/actor/shadow
    ${OPENSEES_SRC_DIR}/analysis
    ${OPENSEES_SRC_DIR}/analysis/algorithm
    ${OPENSEES_SRC_DIR}/analysis/algorithm/domainDecompAlgo
    ${OPENSEES_SRC_DIR}/analysis/algorithm/eigenAlgo
    ${OPENSEES_SRC_DIR}/analysis/algorithm/equiSolnAlgo
    ${OPENSEES_SRC_DIR}/analysis/algorithm/equiSolnAlgo/accelerator
    ${OPENSEES_SRC_DIR}/analysis/analysis
    ${OPENSEES_SRC_DIR}/analysis/dof_grp
    ${OPENSEES_SRC_DIR}/analysis/fe_ele
    ${OPENSEES_SRC_DIR}/analysis/fe_ele/lagrange
    ${OPENSEES_SRC_DIR}/analysis/fe_ele/penalty
    ${OPENSEES_SRC_DIR}/analysis/fe_ele/transformation
    ${OPENSEES_SRC_DIR}/analysis/handler
    ${OPENSEES_SRC_DIR}/analysis/integrator
    ${OPENSEES_SRC_DIR}/analysis/model
    ${OPENSEES_SRC_DIR}/analysis/model/simple
    ${OPENSEES_SRC_DIR}/analysis/numberer
    ${OPENSEES_SRC_DIR}/api
    ${OPENSEES_SRC_DIR}/convergenceTest
    ${OPENSEES_SRC_DIR}/coordTransformation
    ${OPENSEES_SRC_DIR}/damage
    ${OPENSEES_SRC_DIR}/damping
    ${OPENSEES_SRC_DIR}/database
    ${OPENSEES_SRC_DIR}/doc
    ${OPENSEES_SRC_DIR}/doc/fig_files
    ${OPENSEES_SRC_DIR}/domain
    ${OPENSEES_SRC_DIR}/domain/component
    ${OPENSEES_SRC_DIR}/domain/constraints
    ${OPENSEES_SRC_DIR}/domain/domain
    ${OPENSEES_SRC_DIR}/domain/domain/partitioned
    ${OPENSEES_SRC_DIR}/domain/domain/single
    ${OPENSEES_SRC_DIR}/domain/groundMotion
    ${OPENSEES_SRC_DIR}/domain/IGA
    ${OPENSEES_SRC_DIR}/domain/load
    ${OPENSEES_SRC_DIR}/domain/loadBalancer
    ${OPENSEES_SRC_DIR}/domain/node
    ${OPENSEES_SRC_DIR}/domain/partitioner
    ${OPENSEES_SRC_DIR}/domain/pattern
    ${OPENSEES_SRC_DIR}/domain/pattern/drm
    ${OPENSEES_SRC_DIR}/domain/region
    ${OPENSEES_SRC_DIR}/domain/subdomain
    ${OPENSEES_SRC_DIR}/element
    ${OPENSEES_SRC_DIR}/element/absorbentBoundaries
    ${OPENSEES_SRC_DIR}/element/adapter
    ${OPENSEES_SRC_DIR}/element/beam2d
    ${OPENSEES_SRC_DIR}/element/beam3d
    ${OPENSEES_SRC_DIR}/element/beamWithHinges
    ${OPENSEES_SRC_DIR}/element/brick
    ${OPENSEES_SRC_DIR}/element/catenaryCable
    ${OPENSEES_SRC_DIR}/element/CEqElement
    ${OPENSEES_SRC_DIR}/element/componentElement
    ${OPENSEES_SRC_DIR}/element/dispBeamColumn
    ${OPENSEES_SRC_DIR}/element/dispBeamColumnInt
    ${OPENSEES_SRC_DIR}/element/dmglib
    ${OPENSEES_SRC_DIR}/element/elasticBeamColumn
    ${OPENSEES_SRC_DIR}/element/elastomericBearing
    ${OPENSEES_SRC_DIR}/element/feap
    ${OPENSEES_SRC_DIR}/element/fmkPlanarTruss
    ${OPENSEES_SRC_DIR}/element/forceBeamColumn
    ${OPENSEES_SRC_DIR}/element/fourNodeQuad
    ${OPENSEES_SRC_DIR}/element/frictionBearing
    ${OPENSEES_SRC_DIR}/element/frictionBearing/frictionModel
    ${OPENSEES_SRC_DIR}/element/generic
    ${OPENSEES_SRC_DIR}/element/gradientInelasticBeamColumn
    ${OPENSEES_SRC_DIR}/element/HUelements
    ${OPENSEES_SRC_DIR}/element/IGA
    ${OPENSEES_SRC_DIR}/element/joint
    ${OPENSEES_SRC_DIR}/element/LHMYS
    ${OPENSEES_SRC_DIR}/element/masonry
    ${OPENSEES_SRC_DIR}/element/mefi
    ${OPENSEES_SRC_DIR}/element/mixedBeamColumn
    ${OPENSEES_SRC_DIR}/element/mvlem
    ${OPENSEES_SRC_DIR}/element/nonlinearBeamColumn
    ${OPENSEES_SRC_DIR}/element/nonlinearBeamColumn/matrixutil
    ${OPENSEES_SRC_DIR}/element/nonlinearBeamColumn/tcl
    ${OPENSEES_SRC_DIR}/element/PFEMElement
    ${OPENSEES_SRC_DIR}/element/pipe
    ${OPENSEES_SRC_DIR}/element/PML
    ${OPENSEES_SRC_DIR}/element/pyMacro
    ${OPENSEES_SRC_DIR}/element/RockingBC
    ${OPENSEES_SRC_DIR}/element/shell
    ${OPENSEES_SRC_DIR}/element/surfaceLoad
    ${OPENSEES_SRC_DIR}/element/tetrahedron
    ${OPENSEES_SRC_DIR}/element/triangle
    ${OPENSEES_SRC_DIR}/element/truss
    ${OPENSEES_SRC_DIR}/element/twoNodeLink
    ${OPENSEES_SRC_DIR}/element/updatedLagrangianBeamColumn
    ${OPENSEES_SRC_DIR}/element/UP-ucsd
    ${OPENSEES_SRC_DIR}/element/UWelements
    ${OPENSEES_SRC_DIR}/element/XMUelements
    ${OPENSEES_SRC_DIR}/element/zeroLength
    ${OPENSEES_SRC_DIR}/graph
    ${OPENSEES_SRC_DIR}/graph/graph
    ${OPENSEES_SRC_DIR}/graph/numberer
    ${OPENSEES_SRC_DIR}/graph/partitioner
    ${OPENSEES_SRC_DIR}/handler
    ${OPENSEES_SRC_DIR}/interpreter
    ${OPENSEES_SRC_DIR}/java
    ${OPENSEES_SRC_DIR}/java/kepler
    ${OPENSEES_SRC_DIR}/java/kepler/actors
    ${OPENSEES_SRC_DIR}/java/kepler/actors/OpenSeesEvaluator
    ${OPENSEES_SRC_DIR}/java/kepler/actors/OpenSeesParameter
    ${OPENSEES_SRC_DIR}/java/kepler/actors/StringToDouble
    ${OPENSEES_SRC_DIR}/java/kepler/opensees
    ${OPENSEES_SRC_DIR}/machine
    ${OPENSEES_SRC_DIR}/material
    ${OPENSEES_SRC_DIR}/material/nD
    ${OPENSEES_SRC_DIR}/material/nD/ASDPlasticMaterial3D
    ${OPENSEES_SRC_DIR}/material/nD/ASDPlasticMaterial3D/ElasticityModels
    ${OPENSEES_SRC_DIR}/material/nD/ASDPlasticMaterial3D/PlasticFlowDirections
    ${OPENSEES_SRC_DIR}/material/nD/ASDPlasticMaterial3D/YieldFunctions
    ${OPENSEES_SRC_DIR}/material/nD/cyclicSoil
    ${OPENSEES_SRC_DIR}/material/nD/feap
    ${OPENSEES_SRC_DIR}/material/nD/matCMM
    ${OPENSEES_SRC_DIR}/material/nD/OrthotropicRotatingAngleConcreteT2DMaterial01
    ${OPENSEES_SRC_DIR}/material/nD/reinforcedConcretePlaneStress
    ${OPENSEES_SRC_DIR}/material/nD/SmearedSteelDoubleLayerT2DMaterial01
    ${OPENSEES_SRC_DIR}/material/nD/soil
    ${OPENSEES_SRC_DIR}/material/nD/stressDensityModel
    ${OPENSEES_SRC_DIR}/material/nD/UANDESmaterials
    ${OPENSEES_SRC_DIR}/material/nD/UWmaterials
    ${OPENSEES_SRC_DIR}/material/section
    ${OPENSEES_SRC_DIR}/material/section/fiber
    ${OPENSEES_SRC_DIR}/material/section/integration
    ${OPENSEES_SRC_DIR}/material/section/LayeredMembraneSection
    ${OPENSEES_SRC_DIR}/material/section/repres
    ${OPENSEES_SRC_DIR}/material/section/repres/cell
    ${OPENSEES_SRC_DIR}/material/section/repres/patch
    ${OPENSEES_SRC_DIR}/material/section/repres/reinfBar
    ${OPENSEES_SRC_DIR}/material/section/repres/reinfLayer
    ${OPENSEES_SRC_DIR}/material/section/repres/section
    ${OPENSEES_SRC_DIR}/material/section/yieldSurface
    ${OPENSEES_SRC_DIR}/material/uniaxial
    ${OPENSEES_SRC_DIR}/material/uniaxial/backbone
    ${OPENSEES_SRC_DIR}/material/uniaxial/drain
    ${OPENSEES_SRC_DIR}/material/uniaxial/fedeas
    ${OPENSEES_SRC_DIR}/material/uniaxial/limitState
    ${OPENSEES_SRC_DIR}/material/uniaxial/limitState/limitCurve
    ${OPENSEES_SRC_DIR}/material/uniaxial/PY
    ${OPENSEES_SRC_DIR}/material/uniaxial/snap
    ${OPENSEES_SRC_DIR}/material/uniaxial/stiffness
    ${OPENSEES_SRC_DIR}/material/uniaxial/strength
    ${OPENSEES_SRC_DIR}/material/uniaxial/unloading
    ${OPENSEES_SRC_DIR}/material/yieldSurface
    ${OPENSEES_SRC_DIR}/material/yieldSurface/evolution
    ${OPENSEES_SRC_DIR}/material/yieldSurface/plasticHardeningMaterial
    ${OPENSEES_SRC_DIR}/material/yieldSurface/yieldSurfaceBC
    ${OPENSEES_SRC_DIR}/matrix
    ${OPENSEES_SRC_DIR}/modelbuilder
    ${OPENSEES_SRC_DIR}/modelbuilder/tcl
    ${OPENSEES_SRC_DIR}/optimization
    ${OPENSEES_SRC_DIR}/optimization/domain
    ${OPENSEES_SRC_DIR}/optimization/domain/component
    ${OPENSEES_SRC_DIR}/optimization/SNOPT
    ${OPENSEES_SRC_DIR}/optimization/tcl
    ${OPENSEES_SRC_DIR}/recorder
    ${OPENSEES_SRC_DIR}/recorder/response
    ${OPENSEES_SRC_DIR}/reliability
    ${OPENSEES_SRC_DIR}/reliability/analysis
    ${OPENSEES_SRC_DIR}/reliability/analysis/analysis
    ${OPENSEES_SRC_DIR}/reliability/analysis/analysis/system
    ${OPENSEES_SRC_DIR}/reliability/analysis/convergenceCheck
    ${OPENSEES_SRC_DIR}/reliability/analysis/curvature
    ${OPENSEES_SRC_DIR}/reliability/analysis/designPoint
    ${OPENSEES_SRC_DIR}/reliability/analysis/direction
    ${OPENSEES_SRC_DIR}/reliability/analysis/gFunction
    ${OPENSEES_SRC_DIR}/reliability/analysis/gradient
    ${OPENSEES_SRC_DIR}/reliability/analysis/hessian
    ${OPENSEES_SRC_DIR}/reliability/analysis/hessianApproximation
    ${OPENSEES_SRC_DIR}/reliability/analysis/meritFunction
    ${OPENSEES_SRC_DIR}/reliability/analysis/misc
    ${OPENSEES_SRC_DIR}/reliability/analysis/randomNumber
    ${OPENSEES_SRC_DIR}/reliability/analysis/rootFinding
    ${OPENSEES_SRC_DIR}/reliability/analysis/sensitivity
    ${OPENSEES_SRC_DIR}/reliability/analysis/stepSize
    ${OPENSEES_SRC_DIR}/reliability/analysis/transformation
    ${OPENSEES_SRC_DIR}/reliability/domain
    ${OPENSEES_SRC_DIR}/reliability/domain/components
    ${OPENSEES_SRC_DIR}/reliability/domain/distributions
    ${OPENSEES_SRC_DIR}/reliability/domain/filter
    ${OPENSEES_SRC_DIR}/reliability/domain/functionEvaluator
    ${OPENSEES_SRC_DIR}/reliability/domain/modulatingFunction
    ${OPENSEES_SRC_DIR}/reliability/domain/performanceFunction
    ${OPENSEES_SRC_DIR}/reliability/domain/spectrum
    ${OPENSEES_SRC_DIR}/reliability/domain/storage
    ${OPENSEES_SRC_DIR}/reliability/FEsensitivity
    ${OPENSEES_SRC_DIR}/reliability/tcl
    ${OPENSEES_SRC_DIR}/remote
    ${OPENSEES_SRC_DIR}/renderer
    ${OPENSEES_SRC_DIR}/runtime
    ${OPENSEES_SRC_DIR}/runtime/commands
    ${OPENSEES_SRC_DIR}/runtime/commands/analysis
    ${OPENSEES_SRC_DIR}/runtime/commands/domain
    ${OPENSEES_SRC_DIR}/runtime/commands/domain/database
    ${OPENSEES_SRC_DIR}/runtime/commands/domain/loading
    ${OPENSEES_SRC_DIR}/runtime/commands/domain/loading/drm
    ${OPENSEES_SRC_DIR}/runtime/commands/modeling
    ${OPENSEES_SRC_DIR}/runtime/commands/modeling/invoking
    ${OPENSEES_SRC_DIR}/runtime/commands/parallel
    ${OPENSEES_SRC_DIR}/runtime/commands/utilities
    ${OPENSEES_SRC_DIR}/runtime/parsing
    ${OPENSEES_SRC_DIR}/runtime/python
    ${OPENSEES_SRC_DIR}/runtime/runtime
    ${OPENSEES_SRC_DIR}/runtime/runtime/SectionBuilder
    ${OPENSEES_SRC_DIR}/runtime/runtime/SectionBuilder/cell
    ${OPENSEES_SRC_DIR}/runtime/runtime/SectionBuilder/patch
    ${OPENSEES_SRC_DIR}/runtime/runtime/SectionBuilder/reinfBar
    ${OPENSEES_SRC_DIR}/runtime/runtime/SectionBuilder/reinfLayer
    ${OPENSEES_SRC_DIR}/scripts
    ${OPENSEES_SRC_DIR}/string
    ${OPENSEES_SRC_DIR}/system_of_eqn
    ${OPENSEES_SRC_DIR}/system_of_eqn/eigenSOE
    ${OPENSEES_SRC_DIR}/system_of_eqn/linearSOE
    ${OPENSEES_SRC_DIR}/system_of_eqn/linearSOE/bandGEN
    ${OPENSEES_SRC_DIR}/system_of_eqn/linearSOE/bandSPD
    ${OPENSEES_SRC_DIR}/system_of_eqn/linearSOE/cg
    ${OPENSEES_SRC_DIR}/system_of_eqn/linearSOE/diagonal
    ${OPENSEES_SRC_DIR}/system_of_eqn/linearSOE/fullGEN
    ${OPENSEES_SRC_DIR}/system_of_eqn/linearSOE/itpack
    ${OPENSEES_SRC_DIR}/system_of_eqn/linearSOE/mumps
    ${OPENSEES_SRC_DIR}/system_of_eqn/linearSOE/pardiso
    ${OPENSEES_SRC_DIR}/system_of_eqn/linearSOE/petsc
    ${OPENSEES_SRC_DIR}/system_of_eqn/linearSOE/profileSPD
    ${OPENSEES_SRC_DIR}/system_of_eqn/linearSOE/sparseGEN
    ${OPENSEES_SRC_DIR}/system_of_eqn/linearSOE/sparseSYM
    ${OPENSEES_SRC_DIR}/system_of_eqn/linearSOE/umfGEN
    ${OPENSEES_SRC_DIR}/tagged
    ${OPENSEES_SRC_DIR}/tagged/storage
    ${OPENSEES_SRC_DIR}/tcl
    ${OPENSEES_SRC_DIR}/tcl/include
    ${OPENSEES_SRC_DIR}/unittest
    ${OPENSEES_SRC_DIR}/utility
)


set(OPENSEES_BIN_DIR /home/weldform-pc/Numerico/OpenSees_bin)
link_directories(${OPENSEES_BIN_DIR})

add_executable(shell_example shell_example.cpp)
target_link_libraries(shell_example OpenSeesLIB OPS_InterpTcl G3 CSPARSE)

#include <OPS_Stream.h>
#include <StandardStream.h>
#include <Domain.h>
#include <ElasticMembranePlateSection.h>
#include <SectionForceDeformation.h>

#include <ASDShellT3.h>
#include <SP_Constraint.h>
#include <LoadPattern.h>
#include <TimeSeries.h>
#include <LinearSeries.h>
#include <NodalLoad.h>
#include <AnalysisModel.h>
#include <CTestNormDispIncr.h>
#include <NewtonRaphson.h>
#include <StaticAnalysis.h>
#include <TransformationConstraintHandler.h>
#include <RCM.h>
#include <DOF_Numberer.h>
#include <BandGenLinSolver.h>
#include <BandGenLinSOE.h>
#include <BandGenLinLapackSolver.h>  // Alternative 2: LAPACK-based solve
#include <LoadControl.h>
#include <Node.h>
#include <Matrix.h>
#include <Vector.h>
#include <elementAPI.h>
#include <iostream>



// Drilling DOF mode (must match OpenSees definition)
enum DrillingDOFMode {
    NO_DRILLING_DOF = 0,
    BASIC_DRILLING_DOF = 1,
    ENHANCED_DRILLING_DOF = 2
};

int main() {
    // Output
    StandardStream sserr;
    OPS_Stream* opserrPtr = &sserr;
    Domain theDomain;

    // Section
    double E = 200.0e9, nu = 0.3, h = 0.1;
    SectionForceDeformation* theSection = new ElasticMembranePlateSection(1, E, nu, h, 1.0);

    // Nodes
    theDomain.addNode(new Node(1, 6, 0.0, 0.0, 0.0));
    theDomain.addNode(new Node(2, 6, 1.0, 0.0, 0.0));
    theDomain.addNode(new Node(3, 6, 0.0, 1.0, 0.0));

    // Constraints
    for (int dof = 0; dof < 6; ++dof) {
        theDomain.addSP_Constraint(new SP_Constraint(1, dof, 0.0, true));
        theDomain.addSP_Constraint(new SP_Constraint(3, dof, 0.0, true));
    }

    // Element
  // Element
  Vector drillingStiffness(1);
  drillingStiffness(0) = 1.0; // Default drilling stiffness value

  ASDShellT3* theElement = new ASDShellT3(
      1,       // element tag
      1,       // node i
      2,       // node j
      3,       // node k
      theSection,  // section object
      drillingStiffness,  // drilling stiffness vector
      true,    // corotational transformation flag
      false,   // reduced integration flag
      //ASDShellT3::DrillingDOFMode::BASIC_DRILLING_DOF, // drilling mode
      ASDShellT3::DrillingDOF_Elastic,
      nullptr  // damping (optional)
  );
    theDomain.addElement(theElement);

    // Load pattern
    LoadPattern* thePattern = new LoadPattern(1);
    TimeSeries* theSeries = new LinearSeries();
    thePattern->setTimeSeries(theSeries);
    theDomain.addLoadPattern(thePattern);

    Vector load(6); load(2) = -50.0;
    NodalLoad* theLoad = new NodalLoad(1, 2, load, false);
    theDomain.addNodalLoad(theLoad, 1);

    // Static analysis setup
    AnalysisModel* theModel = new AnalysisModel();
    CTestNormDispIncr* theTest = new CTestNormDispIncr(1e-5, 100, 0);
    EquiSolnAlgo* theAlgo = new NewtonRaphson(*theTest);
    StaticIntegrator* theIntegrator = new LoadControl(1.0, 1, 1.0, 1.0);
    ConstraintHandler* theHandler = new TransformationConstraintHandler();
    DOF_Numberer* theNumberer = new DOF_Numberer(*new RCM());
    //LinearSOE* theSOE = new BandGenLinSOE(*new BandGenLinSolver());

    // Use BandSPDLinSolver instead of abstract BandGenLinSolver
    //LinearSOE* theSOE = new BandGenLinSOE(*new BandSPDLinSolver());

    //FullGenLinSolver    *theSolver = new FullGenLinLapackSolver();    
    //LinearSOE* theSOE = new BandSPDLinSOE(theSolver);
  
    #ifdef BLAS_EXISTS
    cout << "BLAS EXISTS"<<endl;
    #endif

/*
    StaticAnalysis theAnalysis(theDomain, *theHandler, *theNumberer, *theModel,
                             *theAlgo, *theSOE, *theIntegrator);
    // Run analysis
    int analysisResult = theAnalysis.analyze(1);
    if (analysisResult < 0) {
        std::cerr << "Analysis failed!" << std::endl;
        return -1;
    }
*/
    // Displacement at node 2
    const Vector& disp = theDomain.getNode(2)->getDisp();
    //std::cout << "Node 2 displacement: " << disp << std::endl;

    // Stiffness matrix
    Matrix K(18, 18);
    //theElement->getTangentStiff(K);
    //std::cout << "\nElement tangent stiffness matrix:\n" << K << std::endl;

    // Internal forces
    const Vector& R = theElement->getResistingForce();
    //std::cout << "\nInternal resisting force:\n" << R << std::endl;

    // Clean up
    delete theSection;
    delete theElement;
    delete thePattern;
    delete theSeries;
    delete theLoad;
    delete theModel;
    delete theTest;
    delete theAlgo;
    delete theIntegrator;
    delete theHandler;
    delete theNumberer;
    
    //delete theSOE;

    return 0;
}

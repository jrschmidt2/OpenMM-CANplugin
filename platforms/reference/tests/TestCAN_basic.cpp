#include "CANReferenceTest.h"

using namespace CustomAnisotropicNonbondedPlugin;
using namespace OpenMM;
using namespace std;
using namespace std::chrono;

//ReferencePlatform platformtest;

extern "C" OPENMM_EXPORT void registerCustomAnisotropicNonbondedReferenceKernelFactories();

const double TOL = 1e-10;

void test_novar() {
	cout << "TEST NOVAR" << endl;
	//system
	System system;
	system.addParticle(1.0);
	system.addParticle(1.0);
	//forceField
	CustomAnisotropicNonbondedForce* forceField = new CustomAnisotropicNonbondedForce("-0.1");
	forceField->addParticle(vector<double>(),5,-1,-1,-1);
	forceField->addParticle(vector<double>(),5,-1,-1,-1);
	system.addForce(forceField);
	//integrator
	VerletIntegrator integrator(0.01);
	//platform
	Platform& platform = Platform::getPlatformByName("Reference");
	//context & positions
	Context context(system, integrator, platform);
	vector<Vec3> positions(2);
	positions[0] = Vec3(0,0,0);	
	positions[1] = Vec3(0,0,1);	
	context.setPositions(positions);
	//state & calculate
	State state = context.getState(State::Forces | State::Energy);
	const vector<Vec3>& forces = state.getForces();
	double refforce = 0;
	ASSERT_EQUAL_VEC(Vec3(0, 0, -refforce), forces[0], TOL);
	ASSERT_EQUAL_VEC(Vec3(0, 0, refforce), forces[1], TOL);
	ASSERT_EQUAL_TOL((-0.1), state.getPotentialEnergy(), TOL);
}

void test_simple() {
	cout << "TEST SIMPLE, Z" << endl;
	//system
	System system;
	system.addParticle(1.0);
	system.addParticle(1.0);
	//forceField
	CustomAnisotropicNonbondedForce* forceField = new CustomAnisotropicNonbondedForce("-0.1*r");
	forceField->addParticle(vector<double>(),5,-1,-1,-1);
	forceField->addParticle(vector<double>(),5,-1,-1,-1);
	system.addForce(forceField);
	//integrator
	VerletIntegrator integrator(0.01);
	//platform
	Platform& platform = Platform::getPlatformByName("Reference");
	//context & positions
	Context context(system, integrator, platform);
	vector<Vec3> positions(2);
	positions[0] = Vec3(0,0,0);	
	positions[1] = Vec3(0,0,1);	
	context.setPositions(positions);
	//state & calculate
	State state = context.getState(State::Forces | State::Energy);
	const vector<Vec3>& forces = state.getForces();
	double refforce = 0.1;
	ASSERT_EQUAL_VEC(Vec3(0, 0, -refforce), forces[0], TOL);
	ASSERT_EQUAL_VEC(Vec3(0, 0, refforce), forces[1], TOL);
	ASSERT_EQUAL_TOL((-0.1), state.getPotentialEnergy(), TOL);
}

void test_der() {
	cout << "TEST DERIVATIVE, Z" << endl;
	//system
	System system;
	system.addParticle(1.0);
	system.addParticle(1.0);
	//forceField
	CustomAnisotropicNonbondedForce* forceField = new CustomAnisotropicNonbondedForce("-0.1/(r^3)");
	forceField->addParticle(vector<double>(),5,-1,-1,-1);
	forceField->addParticle(vector<double>(),5,-1,-1,-1);
	system.addForce(forceField);
	//integrator
	VerletIntegrator integrator(0.01);
	//platform
	Platform& platform = Platform::getPlatformByName("Reference");
	//context & positions
	Context context(system, integrator, platform);
	vector<Vec3> positions(2);
	positions[0] = Vec3(0,0,0);	
	positions[1] = Vec3(0,0,2);	
	context.setPositions(positions);
	//state & calculate
	State state = context.getState(State::Forces | State::Energy);
	const vector<Vec3>& forces = state.getForces();
	double refforce = -0.3/(2.0*2.0*2.0*2.0);
	ASSERT_EQUAL_VEC(Vec3(0, 0, -refforce), forces[0], TOL);
	ASSERT_EQUAL_VEC(Vec3(0, 0, refforce), forces[1], TOL);
	ASSERT_EQUAL_TOL((-0.1/8.0), state.getPotentialEnergy(), TOL);
}


void test_cubic() {
	cout << "TEST 3 PARTICLE" << endl;
	System system;
	system.addParticle(1.0);
	system.addParticle(1.0);
	system.addParticle(1.0);
	VerletIntegrator integrator(0.01);
	CustomAnisotropicNonbondedForce* forceField = new CustomAnisotropicNonbondedForce("-0.5*r^2");
	forceField->addParticle(vector<double>(),4,-1,-1,1);
	forceField->addParticle(vector<double>(),4,-1,-1,0);
	forceField->addParticle(vector<double>(),5,-1,-1,-1);
	system.addForce(forceField);
	Platform& platform = Platform::getPlatformByName("Reference");
	Context context(system, integrator, platform);
	vector<Vec3> positions(3);
	positions[0] = Vec3(0,0,0);
	positions[1] = Vec3(1,0,0);
	positions[2] = Vec3(0.5,sqrt(3)/2,0);
	context.setPositions(positions);
	State state = context.getState(State::Forces | State::Energy);
	const vector<Vec3>& forces = state.getForces();
	double sin3 = sin(M_PI/3);
	ASSERT_EQUAL_VEC(Vec3(-1.5, -sin3, 0), forces[0], TOL);
	ASSERT_EQUAL_VEC(Vec3(1.5, -sin3, 0), forces[1], TOL);
	ASSERT_EQUAL_VEC(Vec3(0, 2.0*sin3, 0), forces[2], TOL);
	ASSERT_EQUAL_TOL((-1.5), state.getPotentialEnergy(), TOL);
}

void test_global() {
	cout << "TEST GLOBAL PARAMETERS" << endl;
	System system;
	system.addParticle(1.0);
	system.addParticle(1.0);
	CustomAnisotropicNonbondedForce* forceField = new CustomAnisotropicNonbondedForce("-a*r");
	forceField->addParticle(vector<double>(),5,-1,-1,-1);
	forceField->addParticle(vector<double>(),5,-1,-1,-1);
	forceField->addGlobalParameter("a",0.1);
	system.addForce(forceField);
	VerletIntegrator integrator(0.01);
	Platform& platform = Platform::getPlatformByName("Reference");
	Context context(system, integrator, platform);
	vector<Vec3> positions(2);
	positions[0] = Vec3(0,0,0);	
	positions[1] = Vec3(1,0,0);	
	context.setPositions(positions);
	State state = context.getState(State::Forces | State::Energy);
	const vector<Vec3>& forces = state.getForces();
	double refforce = 0.1;
	ASSERT_EQUAL_VEC(Vec3(-refforce, 0, 0), forces[0], TOL);
	ASSERT_EQUAL_VEC(Vec3(refforce, 0, 0), forces[1], TOL);
	ASSERT_EQUAL_TOL((-0.1), state.getPotentialEnergy(), TOL);
}

void test_per() {
	cout << "TEST PER PARTICLE PARAMETERS" << endl;
	System system;
	system.addParticle(1.0);
	system.addParticle(1.0);
	CustomAnisotropicNonbondedForce* forceField = new CustomAnisotropicNonbondedForce("-a*r; a=a1*a2");
	forceField->addPerParticleParameter("a");
	vector<double> params(1);
	params[0] = 0.5;
	forceField->addParticle(params,5,-1,-1,-1);
	params[0] = 0.2;
	forceField->addParticle(params,5,-1,-1,-1);
	system.addForce(forceField);
	VerletIntegrator integrator(0.01);
	Platform& platform = Platform::getPlatformByName("Reference");
	Context context(system, integrator, platform);
	vector<Vec3> positions(2);
	positions[0] = Vec3(0,0,0);	
	positions[1] = Vec3(1,0,0);	
	context.setPositions(positions);
	State state = context.getState(State::Forces | State::Energy);
	const vector<Vec3>& forces = state.getForces();
	double refforce = 0.1;
	ASSERT_EQUAL_VEC(Vec3(-refforce, 0, 0), forces[0], TOL);
	ASSERT_EQUAL_VEC(Vec3(refforce, 0, 0), forces[1], TOL);
	ASSERT_EQUAL_TOL((-0.1), state.getPotentialEnergy(), TOL);
}

void test_theta() {
	cout << "TEST 3 PARTICLE THETA" << endl;
	System system;
	system.addParticle(1.0);
	system.addParticle(1.0);
	system.addParticle(1.0);
	VerletIntegrator integrator(0.01);
	CustomAnisotropicNonbondedForce* forceField = new CustomAnisotropicNonbondedForce("sin(theta1)");
	forceField->addParticle(vector<double>(),4,-1,-1,1);
	forceField->addParticle(vector<double>(),4,-1,-1,0);
	forceField->addParticle(vector<double>(),5,-1,-1,-1);
	system.addForce(forceField);
	Platform& platform = Platform::getPlatformByName("Reference");
	Context context(system, integrator, platform);
	vector<Vec3> positions(3);
	positions[0] = Vec3(0,0,0);
	positions[1] = Vec3(1,0,0);
	positions[2] = Vec3(0.5,sqrt(3)/2,0);
	context.setPositions(positions);
	State state = context.getState(State::Forces | State::Energy);
	const vector<Vec3>& forces = state.getForces();
	double sin3 = sin(M_PI/3);
	ASSERT_EQUAL_VEC(Vec3(-0.5*sin3,0.25,0), forces[0], TOL);
	ASSERT_EQUAL_VEC(Vec3(0.5*sin3,0.25,0), forces[1], TOL);
	ASSERT_EQUAL_VEC(Vec3(0,-0.5,0), forces[2], TOL);
	ASSERT_EQUAL_TOL(2*(sin3),state.getPotentialEnergy(),TOL);		
}


void test_phi() {
	cout << "TEST PHI" << endl;
	System system;
	system.addParticle(1.0);
	system.addParticle(1.0);
	VerletIntegrator integrator(0.01);
	CustomAnisotropicNonbondedForce* forceField = new CustomAnisotropicNonbondedForce("sin(phi1)");
	forceField->addParticle(vector<double>(),5,-1,-1,-1);
	forceField->addParticle(vector<double>(),5,-1,-1,-1);
	system.addForce(forceField);
	Platform& platform = Platform::getPlatformByName("Reference");
	Context context(system, integrator, platform);
	vector<Vec3> positions(2);
	positions[0] = Vec3(0,0,0);
	positions[1] = Vec3(1,1,0);
	context.setPositions(positions);
	State state = context.getState(State::Forces | State::Energy);
	const vector<Vec3>& forces = state.getForces();
	double sin4 = sin(M_PI/4.0);
	double cos4 = cos(M_PI/4.0);
	ASSERT_EQUAL_VEC(Vec3(-cos4/2.0,-cos4/2.0,0), forces[0], TOL);
	ASSERT_EQUAL_VEC(Vec3(cos4/2.0,-cos4/2.0,0),forces[1],TOL);
	ASSERT_EQUAL_TOL(sin4 ,state.getPotentialEnergy(),TOL);		
}

void test_theta2() {
	cout << "TEST 3 PARTICLE THETA1 & THETA2" << endl;
	System system;
	system.addParticle(1.0);
	system.addParticle(1.0);
	system.addParticle(1.0);
	VerletIntegrator integrator(0.01);
	CustomAnisotropicNonbondedForce* forceField = new CustomAnisotropicNonbondedForce("sin(theta1) + sin(theta2)");
	forceField->addParticle(vector<double>(),4,-1,-1,1);
	forceField->addParticle(vector<double>(),4,-1,-1,0);
	forceField->addParticle(vector<double>(),5,-1,-1,-1);
	system.addForce(forceField);
	Platform& platform = Platform::getPlatformByName("Reference");
	Context context(system, integrator, platform);
	vector<Vec3> positions(3);
	positions[0] = Vec3(0,0,0);
	positions[1] = Vec3(1,0,0);
	positions[2] = Vec3(0.5,sqrt(3.0)/2.0,0);
	context.setPositions(positions);
	State state = context.getState(State::Forces | State::Energy);
	const vector<Vec3>& forces = state.getForces();
	double sin3 = sin(M_PI/3.0);
	ASSERT_EQUAL_VEC(Vec3(-sin3*0.5,0.25, 0), forces[0], TOL);
	ASSERT_EQUAL_VEC(Vec3(sin3*0.5,0.25, 0), forces[1], TOL);
	ASSERT_EQUAL_VEC(Vec3(0,-0.5, 0), forces[2], TOL);
	ASSERT_EQUAL_TOL(2.0*sin3+2.0,state.getPotentialEnergy(),TOL);		
}
void test_rtheta() {
	cout << "TEST 4 PARTICLE R,THETA1" << endl;
	System system;
	system.addParticle(1.0);
	system.addParticle(1.0);
	system.addParticle(1.0);
	system.addParticle(1.0);
	VerletIntegrator integrator(0.01);
	CustomAnisotropicNonbondedForce* forceField = new CustomAnisotropicNonbondedForce("r*sin(theta1)");
	forceField->addParticle(vector<double>(),4,-1,-1,1);
	forceField->addParticle(vector<double>(),4,-1,-1,0);
	forceField->addParticle(vector<double>(),4,-1,-1,3);
	forceField->addParticle(vector<double>(),4,-1,-1,2);
	system.addForce(forceField);
	Platform& platform = Platform::getPlatformByName("Reference");
	Context context(system, integrator, platform);
	vector<Vec3> positions(4);
	positions[0] = Vec3(0,0,0);
	positions[1] = Vec3(1,0,0);
	positions[2] = Vec3(0.5,1,-0.5);
	positions[3] = Vec3(0.5,1,0.5);
	context.setPositions(positions);
	State state = context.getState(State::Forces | State::Energy);
	const vector<Vec3>& forces = state.getForces();
	ASSERT_EQUAL_VEC(Vec3(0,4.0/sqrt(5.0),0), forces[0], TOL);
	ASSERT_EQUAL_VEC(Vec3(0,4.0/sqrt(5.0),0), forces[1], TOL);
	ASSERT_EQUAL_VEC(Vec3(0,-4.0/sqrt(5.0),2.0/sqrt(5.0)), forces[2], TOL);
	ASSERT_EQUAL_VEC(Vec3(0,-4.0/sqrt(5.0),-2.0/sqrt(5.0)), forces[3], TOL);
	ASSERT_EQUAL_TOL(sqrt(1.5)*sqrt(5.0/6.0)*4.0,state.getPotentialEnergy(),TOL);		
}

void testExclusions() {
	cout << "TEST EXCLUSIONS" << endl;
    System system;
    VerletIntegrator integrator(0.01);
    CustomAnisotropicNonbondedForce* nonbonded = new CustomAnisotropicNonbondedForce("a*r; a=a1+a2");
    nonbonded->addPerParticleParameter("a");
    vector<double> params(1);
    vector<Vec3> positions(4);
    for (int i = 0; i < 4; i++) {
        system.addParticle(1.0);
        params[0] = i+1;
        nonbonded->addParticle(params,5,-1,-1,-1);
        positions[i] = Vec3(i, 0, 0);
    }
    nonbonded->addExclusion(0, 1);
    nonbonded->addExclusion(1, 2);
    nonbonded->addExclusion(2, 3);
    nonbonded->addExclusion(0, 2);
    nonbonded->addExclusion(1, 3);
    system.addForce(nonbonded);
    Platform& platform = Platform::getPlatformByName("Reference");
    Context context(system, integrator, platform);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forceFields = state.getForces();
    ASSERT_EQUAL_VEC(Vec3(1+4, 0, 0), forceFields[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, 0, 0), forceFields[1], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, 0, 0), forceFields[2], TOL);
    ASSERT_EQUAL_VEC(Vec3(-(1+4), 0, 0), forceFields[3], TOL);
    ASSERT_EQUAL_TOL((1+4)*3.0, state.getPotentialEnergy(), TOL);
}

void testCutoff() {
	cout << "TEST CUTOFF" << endl;
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomAnisotropicNonbondedForce* forceFieldField = new CustomAnisotropicNonbondedForce("r");
    forceFieldField->addParticle(vector<double>(),5,-1,-1,-1);
    forceFieldField->addParticle(vector<double>(),5,-1,-1,-1);
    forceFieldField->addParticle(vector<double>(),5,-1,-1,-1);
    forceFieldField->setNonbondedMethod(CustomAnisotropicNonbondedForce::CutoffNonPeriodic);
    forceFieldField->setCutoffDistance(2.5);
    system.addForce(forceFieldField);
    ASSERT(!forceFieldField->usesPeriodicBoundaryConditions());
    ASSERT(!system.usesPeriodicBoundaryConditions());
    Platform& platform = Platform::getPlatformByName("Reference");
    Context context(system, integrator, platform);
    vector<Vec3> positions(3);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(0, 2, 0);
    positions[2] = Vec3(0, 3, 0);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forceFields = state.getForces();
    ASSERT_EQUAL_VEC(Vec3(0, 1, 0), forceFields[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, 0, 0), forceFields[1], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, -1, 0), forceFields[2], TOL);
    ASSERT_EQUAL_TOL(2.0+1.0, state.getPotentialEnergy(), TOL);
}
void testManyParameters() {
	System system;
	system.addParticle(1.0);
	system.addParticle(1.0);
	VerletIntegrator integrator(0.01);
	CustomAnisotropicNonbondedForce* forceField = new CustomAnisotropicNonbondedForce("(a1*a2+b1*b2+c1*c2+d1*d2+e1*e2)*r");
	forceField->addPerParticleParameter("a");
	forceField->addPerParticleParameter("b");
	forceField->addPerParticleParameter("c");
	forceField->addPerParticleParameter("d");
	forceField->addPerParticleParameter("e");
	vector<double> params(5);
	params[0] = 1.0;
	params[1] = 2.0;
	params[2] = 3.0;
	params[3] = 4.0;
	params[4] = 5.0;
	forceField->addParticle(params,5,-1,-1,-1);
	params[0] = 1.1;
	params[1] = 1.2;
	params[2] = 1.3;
	params[3] = 1.4;
	params[4] = 1.5;
	forceField->addParticle(params,5,-1,-1,-1);
	system.addForce(forceField);
	Platform& platform = Platform::getPlatformByName("Reference");
	Context context(system,integrator,platform);
	vector<Vec3> positions(2);
	positions[0] = Vec3(0,0,0);
	positions[1] = Vec3(2,0,0);
	context.setPositions(positions);
	State state = context.getState(State::Forces | State::Energy);
	vector<Vec3> forces = state.getForces();
	double force = 1*1.1 + 2*1.2 + 3*1.3 + 4*1.4 + 5*1.5;
	ASSERT_EQUAL_VEC(Vec3(force,0,0), forces[0], TOL);
	ASSERT_EQUAL_VEC(Vec3(-force,0,0), forces[1], TOL);
	ASSERT_EQUAL_TOL(2*force, state.getPotentialEnergy(), TOL);
}

void testPeriodic() {
	cout << "TEST PERIODIC" << endl;
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    CustomAnisotropicNonbondedForce* forceFieldField = new CustomAnisotropicNonbondedForce("r");
    forceFieldField->addParticle(vector<double>(),5, -1,-1,-1);
    forceFieldField->addParticle(vector<double>(),5,-1,-1,-1);
    forceFieldField->addParticle(vector<double>(),5,-1,-1,-1);
    forceFieldField->setNonbondedMethod(CustomAnisotropicNonbondedForce::CutoffPeriodic);
    forceFieldField->setCutoffDistance(2.0);
    system.setDefaultPeriodicBoxVectors(Vec3(4, 0, 0), Vec3(0, 4, 0), Vec3(0, 0, 4));
    system.addForce(forceFieldField);
    ASSERT(forceFieldField->usesPeriodicBoundaryConditions());
    ASSERT(system.usesPeriodicBoundaryConditions());
    Platform& platform = Platform::getPlatformByName("Reference");
    Context context(system, integrator, platform);
    vector<Vec3> positions(3);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(0, 2.1, 0);
    positions[2] = Vec3(0, 3, 0);
    context.setPositions(positions);
    State state = context.getState(State::Forces | State::Energy);
    const vector<Vec3>& forceFields = state.getForces();
    ASSERT_EQUAL_VEC(Vec3(0, -2, 0), forceFields[0], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, 2, 0), forceFields[1], TOL);
    ASSERT_EQUAL_VEC(Vec3(0, 0, 0), forceFields[2], TOL);
    ASSERT_EQUAL_TOL(1.9+1+0.9, state.getPotentialEnergy(), TOL);
}

void testMultipleCutoffs() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    VerletIntegrator integrator(0.01);
    
    // Add multiple nonbonded forceFields that have different cutoffs.
    
    CustomAnisotropicNonbondedForce* nonbonded1 = new CustomAnisotropicNonbondedForce("2*r");
    nonbonded1->addParticle(vector<double>(),5,-1,-1,-1);
    nonbonded1->addParticle(vector<double>(),5,-1,-1,-1);
    nonbonded1->setNonbondedMethod(CustomAnisotropicNonbondedForce::CutoffNonPeriodic);
    nonbonded1->setCutoffDistance(2.5);
    system.addForce(nonbonded1);
    CustomAnisotropicNonbondedForce* nonbonded2 = new CustomAnisotropicNonbondedForce("3*r");
    nonbonded2->addParticle(vector<double>(),5,-1,-1,-1);
    nonbonded2->addParticle(vector<double>(),5,-1,-1,-1);
    nonbonded2->setNonbondedMethod(CustomAnisotropicNonbondedForce::CutoffNonPeriodic);
    nonbonded2->setCutoffDistance(2.9);
    nonbonded2->setForceGroup(1);
    system.addForce(nonbonded2);
    Platform& platform = Platform::getPlatformByName("Reference");
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    positions[0] = Vec3(0, 0, 0);
    positions[1] = Vec3(0, 0, 0);
    for (double r = 2.4; r < 3.2; r += 0.2) {
        positions[1][1] = r;
        context.setPositions(positions);
        double e1 = (r < 2.5 ? 2.0*r : 0.0);
        double e2 = (r < 2.9 ? 3.0*r : 0.0);
        double f1 = (r < 2.5 ? 2.0 : 0.0);
        double f2 = (r < 2.9 ? 3.0 : 0.0);
        
        // Check the first forceField.
        
        State state = context.getState(State::Forces | State::Energy, false, 1);
        ASSERT_EQUAL_VEC(Vec3(0, f1, 0), state.getForces()[0], TOL);
        ASSERT_EQUAL_VEC(Vec3(0, -f1, 0), state.getForces()[1], TOL);
        ASSERT_EQUAL_TOL(e1, state.getPotentialEnergy(), TOL);
        
        // Check the second forceField.
        
        state = context.getState(State::Forces | State::Energy, false, 2);
        ASSERT_EQUAL_VEC(Vec3(0, f2, 0), state.getForces()[0], TOL);
        ASSERT_EQUAL_VEC(Vec3(0, -f2, 0), state.getForces()[1], TOL);
        ASSERT_EQUAL_TOL(e2, state.getPotentialEnergy(), TOL);
        
        // Check the sum of both forceFields.

        state = context.getState(State::Forces | State::Energy);
        ASSERT_EQUAL_VEC(Vec3(0, f1+f2, 0), state.getForces()[0], TOL);
        ASSERT_EQUAL_VEC(Vec3(0, -f1-f2, 0), state.getForces()[1], TOL);
        ASSERT_EQUAL_TOL(e1+e2, state.getPotentialEnergy(), TOL);
    }
}

void testTriclinic() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    Vec3 a(3.1, 0, 0);
    Vec3 b(0.4, 3.5, 0);
    Vec3 c(-0.1, -0.5, 4.0);
    system.setDefaultPeriodicBoxVectors(a, b, c);
    VerletIntegrator integrator(0.01);
    CustomAnisotropicNonbondedForce* nonbonded = new CustomAnisotropicNonbondedForce("r");
    nonbonded->addParticle(vector<double>(),5,-1,-1,-1);
    nonbonded->addParticle(vector<double>(),5,-1,-1,-1);
    nonbonded->setNonbondedMethod(CustomAnisotropicNonbondedForce::CutoffPeriodic);
    const double cutoff = 1.5;
    nonbonded->setCutoffDistance(cutoff);
    system.addForce(nonbonded);
    Platform& platform = Platform::getPlatformByName("Reference");
    Context context(system, integrator, platform);
    vector<Vec3> positions(2);
    OpenMM_SFMT::SFMT sfmt;
    init_gen_rand(0, sfmt);
    for (int iteration = 0; iteration < 50; iteration++) {
        // Generate random positions for the two particles.

        positions[0] = a*genrand_real2(sfmt) + b*genrand_real2(sfmt) + c*genrand_real2(sfmt);
        positions[1] = a*genrand_real2(sfmt) + b*genrand_real2(sfmt) + c*genrand_real2(sfmt);
        context.setPositions(positions);

        // Loop over all possible periodic copies and find the nearest one.

        Vec3 delta;
        double distance2 = 100.0;
        for (int i = -1; i < 2; i++)
            for (int j = -1; j < 2; j++)
                for (int k = -1; k < 2; k++) {
                    Vec3 d = positions[1]-positions[0]+a*i+b*j+c*k;
                    if (d.dot(d) < distance2) {
                        delta = d;
                        distance2 = d.dot(d);
                    }
                }
        double distance = sqrt(distance2);

        // See if the forceField and energy are correct.

        State state = context.getState(State::Forces | State::Energy);
        if (distance >= cutoff) {
            ASSERT_EQUAL(0.0, state.getPotentialEnergy());
            ASSERT_EQUAL_VEC(Vec3(0, 0, 0), state.getForces()[0], 0);
            ASSERT_EQUAL_VEC(Vec3(0, 0, 0), state.getForces()[1], 0);
        }
        else {
            const Vec3 forceField = delta/sqrt(delta.dot(delta));
            ASSERT_EQUAL_TOL(distance, state.getPotentialEnergy(), TOL);
            ASSERT_EQUAL_VEC(forceField, state.getForces()[0], TOL);
            ASSERT_EQUAL_VEC(-forceField, state.getForces()[1], TOL);
        }
    }
}

int main() {
   steady_clock::time_point start = steady_clock::now();
    try {
	registerCustomAnisotropicNonbondedReferenceKernelFactories();
	test_novar();
	test_simple();
	test_der();
	test_cubic();
	test_global();
	test_per();
	test_theta();
	test_phi();
	test_theta2();
	test_rtheta();
	testExclusions();
	testCutoff();
	testManyParameters();
	testPeriodic();
	testMultipleCutoffs();
	testTriclinic();
    }
    catch(const exception& e) {
        std::cout << "exception: " << e.what() << endl;
        return 1;
    }
    steady_clock::time_point end = steady_clock::now();
    duration<double> timer = duration_cast<duration<double>>(end-start);
    std::cout << timer.count() << endl;
    std::cout << "Done" << endl;
    return 0;
}



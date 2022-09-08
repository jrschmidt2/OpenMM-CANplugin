#ifdef WIN32
  #define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include "sfmt/SFMT.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/Context.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "CustomAnisotropicNonbondedForce.h"
#include "openmm/CustomHbondForce.h"
#include "openmm/CustomNonbondedForce.h"
#include <cmath>
#include <iostream>
#include <set>
#include <vector>
#include <string>
#include <cstring>
#include <ctime>
#include <chrono>

using namespace CustomAnisotropicNonbondedPlugin;
using namespace OpenMM;
using namespace std;
using namespace std::chrono;

extern "C" OPENMM_EXPORT void registerCustomAnisotropicNonbondedReferenceKernelFactories();

const double TOL = 1e-15;

void test_Y00() {
	cout  << "Test YOO" << endl;
	System system;
	system.addParticle(1.0);
	system.addParticle(1.0);
	system.addParticle(1.0);
	VerletIntegrator integrator(0.01);
	CustomAnisotropicNonbondedForce* forceField = new CustomAnisotropicNonbondedForce("y00");
	forceField->addParticle(vector<double>(),4,-1,-1,1);
	forceField->addParticle(vector<double>(),4,-1,-1,2);
	forceField->addParticle(vector<double>(),4,-1,-1,0);
	system.addForce(forceField);
	Platform& platform = Platform::getPlatformByName("Reference");
	Context context(system, integrator, platform);
	vector<Vec3> positions(3);
	positions[0] = Vec3(0,0,0);
	positions[1] = Vec3(0,1,0);
	positions[2] = Vec3(sqrt(3)/2,0.5,0);
	context.setPositions(positions);
	State state = context.getState(State::Forces | State::Energy);
	const vector<Vec3>& forces = state.getForces();

	System refsystem;
	refsystem.addParticle(1.0);
	refsystem.addParticle(1.0);
	refsystem.addParticle(1.0);
	VerletIntegrator refintegrator(0.01);
	CustomAnisotropicNonbondedForce* refFF = new CustomAnisotropicNonbondedForce("sqrt(1/(4*pi))");
	refFF->addGlobalParameter("pi",M_PI);
	refFF->addParticle(vector<double>(),4,-1,-1,1);
	refFF->addParticle(vector<double>(),4,-1,-1,2);
	refFF->addParticle(vector<double>(),4,-1,-1,0);
	refsystem.addForce(refFF);
	Platform& refplatform = Platform::getPlatformByName("Reference");
	Context refcontext(refsystem,refintegrator,refplatform);
	refcontext.setPositions(positions);
	State refstate = refcontext.getState(State::Forces | State::Energy);
	const vector<Vec3>& refforces = refstate.getForces();
	ASSERT_EQUAL_VEC(refforces[0],forces[0],TOL);
	ASSERT_EQUAL_VEC(refforces[1],forces[1],TOL);
	ASSERT_EQUAL_VEC(refforces[2],forces[2],TOL);
	ASSERT_EQUAL_TOL(refstate.getPotentialEnergy(), state.getPotentialEnergy(), TOL);
}

void test_Y10() {
	cout  << "Test Y1O" << endl;
	System system;
	system.addParticle(1.0);
	system.addParticle(1.0);
	system.addParticle(1.0);
	VerletIntegrator integrator(0.01);
	CustomAnisotropicNonbondedForce* forceField = new CustomAnisotropicNonbondedForce("y10_1+y10_2");
	forceField->addParticle(vector<double>(),4,-1,-1,1);
	forceField->addParticle(vector<double>(),4,-1,-1,2);
	forceField->addParticle(vector<double>(),4,-1,-1,0);
	system.addForce(forceField);
	Platform& platform = Platform::getPlatformByName("Reference");
	Context context(system, integrator, platform);
	vector<Vec3> positions(3);
	positions[0] = Vec3(0,0,0);
	positions[1] = Vec3(0,1,0);
	positions[2] = Vec3(sqrt(3)/2,0.5,0);
	context.setPositions(positions);
	State state = context.getState(State::Forces | State::Energy);
	const vector<Vec3>& forces = state.getForces();

	System refsystem;
	refsystem.addParticle(1.0);
	refsystem.addParticle(1.0);
	refsystem.addParticle(1.0);
	VerletIntegrator refintegrator(0.01);
	CustomAnisotropicNonbondedForce* refFF = new CustomAnisotropicNonbondedForce("sqrt(3/(4*pi))*(cos(theta1)+cos(theta2))");
	refFF->addGlobalParameter("pi",M_PI);
	refFF->addParticle(vector<double>(),4,-1,-1,1);
	refFF->addParticle(vector<double>(),4,-1,-1,2);
	refFF->addParticle(vector<double>(),4,-1,-1,0);
	refsystem.addForce(refFF);
	Platform& refplatform = Platform::getPlatformByName("Reference");
	Context refcontext(refsystem,refintegrator,refplatform);
	refcontext.setPositions(positions);
	State refstate = refcontext.getState(State::Forces | State::Energy);
	const vector<Vec3>& refforces = refstate.getForces();
	ASSERT_EQUAL_VEC(refforces[0],forces[0],TOL);
	ASSERT_EQUAL_VEC(refforces[1],forces[1],TOL);
	ASSERT_EQUAL_VEC(refforces[2],forces[2],TOL);
	ASSERT_EQUAL_TOL(refstate.getPotentialEnergy(), state.getPotentialEnergy(), TOL);
}

void test_Y11() {
	cout  << "Test Y11" << endl;
	System system;
	system.addParticle(1.0);
	system.addParticle(1.0);
	system.addParticle(1.0);
	VerletIntegrator integrator(0.01);
	CustomAnisotropicNonbondedForce* forceField = new CustomAnisotropicNonbondedForce("y11c_1+y11c_2+y11s_1+y11s_2");
	forceField->addGlobalParameter("pi",M_PI);
	forceField->addParticle(vector<double>(),4,-1,-1,1);
	forceField->addParticle(vector<double>(),4,-1,-1,2);
	forceField->addParticle(vector<double>(),4,-1,-1,0);
	system.addForce(forceField);
	Platform& platform = Platform::getPlatformByName("Reference");
	Context context(system, integrator, platform);
	vector<Vec3> positions(3);
	positions[0] = Vec3(0,0,0);
	positions[1] = Vec3(0,1,0);
	positions[2] = Vec3(sqrt(3)/2,0.5,0);
	context.setPositions(positions);
	State state = context.getState(State::Forces | State::Energy);
	const vector<Vec3>& forces = state.getForces();

	System refsystem;
	refsystem.addParticle(1.0);
	refsystem.addParticle(1.0);
	refsystem.addParticle(1.0);
	VerletIntegrator refintegrator(0.01);
	CustomAnisotropicNonbondedForce* refFF = new CustomAnisotropicNonbondedForce("sqrt(3/(4*pi))*(sin(theta1)*(cos(phi1)+sin(phi1)) + sin(theta2)*(cos(phi2)+sin(phi2)))");
	refFF->addGlobalParameter("pi",M_PI);
	refFF->addParticle(vector<double>(),4,-1,-1,1);
	refFF->addParticle(vector<double>(),4,-1,-1,2);
	refFF->addParticle(vector<double>(),4,-1,-1,0);
	refsystem.addForce(refFF);
	Platform& refplatform = Platform::getPlatformByName("Reference");
	Context refcontext(refsystem,refintegrator,refplatform);
	refcontext.setPositions(positions);
	State refstate = refcontext.getState(State::Forces | State::Energy);
	const vector<Vec3>& refforces = refstate.getForces();
	ASSERT_EQUAL_VEC(refforces[0],forces[0],TOL);
	ASSERT_EQUAL_VEC(refforces[1],forces[1],TOL);
	ASSERT_EQUAL_VEC(refforces[2],forces[2],TOL);
	ASSERT_EQUAL_TOL(refstate.getPotentialEnergy(), state.getPotentialEnergy(), TOL);
}

void test_Y20() {
	cout  << "Test Y20" << endl;
	System system;
	system.addParticle(1.0);
	system.addParticle(1.0);
	system.addParticle(1.0);
	VerletIntegrator integrator(0.01);
	CustomAnisotropicNonbondedForce* forceField = new CustomAnisotropicNonbondedForce("y20_1+y20_2");
	forceField->addGlobalParameter("pi",M_PI);
	forceField->addParticle(vector<double>(),4,-1,-1,1);
	forceField->addParticle(vector<double>(),4,-1,-1,2);
	forceField->addParticle(vector<double>(),4,-1,-1,0);
	system.addForce(forceField);
	Platform& platform = Platform::getPlatformByName("Reference");
	Context context(system, integrator, platform);
	vector<Vec3> positions(3);
	positions[0] = Vec3(0,0,0);
	positions[1] = Vec3(0,1,0);
	positions[2] = Vec3(sqrt(3)/2,0.5,0);
	context.setPositions(positions);
	State state = context.getState(State::Forces | State::Energy);
	const vector<Vec3>& forces = state.getForces();

	System refsystem;
	refsystem.addParticle(1.0);
	refsystem.addParticle(1.0);
	refsystem.addParticle(1.0);
	VerletIntegrator refintegrator(0.01);
	CustomAnisotropicNonbondedForce* refFF = new CustomAnisotropicNonbondedForce("sqrt(5/(16*pi))*((3*cos(theta1)*cos(theta1)-1)+(3*cos(theta2)*cos(theta2)-1))");
	refFF->addGlobalParameter("pi",M_PI);
	refFF->addParticle(vector<double>(),4,-1,-1,1);
	refFF->addParticle(vector<double>(),4,-1,-1,2);
	refFF->addParticle(vector<double>(),4,-1,-1,0);
	refsystem.addForce(refFF);
	Platform& refplatform = Platform::getPlatformByName("Reference");
	Context refcontext(refsystem,refintegrator,refplatform);
	refcontext.setPositions(positions);
	State refstate = refcontext.getState(State::Forces | State::Energy);
	const vector<Vec3>& refforces = refstate.getForces();
	ASSERT_EQUAL_VEC(refforces[0],forces[0],TOL);
	ASSERT_EQUAL_VEC(refforces[1],forces[1],TOL);
	ASSERT_EQUAL_VEC(refforces[2],forces[2],TOL);
	ASSERT_EQUAL_TOL(refstate.getPotentialEnergy(), state.getPotentialEnergy(), TOL);
}
	
void test_Y21() {
	cout  << "Test Y21" << endl;
	System system;
	system.addParticle(1.0);
	system.addParticle(1.0);
	system.addParticle(1.0);
	VerletIntegrator integrator(0.01);
	CustomAnisotropicNonbondedForce* forceField = new CustomAnisotropicNonbondedForce("y21c_1+y21c_2+y21s_1+y21s_2");
	forceField->addGlobalParameter("pi",M_PI);
	forceField->addParticle(vector<double>(),4,-1,-1,1);
	forceField->addParticle(vector<double>(),4,-1,-1,2);
	forceField->addParticle(vector<double>(),4,-1,-1,0);
	system.addForce(forceField);
	Platform& platform = Platform::getPlatformByName("Reference");
	Context context(system, integrator, platform);
	vector<Vec3> positions(3);
	positions[0] = Vec3(0,0,0);
	positions[1] = Vec3(0,1,0);
	positions[2] = Vec3(sqrt(3)/2,0.5,0);
	context.setPositions(positions);
	State state = context.getState(State::Forces | State::Energy);
	const vector<Vec3>& forces = state.getForces();

	System refsystem;
	refsystem.addParticle(1.0);
	refsystem.addParticle(1.0);
	refsystem.addParticle(1.0);
	VerletIntegrator refintegrator(0.01);
	CustomAnisotropicNonbondedForce* refFF = new CustomAnisotropicNonbondedForce("sqrt(15/(4*pi))*(sin(theta1)*cos(theta1)*(cos(phi1)+sin(phi1)) + sin(theta2)*cos(theta2)*(cos(phi2)+sin(phi2)))");
	refFF->addGlobalParameter("pi",M_PI);
	refFF->addParticle(vector<double>(),4,-1,-1,1);
	refFF->addParticle(vector<double>(),4,-1,-1,2);
	refFF->addParticle(vector<double>(),4,-1,-1,0);
	refsystem.addForce(refFF);
	Platform& refplatform = Platform::getPlatformByName("Reference");
	Context refcontext(refsystem,refintegrator,refplatform);
	refcontext.setPositions(positions);
	State refstate = refcontext.getState(State::Forces | State::Energy);
	const vector<Vec3>& refforces = refstate.getForces();
	ASSERT_EQUAL_VEC(refforces[0],forces[0],TOL);
	ASSERT_EQUAL_VEC(refforces[1],forces[1],TOL);
	ASSERT_EQUAL_VEC(refforces[2],forces[2],TOL);
	ASSERT_EQUAL_TOL(refstate.getPotentialEnergy(), state.getPotentialEnergy(), TOL);
}

void test_Y22() {
	cout  << "Test Y22" << endl;
	System system;
	system.addParticle(1.0);
	system.addParticle(1.0);
	system.addParticle(1.0);
	VerletIntegrator integrator(0.01);
	CustomAnisotropicNonbondedForce* forceField = new CustomAnisotropicNonbondedForce("y22c_1+y22c_2+y22s_1+y22s_2");
	forceField->addGlobalParameter("pi",M_PI);
	forceField->addParticle(vector<double>(),4,-1,-1,1);
	forceField->addParticle(vector<double>(),4,-1,-1,2);
	forceField->addParticle(vector<double>(),4,-1,-1,0);
	system.addForce(forceField);
	Platform& platform = Platform::getPlatformByName("Reference");
	Context context(system, integrator, platform);
	vector<Vec3> positions(3);
	positions[0] = Vec3(0,0,0);
	positions[1] = Vec3(0,1,0);
	positions[2] = Vec3(sqrt(3)/2,0.5,0);
	context.setPositions(positions);
	State state = context.getState(State::Forces | State::Energy);
	const vector<Vec3>& forces = state.getForces();

	System refsystem;
	refsystem.addParticle(1.0);
	refsystem.addParticle(1.0);
	refsystem.addParticle(1.0);
	VerletIntegrator refintegrator(0.01);
	CustomAnisotropicNonbondedForce* refFF = new CustomAnisotropicNonbondedForce("sqrt(15/(16*pi))*(sin(theta1)*sin(theta1)*(cos(2*phi1)+sin(2*phi1)) + sin(theta2)*sin(theta2)*(cos(2*phi2)+sin(2*phi2)))");
	refFF->addGlobalParameter("pi",M_PI);
	refFF->addParticle(vector<double>(),4,-1,-1,1);
	refFF->addParticle(vector<double>(),4,-1,-1,2);
	refFF->addParticle(vector<double>(),4,-1,-1,0);
	refsystem.addForce(refFF);
	Platform& refplatform = Platform::getPlatformByName("Reference");
	Context refcontext(refsystem,refintegrator,refplatform);
	refcontext.setPositions(positions);
	State refstate = refcontext.getState(State::Forces | State::Energy);
	const vector<Vec3>& refforces = refstate.getForces();
	ASSERT_EQUAL_VEC(refforces[0],forces[0],TOL);
	ASSERT_EQUAL_VEC(refforces[1],forces[1],TOL);
	ASSERT_EQUAL_VEC(refforces[2],forces[2],TOL);
	ASSERT_EQUAL_TOL(refstate.getPotentialEnergy(), state.getPotentialEnergy(), TOL);
}

int main() {
   steady_clock::time_point start = steady_clock::now();
    try {
	registerCustomAnisotropicNonbondedReferenceKernelFactories();
	test_Y00();
	test_Y10();
	test_Y11();
	test_Y20();
	test_Y21();
	test_Y22();
    }
    catch(const exception& e) {
	std::cout << "exception: "<< e.what()<< endl;
	return 1;
    }
    steady_clock::time_point end = steady_clock::now();
    duration<double> timer = duration_cast<duration<double>>(end-start);
    std::cout << timer.count() << endl;
    std::cout << "Done" << endl;
    return 0;
}


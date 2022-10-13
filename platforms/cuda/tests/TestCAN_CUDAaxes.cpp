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
#include <chrono>

using namespace CustomAnisotropicNonbondedPlugin;
using namespace OpenMM;
using namespace std;
using namespace std::chrono;

//CUDAPlatform platformtest;

extern "C" OPENMM_EXPORT void registerCustomAnisotropicNonbondedCudaKernelFactories();

const double TOL = 1e-6;

void test_thetaZThenX() {
	cout << "TEST THETA ZThenX" << endl;
	vector<Vec3> positions(6);
	positions[0] = Vec3(0,0,0);
	positions[1] = Vec3(1,0,0);
	positions[2] = Vec3(0.5,sqrt(3)/2,0);
	positions[3] = Vec3(2,0,2);
	positions[4] = Vec3(2,-sqrt(3)/2,1);
	positions[5] = Vec3(2,0,0);
	System system;
	system.addParticle(1.0);
	system.addParticle(1.0);
	system.addParticle(1.0);
	system.addParticle(1.0);
	system.addParticle(1.0);
	system.addParticle(1.0);
	VerletIntegrator integrator(0.01);
	CustomAnisotropicNonbondedForce* forceField = new CustomAnisotropicNonbondedForce("(theta2+theta2)");
	forceField->addParticle(vector<double>(),0,1,-1,2);
	forceField->addParticle(vector<double>(),0,2,-1,0);
	forceField->addParticle(vector<double>(),0,0,-1,1);
	forceField->addParticle(vector<double>(),0,4,-1,5);
	forceField->addParticle(vector<double>(),0,5,-1,3);
	forceField->addParticle(vector<double>(),0,2,-1,4);
	forceField->addExclusion(0,1);
	forceField->addExclusion(0,2);
	forceField->addExclusion(2,1);
	forceField->addExclusion(3,4);
	forceField->addExclusion(3,5);
	forceField->addExclusion(4,4);
	system.addForce(forceField);
	Platform& platform = Platform::getPlatformByName("CUDA");
	Context context(system, integrator, platform);
	context.setPositions(positions);
	State state = context.getState(State::Forces | State::Energy);
	const vector<Vec3>& forces = state.getForces();
	
	cout << state.getPotentialEnergy() << endl;
}

void test_phiZThenX() {
	cout << "TEST PHI ZThenX" << endl;
	vector<Vec3> positions(6);
	positions[0] = Vec3(0,0,0);
	positions[1] = Vec3(1,0,0);
	positions[2] = Vec3(0.5,sqrt(3)/2,0);
	positions[3] = Vec3(2,0,2);
	positions[4] = Vec3(2,-sqrt(3)/2,1);
	positions[5] = Vec3(2,0,0);

	System system;
	system.addParticle(1.0);
	system.addParticle(1.0);
	system.addParticle(1.0);
	system.addParticle(1.0);
	system.addParticle(1.0);
	system.addParticle(1.0);
	VerletIntegrator integrator(0.01);
	CustomAnisotropicNonbondedForce* forceField = new CustomAnisotropicNonbondedForce("(phi1+phi2)");
	forceField->addParticle(vector<double>(),0,1,-1,2);
	forceField->addParticle(vector<double>(),0,2,-1,0);
	forceField->addParticle(vector<double>(),0,0,-1,1);
	forceField->addParticle(vector<double>(),0,4,-1,5);
	forceField->addParticle(vector<double>(),0,5,-1,3);
	forceField->addParticle(vector<double>(),0,2,-1,4);
	forceField->addExclusion(0,1);
	forceField->addExclusion(0,2);
	forceField->addExclusion(2,1);
	forceField->addExclusion(3,4);
	forceField->addExclusion(3,5);
	forceField->addExclusion(4,4);
	system.addForce(forceField);
	Platform& platform = Platform::getPlatformByName("CUDA");
	Context context(system, integrator, platform);
	context.setPositions(positions);
	State state = context.getState(State::Forces | State::Energy);
	const vector<Vec3>& forces = state.getForces();
	
	cout << state.getPotentialEnergy() << endl;

}

int main() {
    try {
	registerCustomAnisotropicNonbondedCudaKernelFactories();
	test_thetaZThenX();
	test_phiZThenX();
    }
    catch(const exception& e) {
        return 1;
    }
    return 0;
}



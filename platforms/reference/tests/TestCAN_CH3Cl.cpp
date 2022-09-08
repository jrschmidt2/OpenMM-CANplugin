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
//#include "openmm/CustomNonbondedForce.h"
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

const double TOL = 1e-10;

void test_CH3Cl() {
	cout << "TEST CH3Cl" << endl;
	//Define all constants
	std::vector<double> C_params(25);
	C_params[0]=1.555722e+02;
	C_params[1]=0.000000e+00;
	C_params[2]=0.000000e+00;
	C_params[3]=0.000000e+00;
	C_params[4]=9.604317e+01;
	C_params[5]=0.000000e+00;
	C_params[6]=0.000000e+00;
	C_params[7]=0.000000e+00;
	C_params[8]=4.680380e+00;
	C_params[9]=0.000000e+00;
	C_params[10] =0.000000e+00;
	C_params[11] =0.000000e+00;
	C_params[12] =0.000000e+00;
	C_params[13] =0.000000e+00;
	C_params[14]=0.000000e+00;
	C_params[15] =0.000000e+00;
	C_params[16] =1.000000e+00;
	C_params[17] =0.000000e+00;
	C_params[18] =0.000000e+00;
	C_params[19]=0.000000e+00;
	C_params[20]=3.942722e+01;
	C_params[21] =1.354505e-03;
	C_params[22]=5.705004e-05;
	C_params[23] =1.808972e-06;
	C_params[24] =5.209732e-08;
	std::vector<double> Cl_params(25);
	Cl_params[0]=2.378069e+02;
	Cl_params[1]=-5.479380e-01;
	Cl_params[2]=-6.132030e-01;
	Cl_params[3]=0.000000e+00;
	Cl_params[4]=1.648655e+02;
	Cl_params[5]=-3.764330e-01;
	Cl_params[6]=-5.474660e-01;
	Cl_params[7]=0.000000e+00;
	Cl_params[8]=7.462998e+01;
	Cl_params[9]=8.390300e-02;
	Cl_params[10]=-1.626640e-01;
	Cl_params[11]=0.000000e+00;
	Cl_params[12]=4.201879e+01;
	Cl_params[13]=-6.476630e-01;
	Cl_params[14]=1.576530e-01;
	Cl_params[15]=0.000000e+00;
	Cl_params[16]=1.000000e+00;
	Cl_params[17]=-3.329910e-01;
	Cl_params[18]=-2.855520e-01;
	Cl_params[19]=0.000000e+00;
	Cl_params[20]=3.849530e+01;
	Cl_params[21]=5.375530e-03;
	Cl_params[22]=2.608683e-04;
	Cl_params[23]=2.970926e-05;
	Cl_params[24]=1.659470e-06;
	std::vector<double> H_params(25);
	H_params[0]=2.773508e+01;
	H_params[1]=0.000000e+00;
	H_params[2]=0.000000e+00;
	H_params[3]=0.000000e+00;
	H_params[4]=1.191792e+01;
	H_params[5]=0.000000e+00;
	H_params[6]=0.000000e+00;
	H_params[7]=0.000000e+00;
	H_params[8]=5.108128e+00;
	H_params[9]=0.000000e+00;
	H_params[10]=0.000000e+00;
	H_params[11]=0.000000e+00;
	H_params[12]=1.364808e+01;
	H_params[13]=0.000000e+00;
	H_params[14]=0.000000e+00;
	H_params[15]=0.000000e+00;
	H_params[16]=1.000000e+00;
	H_params[17]=0.000000e+00;
	H_params[18]=0.000000e+00;
	H_params[19]=0.000000e+00;
	H_params[20]=4.226471e+01;
	H_params[21]=1.294394e-04;
	H_params[22]=6.796368e-06;
	H_params[23]=4.995272e-07;
	H_params[24]=0.000000e+00;
	vector<Vec3> positions(10);
	positions[0] = Vec3(9.000,	9.000,	9.000);
	positions[1] = Vec3(9.666,	9.666,	9.666);
	positions[2] = Vec3(8.334,	8.334,	8.334);
	positions[3] = Vec3(8.434,	8.334,	8.334);
	positions[4] = Vec3(8.534,	8.334,	8.334);
	positions[5] = Vec3(9.000,	9.000,	3.000);
	positions[6] = Vec3(9.666,	9.666,	3.666);
	positions[7] = Vec3(8.334,	8.334,	2.334);
	positions[8] = Vec3(8.434,	8.334,	2.334);
	positions[9] = Vec3(8.534,	8.334,	2.334);

	steady_clock::time_point start1 = steady_clock::now();
	//define Hbond system
	System refsystem;
	refsystem.addParticle(12.0);
	refsystem.addParticle(35.45);
	refsystem.addParticle(1.0);
	refsystem.addParticle(1.0);
	refsystem.addParticle(1.0);
	refsystem.addParticle(12.0);
	refsystem.addParticle(35.45);
	refsystem.addParticle(1.0);
	refsystem.addParticle(1.0);
	refsystem.addParticle(1.0);
	refsystem.setDefaultPeriodicBoxVectors(Vec3(10.0,0,0),Vec3(0,10.0,0),Vec3(0,0,10.0));
	//define Hbond force
	CustomHbondForce* refforceField = new CustomHbondForce("scale*(A*K2*exBr - (Adi)*(f6*C6/(r^6) + f8*C8/(r^8) + f10*C10/(r^10) + f12*C12/(r^12))); A=Aex-Ael-Ain-Adh; Aex=(Aexch1*Aexch2*Aexch1_sph*Aexch2_sph); Aexch1_sph= 1 + aexch_y101*y101 + aexch_y201*y201 + aexch_y22c1*y22c1; Aexch2_sph= 1 + aexch_y102*y102 + aexch_y202*y202 + aexch_y22c2*y22c2; Ael=(Aelec1*Aelec2*Aelec1_sph*Aelec2_sph); Aelec1_sph= 1 + aelec_y101*y101 + aelec_y201*y201 + aelec_y22c1*y22c1; Aelec2_sph= 1 + aelec_y102*y102 + aelec_y202*y202 + aelec_y22c2*y22c2; Ain=(Aind1*Aind2*Aind1_sph*Aind2_sph); Aind1_sph= 1 + aind_y101*y101 + aind_y201*y201 + aind_y22c1*y22c1; Aind2_sph= 1 + aind_y102*y102 + aind_y202*y202 + aind_y22c2*y22c2; Adh=(Adhf1*Adhf2*Adhf1_sph*Adhf2_sph); Adhf1_sph= 1 + adhf_y101*y101 + adhf_y201*y201 + adhf_y22c1*y22c1; Adhf2_sph= 1 + adhf_y102*y102 + adhf_y202*y202 + adhf_y22c2*y22c2; Adi=(Adisp1*Adisp2*Adisp1_sph*Adisp2_sph); Adisp1_sph= 1 + adisp_y101*y101 + adisp_y201*y201 + adisp_y22c1*y22c1; Adisp2_sph= 1 + adisp_y102*y102 + adisp_y202*y202 + adisp_y22c2*y22c2; K2=(Br^2)/3 + Br + 1; f12 = f10 - exX*((1/39916800)*(X^11)*(1 + X/12)); f10 = f8 - exX*((1/362880)*(X^9)*(1 + X/10)); f8 = f6 - exX*((1/5040)*(X^7)*(1 + X/8)); f6 = 1 - exX*(1 + X * (1 + (1/2)*X*(1 + (1/3)*X*(1 + (1/4)*X*(1 + (1/5)*X*(1 + (1/6)*X)))))); exX = exp(-X); X = Br - r * (2*(B^2)*r + 3*B)/(Br^2 + 3*Br + 3); exBr = exp(-Br); Br = B*r; y101 = cos(theta1); y102 = cos(theta2); y201 = 0.5*(3*cos(theta1)^2 - 1); y202 = 0.5*(3*cos(theta2)^2 - 1); y22c1 = sqrt(0.75)*sin(theta1)^2*cos(2*phi1); y22c2 = sqrt(0.75)*sin(theta2)^2*cos(2*phi2); theta1=angle(d2,d1,a1); theta2=angle(a2,a1,d1); phi1=0; phi2=0; r = distance(a1,d1); B=sqrt(Bexp1*Bexp2); C6=sqrt(C61*C62); C8=sqrt(C81*C82); C10=sqrt(C101*C102); C12=sqrt(C121*C122)");
	//Hbond parameters
	refforceField->addGlobalParameter("scale", 0.5);
	refforceField->addPerDonorParameter("Aexch1");
	refforceField->addPerDonorParameter("aexch_y101");
	refforceField->addPerDonorParameter("aexch_y201");
	refforceField->addPerDonorParameter("aexch_y22c1");
	refforceField->addPerDonorParameter("Aelec1");
	refforceField->addPerDonorParameter("aelec_y101");
	refforceField->addPerDonorParameter("aelec_y201");
	refforceField->addPerDonorParameter("aelec_y22c1");
	refforceField->addPerDonorParameter("Aind1");
	refforceField->addPerDonorParameter("aind_y101");
	refforceField->addPerDonorParameter("aind_y201");
	refforceField->addPerDonorParameter("aind_y22c1");
	refforceField->addPerDonorParameter("Adhf1");
	refforceField->addPerDonorParameter("adhf_y101");
	refforceField->addPerDonorParameter("adhf_y201");
	refforceField->addPerDonorParameter("adhf_y22c1");
	refforceField->addPerDonorParameter("Adisp1");
	refforceField->addPerDonorParameter("adisp_y101");
	refforceField->addPerDonorParameter("adisp_y201");
	refforceField->addPerDonorParameter("adisp_y22c1");
	refforceField->addPerDonorParameter("Bexp1");
	refforceField->addPerDonorParameter("C61");
	refforceField->addPerDonorParameter("C81");
	refforceField->addPerDonorParameter("C101");
	refforceField->addPerDonorParameter("C121");
	refforceField->addPerAcceptorParameter("Aexch2");
	refforceField->addPerAcceptorParameter("aexch_y102");
	refforceField->addPerAcceptorParameter("aexch_y202");
	refforceField->addPerAcceptorParameter("aexch_y22c2");
	refforceField->addPerAcceptorParameter("Aelec2");
	refforceField->addPerAcceptorParameter("aelec_y102");
	refforceField->addPerAcceptorParameter("aelec_y202");
	refforceField->addPerAcceptorParameter("aelec_y22c2");
	refforceField->addPerAcceptorParameter("Aind2");
	refforceField->addPerAcceptorParameter("aind_y102");
	refforceField->addPerAcceptorParameter("aind_y202");
	refforceField->addPerAcceptorParameter("aind_y22c2");
	refforceField->addPerAcceptorParameter("Adhf2");
	refforceField->addPerAcceptorParameter("adhf_y102");
	refforceField->addPerAcceptorParameter("adhf_y202");
	refforceField->addPerAcceptorParameter("adhf_y22c2");
	refforceField->addPerAcceptorParameter("Adisp2");
	refforceField->addPerAcceptorParameter("adisp_y102");
	refforceField->addPerAcceptorParameter("adisp_y202");
	refforceField->addPerAcceptorParameter("adisp_y22c2");
	refforceField->addPerAcceptorParameter("Bexp2");
	refforceField->addPerAcceptorParameter("C62");
	refforceField->addPerAcceptorParameter("C82");
	refforceField->addPerAcceptorParameter("C102");
	refforceField->addPerAcceptorParameter("C122");
        refforceField->setNonbondedMethod(CustomHbondForce::CutoffPeriodic);
	refforceField->setCutoffDistance(5.0);

	int d0 = refforceField->addDonor(0, 1, -1, C_params);
	int d1 = refforceField->addDonor(1, 0, -1, Cl_params);
	int d2 = refforceField->addDonor(2, 0, -1, H_params);
	int d3 = refforceField->addDonor(3, 0, -1, H_params);
	int d4 = refforceField->addDonor(4, 0, -1, H_params);
	int d5 = refforceField->addDonor(5, 6, -1, C_params);
	int d6 = refforceField->addDonor(6, 5, -1, Cl_params);
	int d7 = refforceField->addDonor(7, 5, -1, H_params);
	int d8 = refforceField->addDonor(8, 5, -1, H_params);
	int d9 = refforceField->addDonor(9, 5, -1, H_params);
	int a0 = refforceField->addAcceptor(0, 1, -1, C_params);
	int a1 = refforceField->addAcceptor(1, 0, -1, Cl_params);
	int a2 = refforceField->addAcceptor(2, 0, -1, H_params);
	int a3 = refforceField->addAcceptor(3, 0, -1, H_params);
	int a4 = refforceField->addAcceptor(4, 0, -1, H_params);
	int a5 = refforceField->addAcceptor(5, 6, -1, C_params);
	int a6 = refforceField->addAcceptor(6, 5, -1, Cl_params);
	int a7 = refforceField->addAcceptor(7, 5, -1, H_params);
	int a8 = refforceField->addAcceptor(8, 5, -1, H_params);
	int a9 = refforceField->addAcceptor(9, 5, -1, H_params);

	refforceField->addExclusion(d0,a0);
	refforceField->addExclusion(d0,a1);
	refforceField->addExclusion(d0,a2);
	refforceField->addExclusion(d0,a3);
	refforceField->addExclusion(d0,a4);
	refforceField->addExclusion(d1,a0);
	refforceField->addExclusion(d1,a1);
	refforceField->addExclusion(d1,a2);
	refforceField->addExclusion(d1,a3);
	refforceField->addExclusion(d1,a4);
	refforceField->addExclusion(d2,a0);
	refforceField->addExclusion(d2,a1);
	refforceField->addExclusion(d2,a2);
	refforceField->addExclusion(d2,a3);
	refforceField->addExclusion(d2,a4);
	refforceField->addExclusion(d3,a0);
	refforceField->addExclusion(d3,a1);
	refforceField->addExclusion(d3,a2);
	refforceField->addExclusion(d3,a3);
	refforceField->addExclusion(d3,a4);
	refforceField->addExclusion(d4,a0);
	refforceField->addExclusion(d4,a1);
	refforceField->addExclusion(d4,a2);
	refforceField->addExclusion(d4,a3);
	refforceField->addExclusion(d4,a4);
	refforceField->addExclusion(d5,a5);
	refforceField->addExclusion(d5,a6);
	refforceField->addExclusion(d5,a7);
	refforceField->addExclusion(d5,a8);
	refforceField->addExclusion(d5,a9);
	refforceField->addExclusion(d6,a5);
	refforceField->addExclusion(d6,a6);
	refforceField->addExclusion(d6,a7);
	refforceField->addExclusion(d6,a8);
	refforceField->addExclusion(d6,a9);
	refforceField->addExclusion(d7,a5);
	refforceField->addExclusion(d7,a6);
	refforceField->addExclusion(d7,a7);
	refforceField->addExclusion(d7,a8);
	refforceField->addExclusion(d7,a9);
	refforceField->addExclusion(d8,a5);
	refforceField->addExclusion(d8,a6);
	refforceField->addExclusion(d8,a7);
	refforceField->addExclusion(d8,a8);
	refforceField->addExclusion(d8,a9);
	refforceField->addExclusion(d9,a5);
	refforceField->addExclusion(d9,a6);
	refforceField->addExclusion(d9,a7);
	refforceField->addExclusion(d9,a8);
	refforceField->addExclusion(d9,a9);
	refsystem.addForce(refforceField);

	VerletIntegrator refintegrator(0.01);
	Platform& refplatform = Platform::getPlatformByName("Reference");
	Context refcontext(refsystem, refintegrator, refplatform);
	refcontext.setPositions(positions);
	steady_clock::time_point start1e = steady_clock::now();
	State refstate = refcontext.getState(State::Forces | State::Energy);
	const vector<Vec3>& refforces = refstate.getForces();
	steady_clock::time_point end1 = steady_clock::now();

	//define CAN system
	steady_clock::time_point start2 = steady_clock::now();
	System system;
	system.addParticle(12.0);
	system.addParticle(35.45);
	system.addParticle(1.0);
	system.addParticle(1.0);
	system.addParticle(1.0);
	system.addParticle(12.0);
	system.addParticle(35.45);
	system.addParticle(1.0);
	system.addParticle(1.0);
	system.addParticle(1.0);
	system.setDefaultPeriodicBoxVectors(Vec3(10.0,0,0),Vec3(0,10.0,0),Vec3(0,0,10.0));
	//define CAN force
	CustomAnisotropicNonbondedForce* forceField = new CustomAnisotropicNonbondedForce("(A*K2*exBr - (Adi)*(f6*C6/(r^6) + f8*C8/(r^8) + f10*C10/(r^10) + f12*C12/(r^12))); A=Aex-Ael-Ain-Adh; Aex=(Aexch1*Aexch2*Aexch1_sph*Aexch2_sph); Aexch1_sph= 1 + aexch_y101*y101 + aexch_y201*y201 + aexch_y22c1*y22c1; Aexch2_sph= 1 + aexch_y102*y102 + aexch_y202*y202 + aexch_y22c2*y22c2; Ael=(Aelec1*Aelec2*Aelec1_sph*Aelec2_sph); Aelec1_sph= 1 + aelec_y101*y101 + aelec_y201*y201 + aelec_y22c1*y22c1; Aelec2_sph= 1 + aelec_y102*y102 + aelec_y202*y202 + aelec_y22c2*y22c2; Ain=(Aind1*Aind2*Aind1_sph*Aind2_sph); Aind1_sph= 1 + aind_y101*y101 + aind_y201*y201 + aind_y22c1*y22c1; Aind2_sph= 1 + aind_y102*y102 + aind_y202*y202 + aind_y22c2*y22c2; Adh=(Adhf1*Adhf2*Adhf1_sph*Adhf2_sph); Adhf1_sph= 1 + adhf_y101*y101 + adhf_y201*y201 + adhf_y22c1*y22c1; Adhf2_sph= 1 + adhf_y102*y102 + adhf_y202*y202 + adhf_y22c2*y22c2; Adi=(Adisp1*Adisp2*Adisp1_sph*Adisp2_sph); Adisp1_sph= 1 + adisp_y101*y101 + adisp_y201*y201 + adisp_y22c1*y22c1; Adisp2_sph= 1 + adisp_y102*y102 + adisp_y202*y202 + adisp_y22c2*y22c2; K2=(Br^2)/3 + Br + 1; f12 = f10 - exX*((1/39916800)*(X^11)*(1 + X/12)); f10 = f8 - exX*((1/362880)*(X^9)*(1 + X/10)); f8 = f6 - exX*((1/5040)*(X^7)*(1 + X/8)); f6 = 1 - exX*(1 + X * (1 + (1/2)*X*(1 + (1/3)*X*(1 + (1/4)*X*(1 + (1/5)*X*(1 + (1/6)*X)))))); exX = exp(-X); X = Br - r * (2*(B^2)*r + 3*B)/(Br^2 + 3*Br + 3); exBr = exp(-Br); Br = B*r; y101 = cos(theta1); y102 = cos(theta2); y201 = 0.5*(3*cos(theta1)^2 - 1); y202 = 0.5*(3*cos(theta2)^2 - 1); y22c1 = sqrt(0.75)*sin(theta1)^2*cos(2*phi1); y22c2 = sqrt(0.75)*sin(theta2)^2*cos(2*phi2); B=sqrt(Bexp1*Bexp2); C6=sqrt(C61*C62); C8=sqrt(C81*C82); C10=sqrt(C101*C102); C12=sqrt(C121*C122)");
	//CAN Parameters
	forceField->addPerParticleParameter("Aexch");
	forceField->addPerParticleParameter("aexch_y10");
	forceField->addPerParticleParameter("aexch_y20");
	forceField->addPerParticleParameter("aexch_y22c");
	forceField->addPerParticleParameter("Aelec");
	forceField->addPerParticleParameter("aelec_y10");
	forceField->addPerParticleParameter("aelec_y20");
	forceField->addPerParticleParameter("aelec_y22c");
	forceField->addPerParticleParameter("Aind");
	forceField->addPerParticleParameter("aind_y10");
	forceField->addPerParticleParameter("aind_y20");
	forceField->addPerParticleParameter("aind_y22c");
	forceField->addPerParticleParameter("Adhf");
	forceField->addPerParticleParameter("adhf_y10");
	forceField->addPerParticleParameter("adhf_y20");
	forceField->addPerParticleParameter("adhf_y22c");
	forceField->addPerParticleParameter("Adisp");
	forceField->addPerParticleParameter("adisp_y10");
	forceField->addPerParticleParameter("adisp_y20");
	forceField->addPerParticleParameter("adisp_y22c");
	forceField->addPerParticleParameter("Bexp");
	forceField->addPerParticleParameter("C6");
	forceField->addPerParticleParameter("C8");
	forceField->addPerParticleParameter("C10");
	forceField->addPerParticleParameter("C12");
        forceField->setNonbondedMethod(CustomAnisotropicNonbondedForce::CutoffPeriodic);
	forceField->setCutoffDistance(5.0);

	forceField->addParticle(C_params,0,2,-1,1);
	forceField->addParticle(Cl_params,4,-1,-1,0);
	forceField->addParticle(H_params,4,-1,-1,0);
	forceField->addParticle(H_params,4,-1,-1,0);
	forceField->addParticle(H_params,4,-1,-1,0);
	forceField->addParticle(C_params,0,7,-1,6);
	forceField->addParticle(Cl_params,4,-1,-1,5);
	forceField->addParticle(H_params,4,-1,-1,5);
	forceField->addParticle(H_params,4,-1,-1,5);
	forceField->addParticle(H_params,4,-1,-1,5);
	forceField->addExclusion(1,0);
	forceField->addExclusion(2,0);
	forceField->addExclusion(3,0);
	forceField->addExclusion(4,0);
	forceField->addExclusion(2,1);
	forceField->addExclusion(3,1);
	forceField->addExclusion(4,1);
	forceField->addExclusion(3,2);
	forceField->addExclusion(4,2);
	forceField->addExclusion(4,3);
	forceField->addExclusion(6,5);
	forceField->addExclusion(7,5);
	forceField->addExclusion(8,5);
	forceField->addExclusion(9,5);
	forceField->addExclusion(7,6);
	forceField->addExclusion(8,6);
	forceField->addExclusion(9,6);
	forceField->addExclusion(8,7);
	forceField->addExclusion(9,7);
	forceField->addExclusion(9,8);
	system.addForce(forceField);

	VerletIntegrator integrator(0.01);
	Platform& platform = Platform::getPlatformByName("Reference");
	Context context(system, integrator, platform);
	context.setPositions(positions);
	steady_clock::time_point start2e = steady_clock::now();
	State state = context.getState(State::Forces | State::Energy);
	const vector<Vec3>& forces = state.getForces();
       
	steady_clock::time_point end2 = steady_clock::now();
        duration<double> timer1 = duration_cast<duration<double>>(end1-start1);
        duration<double> timer1e = duration_cast<duration<double>>(end1-start1e);
        duration<double> timer2 = duration_cast<duration<double>>(end2-start2);
        duration<double> timer2e = duration_cast<duration<double>>(end2-start2e);
	cout << timer1.count() << " seconds on customhbondforce" << endl;
	cout << timer1e.count() << " seconds on customhbondforce -- exclusive" << endl;
	cout << timer2.count() << " seconds on customanisotropicnonbondedforce" << endl;
	cout << timer2e.count() << " seconds on customanisotropicnonbondedforce -- exclusive" << endl;

	cout << "Hbond		CAN" << endl;
	cout << refforces[0] << '\t' << forces[0] << '\n';
	cout << refforces[1] << '\t' << forces[1] << '\n';
	cout << refforces[2] << '\t' << forces[2] << '\n';
	cout << refforces[3] << '\t' << forces[3] << '\n';
	cout << refforces[4] << '\t' << forces[4] << '\n';
	cout << refforces[5] << '\t' << forces[5] << '\n';
	cout << refforces[6] << '\t' << forces[6] << '\n';
	cout << refforces[7] << '\t' << forces[7] << '\n';
	cout << refforces[8] << '\t' << forces[8] << '\n';
	cout << refforces[9] << '\t' << forces[9] << '\n';
	cout << refstate.getPotentialEnergy() << '\t' << state.getPotentialEnergy() << '\n';

	ASSERT_EQUAL_VEC(refforces[0], forces[0], TOL);
	ASSERT_EQUAL_VEC(refforces[1], forces[1], TOL);
	ASSERT_EQUAL_VEC(refforces[2], forces[2], TOL);
	ASSERT_EQUAL_VEC(refforces[3], forces[3], TOL);
	ASSERT_EQUAL_VEC(refforces[4], forces[4], TOL);
	ASSERT_EQUAL_VEC(refforces[5], forces[5], TOL);
	ASSERT_EQUAL_VEC(refforces[6], forces[6], TOL);
	ASSERT_EQUAL_VEC(refforces[7], forces[7], TOL);
	ASSERT_EQUAL_VEC(refforces[8], forces[8], TOL);
	ASSERT_EQUAL_VEC(refforces[9], forces[9], TOL);
	ASSERT_EQUAL_TOL(refstate.getPotentialEnergy(), state.getPotentialEnergy(), TOL);
}


int main() {
    try {
	registerCustomAnisotropicNonbondedReferenceKernelFactories();
	test_CH3Cl();
    }
    catch(const exception& e) {
        std::cout << "exception: " << e.what() << endl;
        return 1;
    }
    std::cout << "Done" << endl;
    return 0;
}

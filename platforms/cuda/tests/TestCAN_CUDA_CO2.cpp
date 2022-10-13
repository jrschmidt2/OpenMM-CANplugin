#include "CANCUDATest.h"

using namespace CustomAnisotropicNonbondedPlugin;
using namespace OpenMM;
using namespace std;
using namespace std::chrono;

extern "C" OPENMM_EXPORT void registerCustomAnisotropicNonbondedCudaKernelFactories();

const double TOL = 1e-5;
const double ETOL = 1e-6;

void test_CO2() {
	cout << "TEST CO2" << endl;
	//Define all constants
	vector<double> cparams(25);
	cparams[0]=1.414656e+02;
	cparams[1]=0.000000e+00;
	cparams[2]=4.032290e-01;
	cparams[3]=0.000000e+00;
	cparams[4]=1.421757e+02;
	cparams[5]=0.000000e+00;
	cparams[6]=9.823030e-01;
	cparams[7]=0.000000e+00;
	cparams[8]=1.919441e+01;
	cparams[9]=0.000000e+00;
	cparams[10]=-1.000000e+00;
	cparams[11]=0.000000e+00;
	cparams[12]=1.680096e+01;
	cparams[13]=0.000000e+00;
	cparams[14]=-1.000000e+00;
	cparams[15]=0.000000e+00;
	cparams[16]=1.017705e+00;
	cparams[17]=0.000000e+00;
	cparams[18]=-5.366200e-01;
	cparams[19]=0.000000e+00;
	cparams[20]=4.348231e+01;
	cparams[21]=1.262157e-03;
	cparams[22]=6.962318e-05;
	cparams[23]=3.262527e-06;
	cparams[24]=2.019001e-07;
	vector<double> oparams(25);
	oparams[0]=1.868046e+02;
	oparams[1]=2.045300e-02;
	oparams[2]=-8.284400e-02;
	oparams[3]=0.000000e+00;
	oparams[4]=8.714577e+01;
	oparams[5]=-2.134580e-01;
	oparams[6]=-2.281160e-01;
	oparams[7]=0.000000e+00;
	oparams[8]=3.459401e+01;
	oparams[9]=1.000000e+00;
	oparams[10]=4.099520e-01;
	oparams[11]=0.000000e+00;
	oparams[12]=2.856014e+01;
	oparams[13]=-1.000000e+00;
	oparams[14]=-5.779180e-01;
	oparams[15]=0.000000e+00;
	oparams[16]=7.304340e-01;
	oparams[17]=-1.000000e+01;
	oparams[18]=-3.965300e-01;
	oparams[19]=0.000000e+00;
	oparams[20]=4.450987e+01;
	oparams[21]=9.540046e-04;
	oparams[22]=4.693183e-05;
	oparams[23]=3.163703e-06;
	oparams[24]=1.458385e-07;
	vector<Vec3> positions(9);
	positions[0] = Vec3(0.4752,   0.3752,   0.4752);  
	positions[1] = Vec3(0.5418,   0.4418,   0.5418);
	positions[2] = Vec3(0.4086,   0.3086,   0.4086);
	positions[3] = Vec3(0.4752,   0.3752,   0.5752);  
	positions[4] = Vec3(0.5418,   0.4418,   0.6418);
	positions[5] = Vec3(0.4086,   0.3086,   0.5086);
	positions[6] = Vec3(0.4752,   0.3752,   0.6752);  
	positions[7] = Vec3(0.5418,   0.4418,   0.7418);
	positions[8] = Vec3(0.4086,   0.3086,   0.6086);
	steady_clock::time_point start1 = steady_clock::now();
	//define Hbond system
	System refsystem;
	refsystem.addParticle(12.0);
	refsystem.addParticle(16.0);
	refsystem.addParticle(16.0);
	refsystem.addParticle(12.0);
	refsystem.addParticle(16.0);
	refsystem.addParticle(16.0);
	refsystem.addParticle(12.0);
	refsystem.addParticle(16.0);
	refsystem.addParticle(16.0);
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
        refforceField->setNonbondedMethod(CustomHbondForce::NoCutoff);

	int d0 = refforceField->addDonor(0, 1, -1, cparams);
	int d1 = refforceField->addDonor(1, 0, -1, oparams);
	int d2 = refforceField->addDonor(2, 0, -1, oparams);
	int d3 = refforceField->addDonor(3, 4, -1, cparams);
	int d4 = refforceField->addDonor(4, 3, -1, oparams);
	int d5 = refforceField->addDonor(5, 3, -1, oparams);
	int d6 = refforceField->addDonor(6, 7, -1, cparams);
	int d7 = refforceField->addDonor(7, 6, -1, oparams);
	int d8 = refforceField->addDonor(8, 6, -1, oparams);
	int a0 = refforceField->addAcceptor(0, 1, -1, cparams);
	int a1 = refforceField->addAcceptor(1, 0, -1, oparams);
	int a2 = refforceField->addAcceptor(2, 0, -1, oparams);
	int a3 = refforceField->addAcceptor(3, 4, -1, cparams);
	int a4 = refforceField->addAcceptor(4, 3, -1, oparams);
	int a5 = refforceField->addAcceptor(5, 3, -1, oparams);
	int a6 = refforceField->addAcceptor(6, 7, -1, cparams);
	int a7 = refforceField->addAcceptor(7, 6, -1, oparams);
	int a8 = refforceField->addAcceptor(8, 6, -1, oparams);

	refforceField->addExclusion(d0,a0);
	refforceField->addExclusion(d0,a1);
	refforceField->addExclusion(d0,a2);
	refforceField->addExclusion(d1,a0);
	refforceField->addExclusion(d1,a1);
	refforceField->addExclusion(d1,a2);
	refforceField->addExclusion(d2,a0);
	refforceField->addExclusion(d2,a1);
	refforceField->addExclusion(d2,a2);
	refforceField->addExclusion(d3,a3);
	refforceField->addExclusion(d3,a4);
	refforceField->addExclusion(d3,a5);
	refforceField->addExclusion(d4,a3);
	refforceField->addExclusion(d4,a4);
	refforceField->addExclusion(d4,a5);
	refforceField->addExclusion(d5,a3);
	refforceField->addExclusion(d5,a4);
	refforceField->addExclusion(d5,a5);

	refforceField->addExclusion(d6,a6);
	refforceField->addExclusion(d6,a7);
	refforceField->addExclusion(d6,a8);
	refforceField->addExclusion(d7,a6);
	refforceField->addExclusion(d7,a7);
	refforceField->addExclusion(d7,a8);
	refforceField->addExclusion(d8,a6);
	refforceField->addExclusion(d8,a7);
	refforceField->addExclusion(d8,a8);

	refsystem.addForce(refforceField);

	VerletIntegrator refintegrator(0.01);
	Platform& refplatform = Platform::getPlatformByName("CUDA");
	Context refcontext(refsystem, refintegrator, refplatform);
	refcontext.setPositions(positions);
	refintegrator.step(1);
        steady_clock::time_point start1e = steady_clock::now();
	State refstate = refcontext.getState(State::Forces | State::Energy);
	const vector<Vec3>& refforces = refstate.getForces();
	steady_clock::time_point end1 = steady_clock::now();

	//define CAN system
	steady_clock::time_point start2 = steady_clock::now();
	System system;
	system.addParticle(12.0);
	system.addParticle(16.0);
	system.addParticle(16.0);
	system.addParticle(12.0);
	system.addParticle(16.0);
	system.addParticle(16.0);
	system.addParticle(12.0);
	system.addParticle(16.0);
	system.addParticle(16.0);
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
        forceField->setNonbondedMethod(CustomAnisotropicNonbondedForce::NoCutoff);

	forceField->addParticle(cparams,4,-1,-1,1);
	forceField->addParticle(oparams,4,-1,-1,0);
	forceField->addParticle(oparams,4,-1,-1,0);
	forceField->addParticle(cparams,4,-1,-1,4);
	forceField->addParticle(oparams,4,-1,-1,3);
	forceField->addParticle(oparams,4,-1,-1,3);
	forceField->addParticle(cparams,4,-1,-1,7);
	forceField->addParticle(oparams,4,-1,-1,6);
	forceField->addParticle(oparams,4,-1,-1,6);
	forceField->addExclusion(1,0);
	forceField->addExclusion(2,0);
	forceField->addExclusion(1,2);
	forceField->addExclusion(4,3);
	forceField->addExclusion(4,5);
	forceField->addExclusion(5,3);
	forceField->addExclusion(6,7);
	forceField->addExclusion(6,8);
	forceField->addExclusion(7,8);
	system.addForce(forceField);

	VerletIntegrator integrator(0.01);
	Platform& platform = Platform::getPlatformByName("CUDA");
	Context context(system, integrator, platform);
	context.setPositions(positions);
	integrator.step(1);
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
	ASSERT_EQUAL_TOL(refstate.getPotentialEnergy(), state.getPotentialEnergy(), ETOL);
}


int main() {
    try {
	registerCustomAnisotropicNonbondedCudaKernelFactories();
	test_CO2();
    }
    catch(const exception& e) {
        std::cout << "exception: " << e.what() << endl;
        return 1;
    }
    std::cout << "Done" << endl;
    return 0;
}

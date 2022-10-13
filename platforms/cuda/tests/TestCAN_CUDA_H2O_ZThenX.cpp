#include "CANCUDATest.h"

using namespace CustomAnisotropicNonbondedPlugin;
using namespace OpenMM;
using namespace std;
using namespace std::chrono;

extern "C" OPENMM_EXPORT void registerCustomAnisotropicNonbondedCudaKernelFactories();

const double TOL = 1e-8;

void test_H2O() {
	cout << "TEST H2O" << endl;
	//Define all constants
	std::vector<double> O_params(25);
	O_params[0]=1.742412e+02;
	O_params[1]=1.559280e-01;
	O_params[2]=6.904800e-02;
	O_params[3]=-3.225810e-01;
	O_params[4]=1.162132e+02;
	O_params[5]=1.921310e-01;
	O_params[6]=6.307200e-02;
	O_params[7]=-2.817850e-01;
	O_params[8]=0.000000e+00;
	O_params[9]=2.322250e-01;
	O_params[10]=3.958100e-02;
	O_params[11]=1.827080e-01;
	O_params[12]=3.482895e+01;
	O_params[13]=-2.362900e-01;
	O_params[14]=-2.664300e-02;
	O_params[15]=-5.651100e-02;
	O_params[16]=1.000000e+00;
	O_params[17]=5.338000e-02;
	O_params[18]=2.798300e-02;
	O_params[19]=-8.617900e-02;
	O_params[20]=4.186737e+01;
	O_params[21]=1.461945e-03;
	O_params[22]=9.097242e-05;
	O_params[23]=6.100483e-06;
	O_params[24]=3.076928e-07;
	std::vector<double> H_params(25);
	H_params[0]=1.572426e+01;
	H_params[1]=6.772100e-02;
	H_params[2]=-4.515100e-02;
	H_params[3]=0.000000e+00;
	H_params[4]=4.299823e+00;
	H_params[5]=1.408040e-01;
	H_params[6]=2.020500e-02;
	H_params[7]=0.000000e+00;
	H_params[8]=1.056422e+01;
	H_params[9]=2.661780e-01;
	H_params[10]=7.706700e-02;
	H_params[11]=0.000000e+00;
	H_params[12]=3.861468e+00;
	H_params[13]=-1.000000e+00;
	H_params[14]=2.080820e-01;
	H_params[15]=0.000000e+00;
	H_params[16]=1.000000e+00;
	H_params[17]=-1.061230e-01;
	H_params[18]=3.928500e-02;
	H_params[19]=0.000000e+00;
	H_params[20]=4.280649e+01;
	H_params[21]=4.470987e-05;
	H_params[22]=5.156528e-07;
	H_params[23]=1.344533e-07;
	H_params[24]=0.000000e+00;
	vector<Vec3> positions(10);
	positions[0] = Vec3(0.000,	0.000,	 0.000);
	positions[1] = Vec3(1.037,	0.000,	 0.000);
	positions[2] = Vec3(-0.216,	0.925,	 0.000);
	positions[3] = Vec3(-0.760,	0.668,	-1.061);
	positions[4] = Vec3(-0.447,	1.592,	-1.286);
	positions[5] = Vec3(-0.702,	0.817,	-0.283);
	positions[6] = (positions[1] + positions[2])/2.0;
	positions[7] = (positions[2]-positions[0]).cross(positions[1]-positions[0]) + positions[0];
	positions[8] = (positions[4] + positions[5])/2.0;
	positions[9] = (positions[5]-positions[3]).cross(positions[4]-positions[3]) + positions[3];

	//define Hbond system
	steady_clock::time_point start1 = steady_clock::now();
	System refsystem;
	refsystem.addParticle(16.0);
	refsystem.addParticle(1.0);
	refsystem.addParticle(1.0);
	refsystem.addParticle(16.0);
	refsystem.addParticle(1.0);
	refsystem.addParticle(1.0);
	refsystem.addParticle(0.0);
	refsystem.addParticle(0.0);
	refsystem.addParticle(0.0);
	refsystem.addParticle(0.0);
	refsystem.setVirtualSite(6, new TwoParticleAverageSite(1,2,0.5,0.5));
	refsystem.setVirtualSite(7, new OutOfPlaneSite(0,1,2,0.0,0.5,1.0));
	refsystem.setVirtualSite(8, new TwoParticleAverageSite(4,5,0.5,0.5));
	refsystem.setVirtualSite(9, new OutOfPlaneSite(3,4,5,0.0,0.5,1.0));
	//define Hbond force
	CustomHbondForce* refforceField = new CustomHbondForce("scale*(A*K2*exBr - (Adi)*(f6*C6/(r^6) + f8*C8/(r^8) + f10*C10/(r^10) + f12*C12/(r^12))); A=Aex-Ael-Ain-Adh; Aex=(Aexch1*Aexch2*Aexch1_sph*Aexch2_sph); Aexch1_sph= 1 + aexch_y101*y101 + aexch_y201*y201 + aexch_y22c1*y22c1; Aexch2_sph= 1 + aexch_y102*y102 + aexch_y202*y202 + aexch_y22c2*y22c2; Ael=(Aelec1*Aelec2*Aelec1_sph*Aelec2_sph); Aelec1_sph= 1 + aelec_y101*y101 + aelec_y201*y201 + aelec_y22c1*y22c1; Aelec2_sph= 1 + aelec_y102*y102 + aelec_y202*y202 + aelec_y22c2*y22c2; Ain=(Aind1*Aind2*Aind1_sph*Aind2_sph); Aind1_sph= 1 + aind_y101*y101 + aind_y201*y201 + aind_y22c1*y22c1; Aind2_sph= 1 + aind_y102*y102 + aind_y202*y202 + aind_y22c2*y22c2; Adh=(Adhf1*Adhf2*Adhf1_sph*Adhf2_sph); Adhf1_sph= 1 + adhf_y101*y101 + adhf_y201*y201 + adhf_y22c1*y22c1; Adhf2_sph= 1 + adhf_y102*y102 + adhf_y202*y202 + adhf_y22c2*y22c2; Adi=(Adisp1*Adisp2*Adisp1_sph*Adisp2_sph); Adisp1_sph= 1 + adisp_y101*y101 + adisp_y201*y201 + adisp_y22c1*y22c1; Adisp2_sph= 1 + adisp_y102*y102 + adisp_y202*y202 + adisp_y22c2*y22c2; K2=(Br^2)/3 + Br + 1; f12 = f10 - exX*((1/39916800)*(X^11)*(1 + X/12)); f10 = f8 - exX*((1/362880)*(X^9)*(1 + X/10)); f8 = f6 - exX*((1/5040)*(X^7)*(1 + X/8)); f6 = 1 - exX*(1 + X * (1 + (1/2)*X*(1 + (1/3)*X*(1 + (1/4)*X*(1 + (1/5)*X*(1 + (1/6)*X)))))); exX = exp(-X); X = Br - r * (2*(B^2)*r + 3*B)/(Br^2 + 3*Br + 3); exBr = exp(-Br); Br = B*r; y101 = cos(theta1); y102 = cos(theta2); y201 = 0.5*(3*cos(theta1)^2 - 1); y202 = 0.5*(3*cos(theta2)^2 - 1); y22c1 = sqrt(0.75)*sin(theta1)^2*cos(2*phi1); y22c2 = sqrt(0.75)*sin(theta2)^2*cos(2*phi2); theta1=angle(d2,d1,a1); theta2=angle(a2,a1,d1); phi1=dihedral(d3,d1,d2,a1); phi2=dihedral(a3,a1,a2,d1); r = distance(a1,d1); B=sqrt(Bexp1*Bexp2); C6=sqrt(C61*C62); C8=sqrt(C81*C82); C10=sqrt(C101*C102); C12=sqrt(C121*C122)");
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
	refforceField->setCutoffDistance(10.0);

	int d0 = refforceField->addDonor(0, 6, 7, O_params);
	int d1 = refforceField->addDonor(1, 0, 6, H_params);
	int d2 = refforceField->addDonor(2, 0, 6, H_params);
	int d3 = refforceField->addDonor(3, 8, 9, O_params);
	int d4 = refforceField->addDonor(4, 3, 8, H_params);
	int d5 = refforceField->addDonor(5, 3, 8, H_params);
	int a0 = refforceField->addAcceptor(0, 6, 7, O_params);
	int a1 = refforceField->addAcceptor(1, 0, 6, H_params);
	int a2 = refforceField->addAcceptor(2, 0, 6, H_params);
	int a3 = refforceField->addAcceptor(3, 8, 9, O_params);
	int a4 = refforceField->addAcceptor(4, 3, 8, H_params);
	int a5 = refforceField->addAcceptor(5, 3, 8, H_params);

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
	refsystem.addForce(refforceField);
	VerletIntegrator refintegrator(0.01);
	Platform& refplatform = Platform::getPlatformByName("CUDA");
	Context refcontext(refsystem, refintegrator, refplatform);
	refcontext.setPositions(positions);
	steady_clock::time_point start1e = steady_clock::now();
	State refstate = refcontext.getState(State::Forces | State::Energy);
	const vector<Vec3>& refforces = refstate.getForces();
	steady_clock::time_point end1 = steady_clock::now();

	//define CAN system
	steady_clock::time_point start2 = steady_clock::now();
	System system;
	system.addParticle(16.0);
	system.addParticle(1.0);
	system.addParticle(1.0);
	system.addParticle(16.0);
	system.addParticle(1.0);
	system.addParticle(1.0);
	system.addParticle(0.0);
	system.addParticle(0.0);
	system.addParticle(0.0);
	system.addParticle(0.0);
	system.setVirtualSite(6, new TwoParticleAverageSite(1,2,0.5,0.5));
	system.setVirtualSite(7, new OutOfPlaneSite(0,1,2,0.0,0.5,1.0));
	system.setVirtualSite(8, new TwoParticleAverageSite(4,5,0.5,0.5));
	system.setVirtualSite(9, new OutOfPlaneSite(3,4,5,0.0,0.5,1.0));
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
	forceField->setCutoffDistance(10.0);

	std::vector<double> x_params(25,0);
	std::vector<double> z_params(25,0);
	forceField->addParticle(O_params,0,7,-1,6); //0
	forceField->addParticle(H_params,0,6,-1,0); //1
	forceField->addParticle(H_params,0,6,-1,0); //2
	forceField->addParticle(O_params,0,9,-1,8);  //3
	forceField->addParticle(H_params,0,8,-1,3); //4
	forceField->addParticle(H_params,0,8,-1,3); //5
	forceField->addParticle(z_params,5,-1,-1,-1); //6
	forceField->addParticle(x_params,5,-1,-1,-1); //7
	forceField->addParticle(z_params,5,-1,-1,-1); //8
	forceField->addParticle(x_params,5,-1,-1,-1); //9
	forceField->addExclusion(0,1);
	forceField->addExclusion(0,2);
	forceField->addExclusion(1,2);
	forceField->addExclusion(3,4);
	forceField->addExclusion(3,5);
	forceField->addExclusion(4,5);
	//from offsites
	forceField->addExclusion(6,0);
	forceField->addExclusion(6,1);
	forceField->addExclusion(6,2);
	forceField->addExclusion(6,3);
	forceField->addExclusion(6,4);
	forceField->addExclusion(6,5);
	forceField->addExclusion(6,7);
	forceField->addExclusion(6,8);
	forceField->addExclusion(6,9);
	forceField->addExclusion(7,0);
	forceField->addExclusion(7,1);
	forceField->addExclusion(7,2);
	forceField->addExclusion(7,3);
	forceField->addExclusion(7,4);
	forceField->addExclusion(7,5);
	forceField->addExclusion(7,8);
	forceField->addExclusion(7,9);
	forceField->addExclusion(8,0);
	forceField->addExclusion(8,1);
	forceField->addExclusion(8,2);
	forceField->addExclusion(8,3);
	forceField->addExclusion(8,4);
	forceField->addExclusion(8,5);
	forceField->addExclusion(8,9);
	forceField->addExclusion(9,0);
	forceField->addExclusion(9,1);
	forceField->addExclusion(9,2);
	forceField->addExclusion(9,3);
	forceField->addExclusion(9,4);
	forceField->addExclusion(9,5);

	system.addForce(forceField);

	VerletIntegrator integrator(0.01);
	Platform& platform = Platform::getPlatformByName("CUDA");
	Context context(system, integrator, platform);
	context.setPositions(positions);
	steady_clock::time_point start2e = steady_clock::now();
	State state = context.getState(State::Forces | State::Energy | State::Positions);
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
	registerCustomAnisotropicNonbondedCudaKernelFactories();
	test_H2O();
    }
    catch(const exception& e) {
        std::cout << "exception: " << e.what() << endl;
        return 1;
    }
    std::cout << "Done" << endl;
    return 0;
}

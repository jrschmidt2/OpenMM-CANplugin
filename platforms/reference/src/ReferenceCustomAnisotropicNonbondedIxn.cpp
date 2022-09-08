#include <string.h>
#include <sstream>
#include <iostream>

#include "openmm/reference/SimTKOpenMMUtilities.h"
#include "openmm/reference/ReferenceForce.h"
#include "ReferenceCustomAnisotropicNonbondedIxn.h"

using namespace CustomAnisotropicNonbondedPlugin;
using OpenMM::Vec3;
using OpenMM::ReferenceForce;
using std::map;
using std::pair;
using std::string;
using std::stringstream;
using std::set;
using std::vector;

ReferenceCustomAnisotropicNonbondedIxn::ReferenceCustomAnisotropicNonbondedIxn(const Lepton::CompiledExpression& EExpression, 
const Lepton::CompiledExpression& forceExpR, 
const Lepton::CompiledExpression& forceExpTheta1, 
const Lepton::CompiledExpression& forceExpTheta2, 
const Lepton::CompiledExpression& forceExpPhi1, 
const Lepton::CompiledExpression& forceExpPhi2, 
const Lepton::CompiledExpression& forceExpY101, 
const Lepton::CompiledExpression& forceExpY102, 
const Lepton::CompiledExpression& forceExpY11c1, 
const Lepton::CompiledExpression& forceExpY11c2, 
const Lepton::CompiledExpression& forceExpY11s1, 
const Lepton::CompiledExpression& forceExpY11s2, 
const Lepton::CompiledExpression& forceExpY201, 
const Lepton::CompiledExpression& forceExpY202, 
const Lepton::CompiledExpression& forceExpY21c1, 
const Lepton::CompiledExpression& forceExpY21c2, 
const Lepton::CompiledExpression& forceExpY21s1, 
const Lepton::CompiledExpression& forceExpY21s2, 
const Lepton::CompiledExpression& forceExpY22c1, 
const Lepton::CompiledExpression& forceExpY22c2, 
const Lepton::CompiledExpression& forceExpY22s1, 
const Lepton::CompiledExpression& forceExpY22s2, 
const vector<string>& parameterNames) : cutoff(false), useSwitch(false), periodic(false), 
EExp(EExpression), 
fExpR(forceExpR), 
fExpTheta1(forceExpTheta1), 
fExpTheta2(forceExpTheta2), 
fExpPhi1(forceExpPhi1), 
fExpPhi2(forceExpPhi2), 
fExpY101(forceExpY101), 
fExpY102(forceExpY102), 
fExpY11c1(forceExpY11c1), 
fExpY11c2(forceExpY11c2), 
fExpY11s1(forceExpY11s1), 
fExpY11s2(forceExpY11s2), 
fExpY201(forceExpY201), 
fExpY202(forceExpY202), 
fExpY21c1(forceExpY21c1), 
fExpY21c2(forceExpY21c2), 
fExpY21s1(forceExpY21s1), 
fExpY21s2(forceExpY21s2), 
fExpY22c1(forceExpY22c1), 
fExpY22c2(forceExpY22c2), 
fExpY22s1(forceExpY22s1), 
fExpY22s2(forceExpY22s2), 
paramNames(parameterNames) {

 	expressionSet.registerExpression(this->EExp);
	expressionSet.registerExpression(this->fExpR);
	expressionSet.registerExpression(this->fExpTheta1);
	expressionSet.registerExpression(this->fExpTheta2);
	expressionSet.registerExpression(this->fExpPhi1);
	expressionSet.registerExpression(this->fExpPhi2);
	expressionSet.registerExpression(this->fExpY101);
	expressionSet.registerExpression(this->fExpY102);
	expressionSet.registerExpression(this->fExpY11c1);
	expressionSet.registerExpression(this->fExpY11c2);
	expressionSet.registerExpression(this->fExpY11s1);
	expressionSet.registerExpression(this->fExpY11s2);
	expressionSet.registerExpression(this->fExpY201);
	expressionSet.registerExpression(this->fExpY202);
	expressionSet.registerExpression(this->fExpY21c1);
	expressionSet.registerExpression(this->fExpY21c2);
	expressionSet.registerExpression(this->fExpY21s1);
	expressionSet.registerExpression(this->fExpY21s2);
	expressionSet.registerExpression(this->fExpY22c1);
	expressionSet.registerExpression(this->fExpY22c2);
	expressionSet.registerExpression(this->fExpY22s1);
	expressionSet.registerExpression(this->fExpY22s2);

	rIndex = expressionSet.getVariableIndex("r");
	theta1Index = expressionSet.getVariableIndex("theta1");
	theta2Index = expressionSet.getVariableIndex("theta2");
	phi1Index = expressionSet.getVariableIndex("phi1");
	phi2Index = expressionSet.getVariableIndex("phi2");
	y00i = expressionSet.getVariableIndex("y00");
	y10_1i = expressionSet.getVariableIndex("y10_1");
	y10_2i = expressionSet.getVariableIndex("y10_2");
	y11c_1i = expressionSet.getVariableIndex("y11c_1");
	y11c_2i = expressionSet.getVariableIndex("y11c_2");
	y11s_1i = expressionSet.getVariableIndex("y11s_1");
	y11s_2i = expressionSet.getVariableIndex("y11s_2");
	y20_1i = expressionSet.getVariableIndex("y20_1");
	y20_2i = expressionSet.getVariableIndex("y20_2");
	y21c_1i = expressionSet.getVariableIndex("y21c_1");
	y21c_2i = expressionSet.getVariableIndex("y21c_2");
	y21s_1i = expressionSet.getVariableIndex("y21s_1");
	y21s_2i = expressionSet.getVariableIndex("y21s_2");
	y22c_1i = expressionSet.getVariableIndex("y22c_1");
	y22c_2i = expressionSet.getVariableIndex("y22c_2");
	y22s_1i = expressionSet.getVariableIndex("y22s_1");
	y22s_2i = expressionSet.getVariableIndex("y22s_2");
	for (auto& param : paramNames) {
		for (int j = 1; j < 3; j++) {
			stringstream name;
			name << param << j;
			particleParamIndex.push_back(expressionSet.getVariableIndex(name.str()));
		}
	}

}

ReferenceCustomAnisotropicNonbondedIxn::~ReferenceCustomAnisotropicNonbondedIxn() {
}


//VECTOR FUNCTIONS
//
double ReferenceCustomAnisotropicNonbondedIxn::normalizeVec3(Vec3& vectorToNormalize) const {
	double norm = sqrt(vectorToNormalize.dot(vectorToNormalize));
	if (norm > 0.0) {
		vectorToNormalize *= (1.0/norm);
	}
	return norm;
}


double ReferenceCustomAnisotropicNonbondedIxn::computeAzim(const Vec3& vec1, const Vec3& vec2) {
	double angle;

	double dot = vec1.dot(vec2);
	double norm1 = vec1.dot(vec1);
	double norm2 = vec2.dot(vec2);


	double cosine = dot/sqrt(norm1*norm2);

	if (cosine >= 1)
		angle = 0;
	else if (cosine <= -1)
		angle = M_PI;
	else {
		angle = acos(cosine);
	}

	return angle;
}

double ReferenceCustomAnisotropicNonbondedIxn::computePolar(const Vec3& vec1, const Vec3& vec2, const Vec3& vec3) {
	double angle;
	//Compute orth vec
	Vec3 cross1 = vec1.cross(vec2);
	Vec3 cross2 = vec2.cross(vec3);

	double dot = cross1.dot(cross2);
	Vec3 cross = cross1.cross(cross2);

	//Normalize
	double norm1 = cross1.dot(cross1);
	double norm2 = cross2.dot(cross2);

	if (dot != 0.0)
		dot /= sqrt(norm1*norm2);
	if (dot > 1.0)
		dot = 1.0;
	else if (dot < -1.0)
		dot = -1.0;
	//Check for acos singularity
	if (dot > 0.99 || dot < -0.99) {
		angle = asin(sqrt(cross.dot(cross)/(norm1*norm2)));
		if (dot < 0.0)
			angle = M_PI-angle;
	}
	else{
		angle = acos(dot);
	}
	return angle;
}


//SYSTEM SETUP

void ReferenceCustomAnisotropicNonbondedIxn::setUseCutoff(double distance, const OpenMM::NeighborList& neighbors) {
	cutoff = true;
	cutoffDistance = distance;
	neighborList = &neighbors;
}
/*
void ReferenceCustomAnisotropicNonbondedIxn::setInteractionGroups(const vector<pair<set<int>, set<int> > >& groups) {
	interactionGroups = groups
}
*/

void ReferenceCustomAnisotropicNonbondedIxn::setUseSwitchingFunction(double distance) {
	useSwitch = true;
	switchingDistance = distance;
}

void ReferenceCustomAnisotropicNonbondedIxn::setPeriodic(OpenMM::Vec3* vectors) {
	assert(cutoff);
	assert(vectors[0][0] >= 2.0*cutoffDistance);
	assert(vectors[1][1] >= 2.0*cutoffDistance);
	assert(vectors[2][2] >= 2.0*cutoffDistance);
	periodic = true;
	periodicBoxVectors[0] = vectors[0];
	periodicBoxVectors[1] = vectors[1];
	periodicBoxVectors[2] = vectors[2];
}

//K VECTOR CALCULATION
double ** ReferenceCustomAnisotropicNonbondedIxn::accessAxisParameter(vector<Vec3>& atomCoordinates, int pI, int pX,int pY, int pZ, int axisType) const {
	//Initialize & Identify Axis Type
	Vec3 vectorX, vectorY, vectorZ,vectemp;
	if (axisType == 5) {
		vectorZ = Vec3(0,0,1);
		vectorX = Vec3(1,0,0);
	}
	else {
		vectorZ = atomCoordinates[pZ] - atomCoordinates[pI];
		if (axisType == 4) {
			if (fabs(vectorZ[0]) < 0.866) vectorX = Vec3(1.0,0.0,0.0);
			else vectorX = Vec3(0.0,1.0,0.0);
		}
		else {
			vectorY = atomCoordinates[pY] - atomCoordinates[pI];
			if (axisType == 1) {
				vectemp = vectorZ.cross(vectorY);
				vectorZ += vectorY;
				vectorZ *= 0.5;
				vectorX = vectemp;
			}
			else {
				vectorX = atomCoordinates[pX] - atomCoordinates[pI];
				if (axisType == 2) {
					vectorX += vectorY;
					vectorX *= 0.5;
				}
				else if (axisType == 3) {
					vectorZ += vectorX + vectorY;
					vectorZ *= 1/3.0;
					if (fabs(vectorZ[0]) < 0.866) vectorX = Vec3(1.0,0.0,0.0);
					else vectorX = Vec3(0.0,1.0,0.0);
				}
			}
		}	
	}
	vectorY = vectorZ.cross(vectorX);
	//Store to array
	double** kvecpoint = new double*[3];
	for (int i = 0; i < 3; i++) {
		kvecpoint[i] = new double[3];
	}
	for (int i = 0; i < 3; i++) {
		kvecpoint[0][i] = vectorX[i];
		kvecpoint[1][i] = vectorY[i];
		kvecpoint[2][i] = vectorZ[i];
	}
	return kvecpoint;
}

//INTERACTION CALCULATION

void ReferenceCustomAnisotropicNonbondedIxn::calculateOneIxn(int ii, int jj, vector<Vec3>& atomCoordinates, vector<Vec3>& forces, int* axisTypes, int* atomXs, int* atomYs, int* atomZs, int* calc, Vec3& kx1, Vec3& kx2, Vec3& ky1, Vec3& ky2, Vec3& kz1, Vec3& kz2, double* totalEnergy,int YLM) {

	//Define local coordinates
	//ii,i <-> 1; jj,j <-> 2
	//Distance
	double deltaR[5]; //LastDeltaRIndex
	if (periodic) {
		ReferenceForce::getDeltaRPeriodic(atomCoordinates[ii], atomCoordinates[jj], periodicBoxVectors, deltaR);
	}
	else {
		ReferenceForce::getDeltaR(atomCoordinates[ii], atomCoordinates[jj], deltaR);
	}
	
	double r = deltaR[4]; //RIndex
	double r2 = deltaR[3]; //R2Index
	Vec3 rij = Vec3(deltaR[0],deltaR[1],deltaR[2]); //X,Y,Zindices
	if (cutoff && r >= cutoffDistance) {
		return;
	}
	expressionSet.setVariable(rIndex,r);
	//Angles
	Vec3 rd1,rd2;
	double phi1, phi2, theta1, theta2;
	double Y00,Y10_1,Y10_2,Y11c_1,Y11c_2,Y11s_1,Y11s_2,Y20_1,Y20_2,Y21c_1,Y21c_2,Y21s_1,Y21s_2,Y22c_1,Y22c_2,Y22s_1,Y22s_2;
	double Y10_dt1,Y10_dt2,Y11c_dt1,Y11c_dt2,Y11s_dt1,Y11s_dt2,Y20_dt1,Y20_dt2,Y21c_dt1,Y21c_dt2,Y21s_dt1,Y21s_dt2,Y22c_dt1,Y22c_dt2,Y22s_dt1,Y22s_dt2;
	double Y10_dp1,Y10_dp2,Y11c_dp1,Y11c_dp2,Y11s_dp1,Y11s_dp2,Y20_dp1,Y20_dp2,Y21c_dp1,Y21c_dp2,Y21s_dp1,Y21s_dp2,Y22c_dp1,Y22c_dp2,Y22s_dp1,Y22s_dp2;
	double dEdTheta1=0.0,dEdTheta2=0.0,dEdPhi1=0.0,dEdPhi2=0.0,energy=0.0,dEdR=0.0;
	double dEdY101,dEdY102,dEdY11c1,dEdY11s1,dEdY11c2,dEdY11s2,dEdY201,dEdY202,dEdY21c1,dEdY21s1,dEdY21c2,dEdY21s2,dEdY22c1,dEdY22s1,dEdY22c2,dEdY22s2;

	//IF NO YLM -- PASS BOOL
	if (YLM != 0) {
		if (calc[0] != 1) theta1 = computeAzim(-kz1,-rij);
		else theta1 = 0;
		if (calc[1] != 1) theta2 = computeAzim(-kz2,rij);
		else theta2 = 0;
		if (calc[2] != 1) {
			rd1 = kz1-rij;
			phi1 = computePolar(kx1,kz1,rd1);
		}
		else phi1= 0;
		if (calc[3] != 1) {
			rd2 = kz2+rij;
			phi2 = computePolar(kx2,kz2,rd2);
		}
		else phi2 = 0;	
		expressionSet.setVariable(theta1Index,theta1);
		expressionSet.setVariable(theta2Index,theta2);
		expressionSet.setVariable(phi1Index,phi1);
		expressionSet.setVariable(phi2Index,phi2);

		//Evaluate energy
		energy = EExp.evaluate();
		//R-dependent force
		dEdR = fExpR.evaluate();

		//Angle-dependent forces
		if (calc[0] != 1) dEdTheta1 = fExpTheta1.evaluate();
		else dEdTheta1 = 0.0;
		if (calc[1] != 1) dEdTheta2 = fExpTheta2.evaluate();
		else dEdTheta2 = 0.0;
		if (calc[2] != 1) dEdPhi1 = fExpPhi1.evaluate();
		else dEdPhi1 = 0.0;
		if (calc[3] != 1) dEdPhi2 = fExpPhi2.evaluate();
		else dEdPhi2 = 0.0;

	}

	//IF YLM EXIST -- PASS BOOL
	else {
		//Evaluate angle args
		double argth1, argth2, argph1, argph2;
		double sinth1, sinth2, sinph1, sinph2;
		rd1 = kz1-rij;
		Vec3 cross11 = kx1.cross(kz1);
		Vec3 cross12 = kz1.cross(rd1);
		argph1 = cross11.dot(cross12);
		if (argph1 != 0.0) argph1 /= sqrt(cross11.dot(cross11)*cross12.dot(cross12));
		if (argph1 > 1) argph1 = 1;
		else if (argph1 < -1) argph1 = -1;

		rd2 = kz2+rij;
		Vec3 cross21 = kx2.cross(kz2);
		Vec3 cross22 = kz2.cross(rd2);
		argph2 = cross21.dot(cross22);
		if (argph2 != 0.0) argph2 /= sqrt(cross21.dot(cross21)*cross22.dot(cross22));
		if (argph2 > 1) argph2 = 1;
		else if (argph2 < -1) argph2 = -1;

		argth1 = -kz1.dot(-rij)/(r*sqrt(kz1.dot(kz1)));
		if (argth1 >= 1) argth1 = 1;
		else if (argth1 <= -1) argth1 = -1;

		argth2 = -kz2.dot(rij)/(r*sqrt(kz2.dot(kz2)));
		if (argth2 >= 1) argth2 = 1;
		else if (argth2 <= -1) argth2 = -1;

		//Solve for Ylm
		sinth1 = sqrt(1-argth1*argth1);
		sinth2 = sqrt(1-argth2*argth2);
		sinph1 = sqrt(1-argph1*argph1);
		sinph2 = sqrt(1-argph2*argph2);

		Y00 = cste14;

		Y10_1 = cst34*argth1;
		Y11c_1 = cst34*sinth1*argph1;
		Y11s_1 = cst34*sinth1*sinph1;
		Y20_1 = cst516*(3*argth1*argth1-1);
		Y21c_1 = cst154*argth1*sinth1*argph1;
		Y21s_1 = cst154*argth1*sinth1*sinph1;
		Y22c_1 = cst1516*sinth1*sinth1*(argph1*argph1*2-1);
		Y22s_1 = cst1516*sinth1*sinth1*sinph1*2*argph1;

		Y10_2 = cst34*argth2;
		Y11c_2 = cst34*sinth2*argph2;
		Y11s_2 = cst34*sinth2*sinph2;
		Y20_2 = cst516*(3*argth2*argth2-1);
		Y21c_2 = cst154*argth2*sinth2*argph2;
		Y21s_2 = cst154*argth2*sinth2*sinph2;
		Y22c_2 = cst1516*sinth2*sinth2*(argph2*argph2*2-1);
		Y22s_2 = cst1516*sinth2*sinth2*sinph2*2*argph2;

		expressionSet.setVariable(y00i,Y00);
		expressionSet.setVariable(y10_1i,Y10_1);
		expressionSet.setVariable(y11c_1i,Y11c_1);
		expressionSet.setVariable(y11s_1i,Y11s_1);
		expressionSet.setVariable(y20_1i,Y20_1);
		expressionSet.setVariable(y21c_1i,Y21c_1);
		expressionSet.setVariable(y21s_1i,Y21s_1);
		expressionSet.setVariable(y22c_1i,Y22c_1);
		expressionSet.setVariable(y22s_1i,Y22s_1);
		expressionSet.setVariable(y10_2i,Y10_2);
		expressionSet.setVariable(y11c_2i,Y11c_2);
		expressionSet.setVariable(y11s_2i,Y11s_2);
		expressionSet.setVariable(y20_2i,Y20_2);
		expressionSet.setVariable(y21c_2i,Y21c_2);
		expressionSet.setVariable(y21s_2i,Y21s_2);
		expressionSet.setVariable(y22c_2i,Y22c_2);
		expressionSet.setVariable(y22s_2i,Y22s_2);

		//Evaluate energy
		energy = EExp.evaluate();
		//R-dependent force
		dEdR = fExpR.evaluate();
		//Ylm forces
		dEdY101 = fExpY101.evaluate();
		dEdY102 = fExpY102.evaluate();
		dEdY11c1 = fExpY11c1.evaluate();
		dEdY11c2 = fExpY11c2.evaluate();
		dEdY11s1 = fExpY11s1.evaluate();
		dEdY11s2 = fExpY11s2.evaluate();
		dEdY201 = fExpY201.evaluate();
		dEdY202 = fExpY202.evaluate();
		dEdY21c1 = fExpY21c1.evaluate();
		dEdY21c2 = fExpY21c2.evaluate();
		dEdY21s1 = fExpY21s1.evaluate();
		dEdY21s2 = fExpY21s2.evaluate();
		dEdY22c1 = fExpY22c1.evaluate();
		dEdY22c2 = fExpY22c2.evaluate();
		dEdY22s1 = fExpY22s1.evaluate();
		dEdY22s2 = fExpY22s2.evaluate();
		
		//COMPUTE DEDTHETA1 HERE
		Y10_dt1 = -cst34*sinth1;
		Y11c_dt1 = cst34*argth1*argph1;
		Y11s_dt1 = cst34*argth1*sinph1;
		Y20_dt1 = -cst516*6*sinth1*argth1;
		Y21c_dt1 = cst154*argph1*(2*argth1*argth1-1);
		Y21s_dt1 = cst154*sinph1*(2*argth1*argth1-1);
		Y22c_dt1 = cst1516*(2*argth1*sinth1)*(2*argph1*argph1-1);
		Y22s_dt1 = cst1516*(2*argth1*sinth1)*2*argph1*sinph1;

		dEdTheta1 = dEdY101*Y10_dt1 + dEdY11c1*Y11c_dt1 + dEdY201*Y20_dt1 + dEdY21c1*Y21c_dt1 + dEdY22c1*Y22c_dt1 + dEdY11s1*Y11s_dt1 + dEdY21s1*Y21s_dt1 + dEdY22s1*Y22s_dt1;

		//COMPUTE DEDTHETA2 HERE
		Y10_dt2 = -cst34*sinth2;
		Y11c_dt2 = cst34*argth2*argph2;
		Y11s_dt2 = cst34*argth2*sinph2;
		Y20_dt2 = -cst516*6*sinth2*argth2;
		Y21c_dt2 = cst154*argph2*(2*argth2*argth2-1);
		Y21s_dt2 = cst154*sinph2*(2*argth2*argth2-1);
		Y22c_dt2 = cst1516*(2*argth2*sinth2)*(2*argph2*argph2-1);
		Y22s_dt2 = cst1516*(2*argth2*sinth2)*2*argph2*sinph2;

		dEdTheta2 = dEdY102*Y10_dt2 + dEdY11c2*Y11c_dt2 + dEdY202*Y20_dt2 + dEdY21c2*Y21c_dt2 + dEdY22c2*Y22c_dt2 + dEdY11s2*Y11s_dt2 + dEdY21s2*Y21s_dt2 + dEdY22s2*Y22s_dt2;
		
		//COMPUTE DEDPHI1 HERE
		Y11c_dp1 = -cst34*sinth1*sinph1;
		Y11s_dp1 = cst34*sinth1*argph1;
		Y21c_dp1 = -cst154*sinth1*argth1*sinph1;
		Y21s_dp1 = cst154*argth1*sinth1*argph1;
		Y22c_dp1 = -2*cst1516*sinth1*sinth1*2*argph1*sinph1;
		Y22s_dp1 = 2*cst1516*sinth1*sinth1*(2*argph1*argph1-1);

		dEdPhi1 = dEdY11c1*Y11c_dp1 + dEdY21c1*Y21c_dp1 + dEdY22c1*Y22c_dp1+ dEdY11s1*Y11s_dp1 + dEdY21s1*Y21s_dp1 + dEdY22s1*Y22s_dp1;

		//COMPUTE DEDPHI2 HERE
		Y11c_dp2 = -cst34*sinph2*sinth2;
		Y11s_dp2 = cst34*sinth2*argph2;
		Y21c_dp2 = -cst154*sinph2*sinth2*argth2;
		Y21s_dp2 = cst154*argth2*sinth2*argph2;
		Y22c_dp2 = -2*cst1516*sinth2*sinth2*2*argph2*sinph2;
		Y22s_dp2 = 2*cst1516*sinth2*sinth2*(2*argph2*argph2-1);

		dEdPhi2 = dEdY11c2*Y11c_dp2 + dEdY21c2*Y21c_dp2 + dEdY22c2*Y22c_dp2+ dEdY11s2*Y11s_dp2 + dEdY21s2*Y21s_dp2 + dEdY22s2*Y22s_dp2;
	}
	//Switching force?
	double switchValue = 1.0;
	if (useSwitch) {
		if (r > switchingDistance) {
			double t = (r-switchingDistance)/(cutoffDistance-switchingDistance);
			switchValue = 1+t*t*t*(-10*(15-t*6));
			double switchDeriv = t*t*(-30+t*(60-t*30))/(cutoffDistance-switchingDistance);
			dEdR = switchValue*dEdR + energy*switchDeriv/r;
			energy *= switchValue;
		}
	}

	dEdR /= r;
	for (int k = 0; k < 3; k++) {		
		forces[ii][k] += dEdR*rij[k];
		forces[jj][k] -= dEdR*rij[k];
	}
	//Forces from azimuthal angle
	double lengthphicross,termA,termC;
	Vec3 phicross;
	//Theta1
	if (dEdTheta1*dEdTheta1 > 0.0) {
		phicross = (-kz1).cross(-rij);
		lengthphicross = sqrt(phicross.dot(phicross));
		if (lengthphicross < 1.0e-06) lengthphicross= 1.0e-06;
		termA = dEdTheta1/(kz1.dot(kz1)*lengthphicross);
		termC = -dEdTheta1/(r2*lengthphicross);
		vector<Vec3> deltacrossphi1(3);
		deltacrossphi1[0] = (-kz1).cross(phicross);
		deltacrossphi1[2] = (-rij).cross(phicross);
		for (int k = 0; k < 3; k++) {
			deltacrossphi1[0][k] *=termA;
			deltacrossphi1[2][k] *=termC;
		}
		deltacrossphi1[1] = -(deltacrossphi1[0] + deltacrossphi1[2]);
		forces[ii] += deltacrossphi1[1]; 
		forces[jj] += deltacrossphi1[2];
		if (axisTypes[ii] == 4 || axisTypes[ii] == 0 || axisTypes[ii] == 2) 
			forces[atomZs[ii]] += deltacrossphi1[0];
		else if (axisTypes[ii] == 1) {
			forces[atomYs[ii]] += deltacrossphi1[0]/2.0;
			forces[atomZs[ii]] += deltacrossphi1[0]/2.0;
		}
		else if (axisTypes[ii] == 3) {
			forces[atomZs[ii]] += deltacrossphi1[0]/3.0;
			forces[atomXs[ii]] += deltacrossphi1[0]/3.0;
			forces[atomYs[ii]] += deltacrossphi1[0]/3.0;
		}
	}
	//Theta2
	if (dEdTheta2*dEdTheta2 > 0.0) {
		phicross = (-kz2).cross(rij);
		lengthphicross = sqrt(phicross.dot(phicross));
		if (lengthphicross < 1.0e-06) lengthphicross= 1.0e-06;
		termA = dEdTheta2/(kz2.dot(kz2)*lengthphicross);
		termC = -dEdTheta2/(r2*lengthphicross);
		vector<Vec3> deltacrossphi2(3);
		deltacrossphi2[0] = (-kz2).cross(phicross);
		deltacrossphi2[2] = (rij).cross(phicross);
		for (int k = 0; k < 3; k++) {
			deltacrossphi2[0][k] *=termA;
			deltacrossphi2[2][k] *=termC;
		}
		deltacrossphi2[1] = -(deltacrossphi2[0] + deltacrossphi2[2]);
		forces[jj] += deltacrossphi2[1]; 
		forces[ii] += deltacrossphi2[2];
		if (axisTypes[jj] == 4 || axisTypes[jj] == 0 || axisTypes[jj] == 2) 
			forces[atomZs[jj]] += deltacrossphi2[0];
		else if (axisTypes[jj] == 1) {
			forces[atomZs[jj]] += deltacrossphi2[0]/2.0;
			forces[atomYs[jj]] += deltacrossphi2[0]/2.0;
		}	
		else if (axisTypes[jj] == 3) {
			forces[atomZs[jj]] += deltacrossphi2[0]/3.0;
			forces[atomXs[jj]] += deltacrossphi2[0]/3.0;
			forces[atomYs[jj]] += deltacrossphi2[0]/3.0;
		}
	}
	//Forces from polar angle
	vector<Vec3> internalF(4);
	double kzmag2,normCross1,normCross2,normBC;
	double forceFactors[4];
	Vec3 cross1,cross2,s;
	//Phi1
	if (dEdPhi1*dEdPhi1 > 0.0) {
		cross1 = -ky1;
		cross2 = (kz1).cross(rd1);
		kzmag2 = kz1.dot(kz1);
		normCross1 = (cross1).dot(cross1);
		normCross2 = (cross2).dot(cross2);
		if (normCross2 < 1.0e-06) normCross2 = 1.0e-06;
		normBC= sqrt(kzmag2);
		forceFactors[0] = -dEdPhi1*normBC/normCross1;
		forceFactors[1] = (kx1).dot(kz1)/kzmag2;
		forceFactors[2] = rd1.dot(kz1)/kzmag2;
		forceFactors[3] = dEdPhi1*normBC/normCross2;
		for (int k = 0; k < 3; k++) {
			internalF[0][k] = forceFactors[0]*cross1[k];
			internalF[3][k] = forceFactors[3]*cross2[k];
			s[k] = forceFactors[1] * internalF[0][k] - forceFactors[2]*internalF[3][k];
		}
		internalF[1] = internalF[0] -s;
		internalF[2] = internalF[3] +s;
		forces[ii] -= internalF[1];
		forces[jj] += internalF[3];
		if (axisTypes[ii] == 0) {
			forces[atomXs[ii]] += internalF[0];
			forces[atomZs[ii]] -= internalF[2];
		}
		else if (axisTypes[ii] == 1) {
			forces[atomYs[ii]] -= (internalF[2])/2.0;
			forces[atomZs[ii]] -= (internalF[2])/2.0;
			Vec3 f2 = (atomCoordinates[atomZs[ii]]-atomCoordinates[ii]).cross(internalF[0]);
			Vec3 f3 = (internalF[0]).cross(atomCoordinates[atomYs[ii]]-atomCoordinates[ii]) + Vec3(0.5*(internalF[0][0]),0.5*internalF[0][1],0.5*internalF[0][2]);
			forces[atomYs[ii]] += f2;
			forces[atomZs[ii]] += f3;
			forces[ii] += internalF[0]-f2-f3;
		}
		else if (axisTypes[ii] == 2) {
			forces[atomXs[ii]] += internalF[0]/2.0;
			forces[atomYs[ii]] += internalF[0]/2.0;
			forces[atomZs[ii]] -= internalF[2];
		}
		else if (axisTypes[ii] == 3) {
			forces[atomXs[ii]] -= (internalF[2])/3.0;
			forces[atomYs[ii]] -= (internalF[2])/3.0;
			forces[atomZs[ii]] -= (internalF[2])/3.0;
		}
		else if (axisTypes[ii] == 4) {
			forces[atomZs[ii]] -= internalF[2];
		}
			
	}
	//Phi2
	if (dEdPhi2*dEdPhi2 > 0.0) {
		cross1 = -ky2;
		cross2 = (kz2).cross(rd2);
		kzmag2 = kz2.dot(kz2);
		normCross1 = (cross1).dot(cross1);
		normCross2 = (cross2).dot(cross2);
		if (normCross2 < 1.0e-06) normCross2 = 1.0e-06;
		normBC= sqrt(kzmag2);
		forceFactors[0] = -dEdPhi2*normBC/normCross1;
		forceFactors[1] = (kx2).dot(kz2)/kzmag2;
		forceFactors[2] = rd2.dot(kz2)/kzmag2;
		forceFactors[3] = dEdPhi2*normBC/normCross2;
		for (int k = 0; k < 3; k++) {
			internalF[0][k] = forceFactors[0]*cross1[k];
			internalF[3][k] = forceFactors[3]*cross2[k];
			s[k] = forceFactors[1] * internalF[0][k] - forceFactors[2]*internalF[3][k];
		}
		internalF[1] = internalF[0] - s;
		internalF[2] = internalF[3] + s;
		forces[jj] -= internalF[1];
		forces[ii] += internalF[3];
		if (axisTypes[jj] == 0) {
			forces[atomXs[jj]] += internalF[0];
			forces[atomZs[jj]] -= internalF[2];
		}
		else if (axisTypes[jj] == 1) {
			forces[atomYs[jj]] -= (internalF[2])/2.0;
			forces[atomZs[jj]] -= (internalF[2])/2.0;
			Vec3 f2 = (atomCoordinates[atomZs[jj]]-atomCoordinates[jj]).cross(internalF[0]);
			Vec3 f3 = (internalF[0]).cross(atomCoordinates[atomYs[jj]]-atomCoordinates[jj]) + Vec3(0.5*(internalF[0][0]),0.5*internalF[0][1],0.5*internalF[0][2]);
			forces[atomYs[jj]] += f2;
			forces[atomZs[jj]] += f3;
			forces[jj] += internalF[0]-f2-f3;
		}
		else if (axisTypes[jj] == 2) {
			forces[atomXs[jj]] += internalF[0]/2.0;
			forces[atomYs[jj]] += internalF[0]/2.0;
			forces[atomZs[jj]] -= internalF[2];
		}
		else if (axisTypes[jj] == 3) {
			forces[atomXs[jj]] -= (internalF[2])/3.0;
			forces[atomYs[jj]] -= (internalF[2])/3.0;
			forces[atomZs[jj]] -= (internalF[2])/3.0;
		}
		else if (axisTypes[jj] == 4) {
			forces[atomZs[jj]] -= internalF[2];
		}
	}
	if (totalEnergy)
		*totalEnergy += energy;
}

void ReferenceCustomAnisotropicNonbondedIxn::calculatePairIxn(int numberOfAtoms, vector<Vec3>& atomCoordinates, int* axisTypes, int* atomXs, int* atomYs, int* atomZs, int* calc, double** atomParameters, vector<set<int> >& exclusions, double* fixedParameters, const map<string, double>& globalParameters, vector<Vec3>& forces, double*  totalEnergy, int YLM) {
	//Obtain local axis systems
	vector<Vec3> kxvecs,kyvecs,kzvecs;
	double** kvecpoint = new double*[3];
	for (unsigned int i=0; i< numberOfAtoms; i++) {
		kvecpoint = accessAxisParameter(atomCoordinates,i,atomXs[i],atomYs[i],atomZs[i],axisTypes[i]);
		kxvecs.push_back(Vec3(kvecpoint[0][0],kvecpoint[0][1],kvecpoint[0][2]));
		kyvecs.push_back(Vec3(kvecpoint[1][0],kvecpoint[1][1],kvecpoint[1][2]));
		kzvecs.push_back(Vec3(kvecpoint[2][0],kvecpoint[2][1],kvecpoint[2][2]));
	}

	for (auto& param : globalParameters)
		expressionSet.setVariable(expressionSet.getVariableIndex(param.first),param.second);
	//Calculate energy
	if (cutoff) {
		for (auto& pair : *neighborList) {
			for (int j = 0; j < (int) paramNames.size(); j++) {
				expressionSet.setVariable(particleParamIndex[j*2],atomParameters[pair.first][j]);
				expressionSet.setVariable(particleParamIndex[j*2+1],atomParameters[pair.second][j]);
			}
			calculateOneIxn(pair.first,pair.second,atomCoordinates,forces,axisTypes,atomXs,atomYs,atomZs,calc,kxvecs[pair.first],kxvecs[pair.second],kyvecs[pair.first],kyvecs[pair.second],kzvecs[pair.first],kzvecs[pair.second],totalEnergy,YLM);
		}
	}

	else {
		for (int ii = 0; ii < numberOfAtoms; ii++) {
			for (int jj = ii+1; jj < numberOfAtoms; jj++) {
				if (exclusions[jj].find(ii) == exclusions[jj].end()) {
					for (int j = 0; j < (int) paramNames.size(); j++) {
						expressionSet.setVariable(particleParamIndex[j*2],atomParameters[ii][j]);
						expressionSet.setVariable(particleParamIndex[j*2+1],atomParameters[jj][j]);
					}
					calculateOneIxn(ii,jj,atomCoordinates, forces, axisTypes,atomXs,atomYs,atomZs,calc,kxvecs[ii],kxvecs[jj],kyvecs[ii],kyvecs[jj],kzvecs[ii],kzvecs[jj], totalEnergy,YLM);
				}		
			}
		}
	}
}



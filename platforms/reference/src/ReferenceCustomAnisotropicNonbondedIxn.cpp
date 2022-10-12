#ifdef WIN32
#define _USE_MATH_DEFINES
#endif

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

ReferenceCustomAnisotropicNonbondedIxn::ReferenceCustomAnisotropicNonbondedIxn(const Lepton::CompiledExpression& EExpression, const Lepton::CompiledExpression& forceExpR, 
const Lepton::CompiledExpression& forceExpTheta1, const Lepton::CompiledExpression& forceExpTheta2, 
const Lepton::CompiledExpression& forceExpPhi1, const Lepton::CompiledExpression& forceExpPhi2,
const vector<string>& parameterNames) : cutoff(false), useSwitch(false), periodic(false), 
EExp(EExpression), fExpR(forceExpR), fExpTheta1(forceExpTheta1), fExpTheta2(forceExpTheta2), 
fExpPhi1(forceExpPhi1), fExpPhi2(forceExpPhi2), paramNames(parameterNames) {
    expressionSet.registerExpression(this->EExp);
    expressionSet.registerExpression(this->fExpR);
    expressionSet.registerExpression(this->fExpTheta1);
    expressionSet.registerExpression(this->fExpTheta2);
    expressionSet.registerExpression(this->fExpPhi1);
    expressionSet.registerExpression(this->fExpPhi2);

    rIndex = expressionSet.getVariableIndex("r");
    theta1Index = expressionSet.getVariableIndex("theta1");
    theta2Index = expressionSet.getVariableIndex("theta2");
    phi1Index = expressionSet.getVariableIndex("phi1");
    phi2Index = expressionSet.getVariableIndex("phi2");
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
    //compute
    double dot = vec1.dot(vec2);
    double norm1 = vec1.dot(vec1);
    double norm2 = vec2.dot(vec2);
    double cosine = dot/sqrt(norm1*norm2);
    //check limits
    if (cosine >= 1) angle = 0;
    else if (cosine <= -1) angle = M_PI;
    else angle = acos(cosine);

    return angle;
}

double ReferenceCustomAnisotropicNonbondedIxn::computePolar(const Vec3& vec1, const Vec3& vec2, const Vec3& vec3) {
    double angle;
    //compute
    Vec3 cross1 = vec1.cross(vec2);
    Vec3 cross2 = vec2.cross(vec3);
    double dot = cross1.dot(cross2);
    Vec3 cross = cross1.cross(cross2);
    //normalize
    double norm1 = cross1.dot(cross1);
    double norm2 = cross2.dot(cross2);
    //check limits
    if (dot != 0.0) dot /= sqrt(norm1*norm2);
    if (dot > 1.0) dot = 1.0;
    else if (dot < -1.0) dot = -1.0;
    //Check for acos singularity
    if (dot > 0.99 || dot < -0.99) {
        angle = asin(sqrt(cross.dot(cross)/(norm1*norm2)));
        if (dot < 0.0)angle = M_PI-angle;
    }
    else angle = acos(dot);

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

//INTERACTION CALCULATION

void ReferenceCustomAnisotropicNonbondedIxn::calculateOneIxn(int ii, int jj, vector<Vec3>& atomCoordinates, vector<Vec3>& forces, int* axisTypes, int* atomXs, int* atomYs, int* atomZs, int* calc, Vec3& kx1, Vec3& kx2, Vec3& ky1, Vec3& ky2, Vec3& kz1, Vec3& kz2, double* totalEnergy) {
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
    if (cutoff && r >= cutoffDistance) return;
    expressionSet.setVariable(rIndex,r);
    //Angles
    Vec3 rd1,rd2;
    double phi1, phi2, theta1, theta2;
    double dEdTheta1=0.0,dEdTheta2=0.0,dEdPhi1=0.0,dEdPhi2=0.0,energy=0.0,dEdR=0.0;
    if (calc[0] != 1) {
        theta1 = computeAzim(-kz1,-rij);
        expressionSet.setVariable(theta1Index,theta1);
    }
    if (calc[1] != 1) {
        theta2 = computeAzim(-kz2,rij);
        expressionSet.setVariable(theta2Index,theta2);
    }
    if (calc[2] != 1) {
        rd1 = kz1-rij;
        phi1 = computePolar(kx1,kz1,rd1);
        expressionSet.setVariable(phi1Index,phi1);
    }
    if (calc[3] != 1) {
        rd2 = kz2+rij;
        phi2 = computePolar(kx2,kz2,rd2);
        expressionSet.setVariable(phi2Index,phi2);
    }

    //Angle-dependent forces
    if (calc[0] != 1) dEdTheta1 = fExpTheta1.evaluate();
    else dEdTheta1 = 0.0;
    if (calc[1] != 1) dEdTheta2 = fExpTheta2.evaluate();
    else dEdTheta2 = 0.0;
    if (calc[2] != 1) dEdPhi1 = fExpPhi1.evaluate();
    else dEdPhi1 = 0.0;
    if (calc[3] != 1) dEdPhi2 = fExpPhi2.evaluate();
    else dEdPhi2 = 0.0;

    //Evaluate energy
    energy = EExp.evaluate();
    //R-dependent force
    dEdR = fExpR.evaluate();
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
    double lengththcross,termA,termC;
    Vec3 thcross;
    //Theta1
    if (dEdTheta1*dEdTheta1 > 0.0) {
        thcross = (-kz1).cross(-rij);
        lengththcross = sqrt(thcross.dot(thcross));
        if (lengththcross < 1.0e-06) lengththcross= 1.0e-06;
        termA = dEdTheta1/(kz1.dot(kz1)*lengththcross);
        termC = -dEdTheta1/(r2*lengththcross);
        vector<Vec3> deltacrossth1(3);
        deltacrossth1[0] = (-kz1).cross(thcross);
        deltacrossth1[2] = (-rij).cross(thcross);
        for (int k = 0; k < 3; k++) {
            deltacrossth1[0][k] *=termA;
            deltacrossth1[2][k] *=termC;
        }
        deltacrossth1[1] = -(deltacrossth1[0] + deltacrossth1[2]);
        forces[ii] += deltacrossth1[1]; 
        forces[jj] += deltacrossth1[2];
        if (axisTypes[ii] == 4 || axisTypes[ii] == 0 || axisTypes[ii] == 2) 
            forces[atomZs[ii]] += deltacrossth1[0];
        else if (axisTypes[ii] == 1) {
            forces[atomYs[ii]] += deltacrossth1[0]/2.0;
            forces[atomZs[ii]] += deltacrossth1[0]/2.0;
        }
        else if (axisTypes[ii] == 3) {
            forces[atomZs[ii]] += deltacrossth1[0]/3.0;
            forces[atomXs[ii]] += deltacrossth1[0]/3.0;
            forces[atomYs[ii]] += deltacrossth1[0]/3.0;
        }
    }
    //Theta2
    if (dEdTheta2*dEdTheta2 > 0.0) {
        thcross = (-kz2).cross(rij);
        lengththcross = sqrt(thcross.dot(thcross));
        if (lengththcross < 1.0e-06) lengththcross= 1.0e-06;
        termA = dEdTheta2/(kz2.dot(kz2)*lengththcross);
        termC = -dEdTheta2/(r2*lengththcross);
        vector<Vec3> deltacrossth2(3);
        deltacrossth2[0] = (-kz2).cross(thcross);
        deltacrossth2[2] = (rij).cross(thcross);
        for (int k = 0; k < 3; k++) {
            deltacrossth2[0][k] *=termA;
            deltacrossth2[2][k] *=termC;
        }
        deltacrossth2[1] = -(deltacrossth2[0] + deltacrossth2[2]);
        for (int k = 0; k < 3; k++) {
        }
        forces[jj] += deltacrossth2[1]; 
        forces[ii] += deltacrossth2[2];
        if (axisTypes[jj] == 4 || axisTypes[jj] == 0 || axisTypes[jj] == 2) 
            forces[atomZs[jj]] += deltacrossth2[0];
        else if (axisTypes[jj] == 1) {
            forces[atomZs[jj]] += deltacrossth2[0]/2.0;
            forces[atomYs[jj]] += deltacrossth2[0]/2.0;
        }	
        else if (axisTypes[jj] == 3) {
            forces[atomZs[jj]] += deltacrossth2[0]/3.0;
            forces[atomXs[jj]] += deltacrossth2[0]/3.0;
            forces[atomYs[jj]] += deltacrossth2[0]/3.0;
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

void ReferenceCustomAnisotropicNonbondedIxn::calculatePairIxn(int numberOfAtoms, vector<Vec3>& atomCoordinates, int* axisTypes, int* atomXs, int* atomYs, int* atomZs, int* calc, double** atomParameters, vector<set<int> >& exclusions, double* fixedParameters, const map<string, double>& globalParameters, vector<Vec3>& forces, double*  totalEnergy) {
    //Obtain local axis systems
    vector<Vec3> kxvecs,kyvecs,kzvecs;
//    double* kvecpoint = new double[3];
    for (unsigned int i=0; i< numberOfAtoms; i++) {
        int pI=i; int pX=atomXs[i]; int pY=atomYs[i]; int pZ=atomZs[i]; int axisType=axisTypes[i];
        //Initialize & Identify Axis Type
        Vec3 vectorX, vectorY, vectorZ,vectemp;
        if (axisType == 5) {
            vectorZ = Vec3(0,0,1);
            vectorX = Vec3(1,0,0);
        }
        else {
            double deltaIZ[5]; //LastDeltaRIndex
            if (periodic) {
                ReferenceForce::getDeltaRPeriodic(atomCoordinates[pI], atomCoordinates[pZ], periodicBoxVectors, deltaIZ);
            }
            else {
                ReferenceForce::getDeltaR(atomCoordinates[pI], atomCoordinates[pZ], deltaIZ);
            }	
            vectorZ = Vec3(deltaIZ[0],deltaIZ[1],deltaIZ[2]); //X,Y,Zindices
            if (axisType == 4) {
                if (fabs(vectorZ[0]) < 0.866) vectorX = Vec3(1.0,0.0,0.0);
                else vectorX = Vec3(0.0,1.0,0.0);
            }
            else {
                double deltaIY[5]; //LastDeltaRIndex
                if (periodic) {
                    ReferenceForce::getDeltaRPeriodic(atomCoordinates[pI], atomCoordinates[pY], periodicBoxVectors, deltaIY);
                }
                else {
                    ReferenceForce::getDeltaR(atomCoordinates[pI], atomCoordinates[pY], deltaIY);
                }	
                vectorY = Vec3(deltaIY[0],deltaIY[1],deltaIY[2]); //X,Y,Zindices
                if (axisType == 1) {
                    vectemp = vectorZ.cross(vectorY);
                    vectorZ += vectorY;
                    vectorZ *= 0.5;
                    vectorX = vectemp;
                }
                else {
                    double deltaIX[5]; //LastDeltaRIndex
                    if (periodic) {
                        ReferenceForce::getDeltaRPeriodic(atomCoordinates[pI], atomCoordinates[pX], periodicBoxVectors, deltaIX);
                    }
                    else {
                        ReferenceForce::getDeltaR(atomCoordinates[pI], atomCoordinates[pX], deltaIX);
                    }	
                    vectorX = Vec3(deltaIX[0],deltaIX[1],deltaIX[2]); //X,Y,Zindices
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
        kxvecs.push_back(vectorX); kyvecs.push_back(vectorY); kzvecs.push_back(vectorZ);
    }
    for (auto& param : globalParameters) { 
        expressionSet.setVariable(expressionSet.getVariableIndex(param.first),param.second);
    }
    //Calculate energy
    if (cutoff) {
        for (auto& pair : *neighborList) {
            for (int j = 0; j < (int) paramNames.size(); j++) {
                expressionSet.setVariable(particleParamIndex[j*2],atomParameters[pair.first][j]);
                expressionSet.setVariable(particleParamIndex[j*2+1],atomParameters[pair.second][j]);
            }
            calculateOneIxn(pair.first,pair.second,atomCoordinates,forces,axisTypes,atomXs,atomYs,atomZs,calc,kxvecs[pair.first],kxvecs[pair.second],kyvecs[pair.first],kyvecs[pair.second],kzvecs[pair.first],kzvecs[pair.second],totalEnergy);
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
                    calculateOneIxn(ii,jj,atomCoordinates, forces, axisTypes,atomXs,atomYs,atomZs,calc,kxvecs[ii],kxvecs[jj],kyvecs[ii],kyvecs[jj],kzvecs[ii],kzvecs[jj], totalEnergy);
                }		
            }
        }
    }
    kxvecs.clear(); kxvecs.shrink_to_fit();
    kyvecs.clear(); kyvecs.shrink_to_fit();
    kzvecs.clear(); kzvecs.shrink_to_fit();
}



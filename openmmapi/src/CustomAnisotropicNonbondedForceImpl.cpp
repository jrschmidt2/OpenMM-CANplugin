/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2016 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------*
 * Adaptations listed in corresponding header file			      *
 * Check axis types, whether Ylm used					      *
 *									      *
 * Last Edit 07/2018 tjanicki						      *
 * ---------------------------------------------------------------------------*/


#ifdef WIN32
  #define _USE_MATH_DEFINES // Needed to get M_PI
#endif

#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/SplineFitter.h"
#include "openmm/reference/ReferenceTabulatedFunction.h"

#include "internal/CustomAnisotropicNonbondedForceImpl.h"
#include "CustomAnisotropicNonbondedKernels.h"

#include "lepton/ParsedExpression.h"
#include "lepton/Parser.h"

#include <cmath>
#include <sstream>
#include <utility>
#include <algorithm>
#include <iostream>
#include <cstring>

using namespace CustomAnisotropicNonbondedPlugin;
using namespace OpenMM;
using namespace std;
using Lepton::CustomFunction;
using Lepton::ExpressionTreeNode;
using Lepton::Operation;
using Lepton::ParsedExpression;

CustomAnisotropicNonbondedForceImpl::CustomAnisotropicNonbondedForceImpl(const CustomAnisotropicNonbondedForce& owner) : owner(owner) {
}

CustomAnisotropicNonbondedForceImpl::~CustomAnisotropicNonbondedForceImpl() {
}

void CustomAnisotropicNonbondedForceImpl::initialize(ContextImpl& context) {
    kernel = context.getPlatform().createKernel(CalcCustomAnisotropicNonbondedForceKernel::Name(), context);

    // Check for errors in the specification of parameters and exclusions.
    int numParticles = owner.getNumParticles();
    const OpenMM::System& system = context.getSystem();
    if (owner.getNumParticles() != system.getNumParticles())
        throw OpenMMException("CustomAnisotropicNonbondedForce must have exactly as many particles as the System it belongs to.");
    if (owner.getUseSwitchingFunction()) {
        if (owner.getSwitchingDistance() < 0 || owner.getSwitchingDistance() >= owner.getCutoffDistance())
            throw OpenMMException("CustomAnisotropicNonbondedForce: Switching distance must satisfy 0 <= r_switch < r_cutoff");
    }
    vector<set<int> > exclusions(numParticles);
    vector<double> parameters;
    int numParameters = owner.getNumPerParticleParameters();
    for (int i = 0; i < numParticles; i++) {
    	int axisType, atomZ, atomX, atomY;
        owner.getParticleParameters(i, parameters,axisType,atomX,atomY,atomZ);
        if (parameters.size() != numParameters) {
            stringstream msg;
            msg << "CustomAnisotropicNonbondedForce: Wrong number of parameters for particle ";
            msg << i;
            throw OpenMMException(msg.str());
	}
        if (axisType != CustomAnisotropicNonbondedForce::ZThenX && axisType != CustomAnisotropicNonbondedForce::Bisector&&
		axisType != CustomAnisotropicNonbondedForce::ZBisect && axisType != CustomAnisotropicNonbondedForce::ThreeFold &&
		axisType != CustomAnisotropicNonbondedForce::ZOnly && axisType != CustomAnisotropicNonbondedForce::NoAxisType) {
			std::stringstream msg;
			msg << "CustomAnisotropicNonbondedForce: axis type=" << axisType;
			msg << " not currently handled - only AxisTypes[ ";
			msg << CustomAnisotropicNonbondedForce::ZThenX << "," << CustomAnisotropicNonbondedForce::Bisector << ",";
			msg << CustomAnisotropicNonbondedForce::ZBisect << "," << CustomAnisotropicNonbondedForce::ThreeFold << ",";
			msg << CustomAnisotropicNonbondedForce::NoAxisType;
			msg << "] (ZthenX, Bisector, Z-Bisect, ThreeFold, ZOnly, NoAxisType) currently handled.";
			throw OpenMMException(msg.str());
	}
	if (axisType != CustomAnisotropicNonbondedForce::NoAxisType && (atomZ < 0 || atomZ >= numParticles)) {
		std::stringstream msg;
		msg << "CustomAnisotropicNonbondedForce: invalid z axis particle: " << atomZ;
		throw OpenMMException(msg.str());
	}
	if (axisType != CustomAnisotropicNonbondedForce::NoAxisType && axisType != CustomAnisotropicNonbondedForce::ZOnly && axisType != CustomAnisotropicNonbondedForce::Bisector && (atomX < 0 || atomX >= numParticles)) {
		std::stringstream msg;
		msg << "CustomAnisotropicNonbondedForce: invalid x axis particle: " << atomX;
		throw OpenMMException(msg.str());
	}
	if ((axisType == CustomAnisotropicNonbondedForce::ZBisect || axisType == CustomAnisotropicNonbondedForce::ThreeFold || axisType == CustomAnisotropicNonbondedForce::Bisector) && (atomY < 0 || atomY >= numParticles)) {
		std::stringstream msg;
		msg << "CustomAnisotropicNonbondedForce: invalid y axis particle; " << atomY;
		throw OpenMMException(msg.str());
	}
    }
    for (int i = 0; i < owner.getNumExclusions(); i++) {
        int particle1, particle2;
        owner.getExclusionParticles(i, particle1, particle2);
        if (particle1 < 0 || particle1 >= owner.getNumParticles()) {
            stringstream msg;
            msg << "CustomAnisotropicNonbondedForce: Illegal particle index for an exclusion: ";
            msg << particle1;
            throw OpenMMException(msg.str());
        }
        if (particle2 < 0 || particle2 >= owner.getNumParticles()) {
            stringstream msg;
            msg << "CustomAnisotropicNonbondedForce: Illegal particle index for an exclusion: ";
            msg << particle2;
            throw OpenMMException(msg.str());
        }
        if (exclusions[particle1].count(particle2) > 0 || exclusions[particle2].count(particle1) > 0) {
            stringstream msg;
            msg << "CustomAnisotropicNonbondedForce: Multiple exclusions are specified for particles ";
            msg << particle1;
            msg << " and ";
            msg << particle2;
            throw OpenMMException(msg.str());
        }
        exclusions[particle1].insert(particle2);
        exclusions[particle2].insert(particle1);
    }
    string inputExp = owner.getEnergyFunction();
	int atest = 0; 
	int Ytest = 0; 
	if (inputExp.find("theta1") != -1) atest +=1;
	if (inputExp.find("theta2") != -1) atest +=1;
 	if (inputExp.find("phi1") != -1) atest +=1;
	if (inputExp.find("phi2") != -1) atest +=1;
	if (inputExp.find("y00") != -1) Ytest +=1;
	if (inputExp.find("y10_1") != -1) Ytest +=1;
	if (inputExp.find("y10_2") != -1) Ytest +=1;
	if (inputExp.find("y11c_1") != -1) Ytest +=1;
	if (inputExp.find("y11c_2") != -1) Ytest +=1;
	if (inputExp.find("y11s_1") != -1) Ytest +=1;
	if (inputExp.find("y11s_2") != -1) Ytest +=1;
	if (inputExp.find("y20_1") != -1) Ytest +=1;
	if (inputExp.find("y20_2") != -1) Ytest +=1;
	if (inputExp.find("y21c_1") != -1) Ytest +=1;
	if (inputExp.find("y21c_2") != -1) Ytest +=1;
	if (inputExp.find("y21s_1") != -1) Ytest +=1;
	if (inputExp.find("y21s_2") != -1) Ytest +=1;
	if (inputExp.find("y22c_1") != -1) Ytest +=1;
	if (inputExp.find("y22c_2") != -1) Ytest +=1;
	if (inputExp.find("y22s_1") != -1) Ytest +=1;
	if (inputExp.find("y22s_2") != -1) Ytest +=1;
	
    if (atest > 0 && Ytest > 0) throw OpenMMException("CustomAnisotropicNonbondedForce: Mixed methods! Define force in terms of Ylm OR thetas and phis.");
    else if (atest == 0 && Ytest == 0) cout << "CustomAnisotropicNonbondedForce: Warning! No angular components identified. Switch to CustomNonbondedForce for improved efficiency.'\n'";
    if (Ytest > 1) cout << "CustomAnisotropicNonbondedForce: Warning! Multiple spherical harmonics identified.  Ensure that no two spherical harmonics are combined in multiplication or division.'\n'";
   
    if (owner.getNonbondedMethod() == CustomAnisotropicNonbondedForce::CutoffPeriodic) {
        Vec3 boxVectors[3];
        system.getDefaultPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        double cutoff = owner.getCutoffDistance();
        if (cutoff > 0.5*boxVectors[0][0] || cutoff > 0.5*boxVectors[1][1] || cutoff > 0.5*boxVectors[2][2])
            throw OpenMMException("CustomAnisotropicNonbondedForce: The cutoff distance cannot be greater than half the periodic box size.");
    }
    for (int group = 0; group < owner.getNumInteractionGroups(); group++) {
    	set<int> set1, set2;
        owner.getInteractionGroupParameters(group, set1, set2);
        for (set<int>::iterator it = set1.begin(); it != set1.end(); ++it)
            if ((*it < 0) || (*it >= owner.getNumParticles())) {
                stringstream msg;
                msg << "CustomAnisotropicNonbondedForce: Interaction group " << group << " set1 contains a particle index (" << *it << ") "
                    << "not present in system (" << owner.getNumParticles() << " particles).";
                throw OpenMMException(msg.str());
            }
        for (set<int>::iterator it = set2.begin(); it != set2.end(); ++it)
            if ((*it < 0) || (*it >= owner.getNumParticles())) {
                stringstream msg;
                msg << "CustomAnisotropicNonbondedForce: Interaction group " << group << " set2 contains a particle index (" << *it << ") "
                    << "not present in system (" << owner.getNumParticles() << " particles).";
                throw OpenMMException(msg.str());
            }
    }
    kernel.getAs<CalcCustomAnisotropicNonbondedForceKernel>().initialize(context.getSystem(), owner);
}

double CustomAnisotropicNonbondedForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups&(1<<owner.getForceGroup())) != 0)
        return kernel.getAs<CalcCustomAnisotropicNonbondedForceKernel>().execute(context, includeForces, includeEnergy);
    return 0.0;
}

vector<string> CustomAnisotropicNonbondedForceImpl::getKernelNames() {
    vector<string> names;
    names.push_back(CalcCustomAnisotropicNonbondedForceKernel::Name());
    return names;
}

map<string, double> CustomAnisotropicNonbondedForceImpl::getDefaultParameters() {
    map<string, double> parameters;
    for (int i = 0; i < owner.getNumGlobalParameters(); i++)
        parameters[owner.getGlobalParameterName(i)] = owner.getGlobalParameterDefaultValue(i);
    return parameters;
}

void CustomAnisotropicNonbondedForceImpl::updateParametersInContext(ContextImpl& context) {
    kernel.getAs<CalcCustomAnisotropicNonbondedForceKernel>().copyParametersToContext(context, owner);
    context.systemChanged();
}



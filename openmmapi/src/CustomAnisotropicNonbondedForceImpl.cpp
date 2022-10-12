/* -------------------------------------------------------------------------- *
 *                                   OpenMM  CAN                                 *
 * -------------------------------------------------------------------------- */
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
        if (axisType != CustomAnisotropicNonbondedForce::ZThenX && axisType != CustomAnisotropicNonbondedForce::Bisector &&
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
        if (axisType != CustomAnisotropicNonbondedForce::NoAxisType && axisType != CustomAnisotropicNonbondedForce::ZOnly && 
            axisType != CustomAnisotropicNonbondedForce::Bisector && (atomX < 0 || atomX >= numParticles)) {
                std::stringstream msg;
                msg << "CustomAnisotropicNonbondedForce: invalid x axis particle: " << atomX;
                throw OpenMMException(msg.str());
        }
        if ((axisType == CustomAnisotropicNonbondedForce::ZBisect || axisType == CustomAnisotropicNonbondedForce::ThreeFold || 
            axisType == CustomAnisotropicNonbondedForce::Bisector) && (atomY < 0 || atomY >= numParticles)) {
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

    int atest = owner.checkAngles();
    if (atest == 0) {
        cout << "CustomAnisotropicNonbondedForce: Warning! No angular components identified. Switch to CustomNonbondedForce for improved efficiency.'\n'";
    }

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


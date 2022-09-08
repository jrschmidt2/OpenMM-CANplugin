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

/*----------------------------------------------------------------------------*
 *Adaptations specified in corresponding header file			      *
 *									      *
 *
 *Last Edit 07/2018 tjanicki
 *----------------------------------------------------------------------------*/

#include "openmm/Force.h"
#include "openmm/OpenMMException.h"
#include "CustomAnisotropicNonbondedForce.h"
#include "openmm/internal/AssertionUtilities.h"
#include "internal/CustomAnisotropicNonbondedForceImpl.h"
#include <cmath>
#include <map>
#include <sstream>
#include <utility>

using namespace CustomAnisotropicNonbondedPlugin;
using namespace OpenMM;

using std::map;
using std::pair;
using std::set;
using std::string;
using std::stringstream;
using std::vector;

CustomAnisotropicNonbondedForce::CustomAnisotropicNonbondedForce(const string& energy) : energyExpression(energy), nonbondedMethod(NoCutoff), cutoffDistance(1.0),
    switchingDistance(-1.0), useSwitchingFunction(false) {
}

CustomAnisotropicNonbondedForce::CustomAnisotropicNonbondedForce(const CustomAnisotropicNonbondedForce& rhs) {
    // Copy everything and deep copy the tabulated functions
    energyExpression = rhs.energyExpression;
    nonbondedMethod = rhs.nonbondedMethod;
    cutoffDistance = rhs.cutoffDistance;
    switchingDistance = rhs.switchingDistance;
    useSwitchingFunction = rhs.useSwitchingFunction;
    parameters = rhs.parameters;
    globalParameters = rhs.globalParameters;
    energyParameterDerivatives = rhs.energyParameterDerivatives;
    particles = rhs.particles;
    exclusions = rhs.exclusions;
    interactionGroups = rhs.interactionGroups;
    for (vector<FunctionInfo>::const_iterator it = rhs.functions.begin(); it != rhs.functions.end(); it++)
        functions.push_back(FunctionInfo(it->name, it->function->Copy()));
}

CustomAnisotropicNonbondedForce::~CustomAnisotropicNonbondedForce() {
    for (int i = 0; i < (int) functions.size(); i++)
	delete functions[i].function;
}

const string& CustomAnisotropicNonbondedForce::getEnergyFunction() const {
    return energyExpression;
}

void CustomAnisotropicNonbondedForce::setEnergyFunction(const std::string& energy) {
    energyExpression = energy;
}

CustomAnisotropicNonbondedForce::NonbondedMethod CustomAnisotropicNonbondedForce::getNonbondedMethod() const {
    return nonbondedMethod;
}

void CustomAnisotropicNonbondedForce::setNonbondedMethod(NonbondedMethod method) {
    if (method < 0 || method > 2)
        throw OpenMMException("CustomAnisotropicNonbondedForce: Illegal value for nonbonded method");
    nonbondedMethod = method;
}

double CustomAnisotropicNonbondedForce::getCutoffDistance() const {
    return cutoffDistance;
}

void CustomAnisotropicNonbondedForce::setCutoffDistance(double distance) {
    cutoffDistance = distance;
}

bool CustomAnisotropicNonbondedForce::getUseSwitchingFunction() const {
    return useSwitchingFunction;
}

void CustomAnisotropicNonbondedForce::setUseSwitchingFunction(bool use) {
    useSwitchingFunction = use;
}

double CustomAnisotropicNonbondedForce::getSwitchingDistance() const {
    return switchingDistance;
}

void CustomAnisotropicNonbondedForce::setSwitchingDistance(double distance) {
    switchingDistance = distance;
}


int CustomAnisotropicNonbondedForce::addPerParticleParameter(const string& name) {
    parameters.push_back(PerParticleParameterInfo(name));
    return parameters.size()-1;
}

const string& CustomAnisotropicNonbondedForce::getPerParticleParameterName(int index) const {

    ASSERT_VALID_INDEX(index, parameters);
    return parameters[index].name;
}

void CustomAnisotropicNonbondedForce::setPerParticleParameterName(int index, const string& name) {
    ASSERT_VALID_INDEX(index, parameters);
    parameters[index].name = name;
}

int CustomAnisotropicNonbondedForce::addGlobalParameter(const string& name, double defaultValue) {
    globalParameters.push_back(GlobalParameterInfo(name, defaultValue));
    return globalParameters.size()-1;
}

const string& CustomAnisotropicNonbondedForce::getGlobalParameterName(int index) const {
    ASSERT_VALID_INDEX(index, globalParameters);
    return globalParameters[index].name;
}

void CustomAnisotropicNonbondedForce::setGlobalParameterName(int index, const string& name) {
    ASSERT_VALID_INDEX(index, globalParameters);
    globalParameters[index].name = name;
}

double CustomAnisotropicNonbondedForce::getGlobalParameterDefaultValue(int index) const {
    ASSERT_VALID_INDEX(index, globalParameters);
    return globalParameters[index].defaultValue;
}

void CustomAnisotropicNonbondedForce::setGlobalParameterDefaultValue(int index, double defaultValue) {
    ASSERT_VALID_INDEX(index, globalParameters);
    globalParameters[index].defaultValue = defaultValue;
}

void CustomAnisotropicNonbondedForce::addEnergyParameterDerivative(const string& name) {
    for (int i = 0; i < globalParameters.size(); i++)
        if (name == globalParameters[i].name) {
            energyParameterDerivatives.push_back(i);
            return;
        }
    throw OpenMMException(string("addEnergyParameterDerivative: Unknown global parameter '"+name+"'"));
}

const string& CustomAnisotropicNonbondedForce::getEnergyParameterDerivativeName(int index) const {
    ASSERT_VALID_INDEX(index, energyParameterDerivatives);
    return globalParameters[energyParameterDerivatives[index]].name;
}

int CustomAnisotropicNonbondedForce::addParticle(const vector<double>& parameters, int axisType, int atomX, int atomY, int atomZ) {
    particles.push_back(ParticleInfo(parameters, axisType, atomX, atomY, atomZ));
    return particles.size()-1;
}

void CustomAnisotropicNonbondedForce::getParticleParameters(int index, vector<double>& parameters, int& axisType, int& atomX, int& atomY, int& atomZ) const {
    ASSERT_VALID_INDEX(index, particles);
    parameters = particles[index].parameters;
    axisType = particles[index].axisType;
    atomX = particles[index].atomX;
    atomY = particles[index].atomY;
    atomZ = particles[index].atomZ;
}

void CustomAnisotropicNonbondedForce::setParticleParameters(int index, const vector<double>& parameters, int axisType, int atomX, int atomY, int atomZ) {
    ASSERT_VALID_INDEX(index, particles);
    particles[index].parameters = parameters;
    particles[index].axisType = axisType;
    particles[index].atomX = atomX;
    particles[index].atomY = atomY;
    particles[index].atomZ = atomZ;
}

int CustomAnisotropicNonbondedForce::addExclusion(int particle1, int particle2) {
    exclusions.push_back(ExclusionInfo(particle1, particle2));
    return exclusions.size()-1;
}
void CustomAnisotropicNonbondedForce::getExclusionParticles(int index, int& particle1, int& particle2) const {
    ASSERT_VALID_INDEX(index, exclusions);
    particle1 = exclusions[index].particle1;
    particle2 = exclusions[index].particle2;
}

void CustomAnisotropicNonbondedForce::setExclusionParticles(int index, int particle1, int particle2) {
    ASSERT_VALID_INDEX(index, exclusions);
    exclusions[index].particle1 = particle1;
    exclusions[index].particle2 = particle2;
}

void CustomAnisotropicNonbondedForce::createExclusionsFromBonds(const vector<pair<int, int> >& bonds, int bondCutoff) {
    if (bondCutoff < 1)
        return;
    for (auto& bond : bonds)
        if (bond.first < 0 || bond.second < 0 || bond.first >= particles.size() || bond.second >= particles.size())
            throw OpenMMException("createExclusionsFromBonds: Illegal particle index in list of bonds");
    vector<set<int> > exclusions(particles.size());
    vector<set<int> > bonded12(exclusions.size());
    for (auto& bond : bonds) {
        int p1 = bond.first;
        int p2 = bond.second;
        exclusions[p1].insert(p2);
        exclusions[p2].insert(p1);
        bonded12[p1].insert(p2);
        bonded12[p2].insert(p1);
    }
    for (int level = 0; level < bondCutoff-1; level++) {
        vector<set<int> > currentExclusions = exclusions;
        for (int i = 0; i < (int) particles.size(); i++)
            for (int j : currentExclusions[i])
                exclusions[j].insert(bonded12[i].begin(), bonded12[i].end());
    }
    for (int i = 0; i < (int) exclusions.size(); ++i)
        for (int j : exclusions[i])
            if (j < i) 
                addExclusion(j, i);
}

int CustomAnisotropicNonbondedForce::addTabulatedFunction(const std::string& name, TabulatedFunction* function) {
    functions.push_back(FunctionInfo(name, function));
    return functions.size()-1;
}

const TabulatedFunction& CustomAnisotropicNonbondedForce::getTabulatedFunction(int index) const {
    ASSERT_VALID_INDEX(index, functions);
    return *functions[index].function;
}

TabulatedFunction& CustomAnisotropicNonbondedForce::getTabulatedFunction(int index) {
    ASSERT_VALID_INDEX(index, functions);
    return *functions[index].function;
}

const string& CustomAnisotropicNonbondedForce::getTabulatedFunctionName(int index) const {
    ASSERT_VALID_INDEX(index, functions);
    return functions[index].name;
}

int CustomAnisotropicNonbondedForce::addFunction(const std::string& name, const std::vector<double>& values, double min, double max) {
    functions.push_back(FunctionInfo(name, new Continuous1DFunction(values, min, max)));
    return functions.size()-1;
}

void CustomAnisotropicNonbondedForce::getFunctionParameters(int index, std::string& name, std::vector<double>& values, double& min, double& max) const {
    ASSERT_VALID_INDEX(index, functions);
    Continuous1DFunction* function = dynamic_cast<Continuous1DFunction*>(functions[index].function);
    if (function == NULL)
        throw OpenMMException("CustomAnisotropicNonbondedForce: function is not a Continuous1DFunction");
    name = functions[index].name;
    function->getFunctionParameters(values, min, max);
}

void CustomAnisotropicNonbondedForce::setFunctionParameters(int index, const std::string& name, const std::vector<double>& values, double min, double max) {
    ASSERT_VALID_INDEX(index, functions);
    Continuous1DFunction* function = dynamic_cast<Continuous1DFunction*>(functions[index].function);
    if (function == NULL)
        throw OpenMMException("CustomAnisotropicNonbondedForce: function is not a Continuous1DFunction");
    functions[index].name = name;
    function->setFunctionParameters(values, min, max);
}
int CustomAnisotropicNonbondedForce::addInteractionGroup(const std::set<int>& set1, const std::set<int>& set2) {
    for (set<int>::iterator it = set1.begin(); it != set1.end(); ++it)
        ASSERT(*it >= 0);
    for (set<int>::iterator it = set2.begin(); it != set2.end(); ++it)
        ASSERT(*it >= 0);
    interactionGroups.push_back(InteractionGroupInfo(set1, set2));
    return interactionGroups.size()-1;
}

void CustomAnisotropicNonbondedForce::getInteractionGroupParameters(int index, std::set<int>& set1, std::set<int>& set2) const {
    ASSERT_VALID_INDEX(index, interactionGroups);
    set1 = interactionGroups[index].set1;
    set2 = interactionGroups[index].set2;
}

void CustomAnisotropicNonbondedForce::setInteractionGroupParameters(int index, const std::set<int>& set1, const std::set<int>& set2) {
    ASSERT_VALID_INDEX(index, interactionGroups);
    for (set<int>::iterator it = set1.begin(); it != set1.end(); ++it)
        ASSERT_VALID_INDEX(*it, particles);
    for (set<int>::iterator it = set2.begin(); it != set2.end(); ++it)
        ASSERT_VALID_INDEX(*it, particles);
    interactionGroups[index].set1 = set1;
    interactionGroups[index].set2 = set2;
}

ForceImpl* CustomAnisotropicNonbondedForce::createImpl() const {
    return new CustomAnisotropicNonbondedForceImpl(*this);
}

void CustomAnisotropicNonbondedForce::updateParametersInContext(Context& context) {
    dynamic_cast<CustomAnisotropicNonbondedForceImpl&>(getImplInContext(context)).updateParametersInContext(getContextImpl(context));
}

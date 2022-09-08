#ifndef OPENMM_CUSTOMANISOTROPICNONBONDEDFORCE_H_
#define OPENMM_CUSTOMANISOTROPICNONBONDEDFORCE_H_

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

/*------------------------------------------------------------------------------*
 * The following adaptations were made to suit CustomAnisotropicNonbondedPlugin *
 *										*
 * axisType, atomX, atomY, atomZ parameters added to addParticle		*								
 * axisType, atomX, atomY, atomZ parameters added to setParticleParameters	*								
 * axisType, atomX, atomY, atomZ parameters added to getParticleParameters	*								
 * axisType, atomX, atomY, atomZ parameters added to class ParticleInfo		*								
 *										*
 * public axisTypes(ZthenX,Bisector,ZBisect,Threefold,ZOnly,NoAxisType) added	*
 *										*
 * Last Edit 07/2018 tjanicki							*
 * -----------------------------------------------------------------------------*/



#include "openmm/TabulatedFunction.h"
#include "openmm/Force.h"
#include "openmm/Vec3.h"
#include "openmm/Context.h"
#include <map>
#include <set>
#include <utility>
#include <vector>
#include "internal/windowsExportCustomAnisotropicNonbonded.h"

namespace CustomAnisotropicNonbondedPlugin {


class OPENMM_EXPORT_CUSTOMANISOTROPICNONBONDED CustomAnisotropicNonbondedForce : public OpenMM::Force {
public:
    /**
     * This is an enumeration of the different methods that may be used for handling long range nonbonded forces.
     */
	
    enum NonbondedMethod {
        /**
         * No cutoff is applied to nonbonded interactions.  The full set of N^2 interactions is computed exactly.
         * This necessarily means that periodic boundary conditions cannot be used.  This is the default.
         */
        NoCutoff = 0,
        /**
         * Interactions beyond the cutoff distance are ignored.
         */
        CutoffNonPeriodic = 1,
        /**
         * Periodic boundary conditions are used, so that each particle interacts only with the nearest periodic copy of
         * each other particle.  Interactions beyond the cutoff distance are ignored.
         */
        CutoffPeriodic = 2,
    };
    /*
 *	Define axis types
 *
 * */
    enum AxisTypes {
	ZThenX = 0,
	Bisector = 1,
	ZBisect = 2,
	ThreeFold = 3,
	ZOnly = 4,
	NoAxisType = 5,
	LastAxisTypeIndex = 6 };
    /**
     * Create a CustomAnisotropicNonbondedForce.
     *
     * @param energy    an algebraic expression giving the interaction energy between two particles as a function
     *                  of r, the distance between them, theta and phi (of local coordinate systems and using the
     *                  physics notation -- theta[0,2pi], phi[0,pi]), as well as any global and per-particle parameters
     */
    explicit CustomAnisotropicNonbondedForce(const std::string& energy);
    CustomAnisotropicNonbondedForce(const CustomAnisotropicNonbondedForce& rhs); // copy constructor
    ~CustomAnisotropicNonbondedForce();
    /**
     * Get the number of particles for which force field parameters have been defined.
     */
    int getNumParticles() const {
        return particles.size();
    }
    /**
     * Get the number of particle pairs whose interactions should be excluded.
     */
    int getNumExclusions() const {
        return exclusions.size();
    }
    /**
     * Get the number of per-particle parameters that the interaction depends on.
     */
    int getNumPerParticleParameters() const {
        return parameters.size();
    }
    /**
     * Get the number of global parameters that the interaction depends on.
     */
    int getNumGlobalParameters() const {
        return globalParameters.size();
    }
    /**
     * Get the number of tabulated functions that have been defined.
     */
    int getNumTabulatedFunctions() const {
        return functions.size();
    }
    /**
     * Get the number of tabulated functions that have been defined.
     *
     * @deprecated This method exists only for backward compatibility.  Use getNumTabulatedFunctions() instead.
     */
    int getNumFunctions() const {
        return functions.size();
    }
    /**
     * Get the number of interaction groups that have been defined.
    */
    int getNumInteractionGroups() const {
        return interactionGroups.size();
     }

    /**
     * Get the number of global parameters with respect to which the derivative of the energy
     * should be computed.
     */
    int getNumEnergyParameterDerivatives() const {
        return energyParameterDerivatives.size();
    }
    /**
     * Get the algebraic expression that gives the interaction energy between two particles
     */
    const std::string& getEnergyFunction() const;
    /**
     * Set the algebraic expression that gives the interaction energy between two particles
     */
    void setEnergyFunction(const std::string& energy);
    /**
     * Get the method used for handling long range nonbonded interactions.
     */
    NonbondedMethod getNonbondedMethod() const;
    /**
     * Set the method used for handling long range nonbonded interactions.
     */
    void setNonbondedMethod(NonbondedMethod method);
    /**
     * Get the cutoff distance (in nm) being used for nonbonded interactions.  If the NonbondedMethod in use
     * is NoCutoff, this value will have no effect.
     *
     * @return the cutoff distance, measured in nm
     */
    double getCutoffDistance() const;
    /**
     * Set the cutoff distance (in nm) being used for nonbonded interactions.  If the NonbondedMethod in use
     * is NoCutoff, this value will have no effect.
     *
     * @param distance    the cutoff distance, measured in nm
     */
    void setCutoffDistance(double distance);
    /**
     * Get whether a switching function is applied to the interaction.  If the nonbonded method is set
     * to NoCutoff, this option is ignored.
     */
    bool getUseSwitchingFunction() const;
    /**
     * Set whether a switching function is applied to the interaction.  If the nonbonded method is set
     * to NoCutoff, this option is ignored.
     */
    void setUseSwitchingFunction(bool use);
    /**
     * Get the distance at which the switching function begins to reduce the interaction.  This must be
     * less than the cutoff distance.
     */
    double getSwitchingDistance() const;
    /**
     * Set the distance at which the switching function begins to reduce the interaction.  This must be
     * less than the cutoff distance.
     */
    void setSwitchingDistance(double distance);
    /**
     * Add a new per-particle parameter that the interaction may depend on.
     *
     * @param name     the name of the parameter
     * @return the index of the parameter that was added
     */
    int addPerParticleParameter(const std::string& name);
    /**
     * Get the name of a per-particle parameter.
     *
     * @param index     the index of the parameter for which to get the name
     * @return the parameter name
     */
    const std::string& getPerParticleParameterName(int index) const;
    /**
     * Set the name of a per-particle parameter.
     *
     * @param index          the index of the parameter for which to set the name
     * @param name           the name of the parameter
     */
    void setPerParticleParameterName(int index, const std::string& name);
    /**
     * Add a new global parameter that the interaction may depend on.
     *
     * @param name             the name of the parameter
     * @param defaultValue     the default value of the parameter
     * @return the index of the parameter that was added
     */
    int addGlobalParameter(const std::string& name, double defaultValue);
    /**
     * Get the name of a global parameter.
     *
     * @param index     the index of the parameter for which to get the name
     * @return the parameter name
     */
    const std::string& getGlobalParameterName(int index) const;
    /**
     * Set the name of a global parameter.
     *
     * @param index          the index of the parameter for which to set the name
     * @param name           the name of the parameter
     */
    void setGlobalParameterName(int index, const std::string& name);
    /**
     * Get the default value of a global parameter.
     *
     * @param index     the index of the parameter for which to get the default value
     * @return the parameter default value
     */
    double getGlobalParameterDefaultValue(int index) const;
    /**
     * Set the default value of a global parameter.
     *
     * @param index          the index of the parameter for which to set the default value
     * @param defaultValue   the default value of the parameter
     */
    void setGlobalParameterDefaultValue(int index, double defaultValue);
    /**
     * Request that this Force compute the derivative of its energy with respect to a global parameter.
     * The parameter must have already been added with addGlobalParameter().
     *
     * @param name             the name of the parameter
     */
    void addEnergyParameterDerivative(const std::string& name);
    /**
     * Get the name of a global parameter with respect to which this Force should compute the
     * derivative of the energy.
     *
     * @param index     the index of the parameter derivative, between 0 and getNumEnergyParameterDerivatives()
     * @return the parameter name
     */
    const std::string& getEnergyParameterDerivativeName(int index) const;
    /**
     * Add the nonbonded force parameters for a particle.  This should be called once for each particle
	*include "openmm/internal/ForceImpl.h"
     * in the System.  When it is called for the i'th time, it specifies the parameters for the i'th particle.
     *
     * @param parameters    the list of parameters for the new particle
     * @param axisType	    integer value assigned to axis type
     * @param atomX, atomY, atomZ	symmetry axes; see docs for details
     * @return the index of the particle that was added
     */
    int addParticle(const std::vector<double>& parameters=std::vector<double>(), int axisType=5, int atomX=-1, int atomY=-1, int atomZ=-1);
    /**
     * Get the nonbonded force parameters for a particle.
     *
     * @param index            the index of the particle for which to get parameters
     * @param axisType	    integer value assigned to axis type
     * @param atomX, atomY, atomZ	symmetry axes; see docs for details
     * @param[out] parameters  the list of parameters for the specified particle
     */
    void getParticleParameters(int index, std::vector<double>& parameters, int& axisType, int& atomX, int& atomY, int& atomZ) const;
    /**
     * Set the nonbonded force parameters for a particle.
     *
     * @param index       the index of the particle for which to set parameters
     * @param axisType	    integer value assigned to axis type
     * @param atomX, atomY, atomZ	symmetry axes; see docs for details
     * @param parameters  the list of parameters for the specified particle
     */
   void setParticleParameters(int index, const std::vector<double>& parameters, int axisType, int atomX, int atomY, int atomZ);
    /**
     * Add a particle pair to the list of interactions that should be excluded.
     *
     * In many cases, you can use createExclusionsFromBonds() rather than adding each exclusion explicitly.
     *
     * @param particle1  the index of the first particle in the pair
     * @param particle2  the index of the second particle in the pair
     * @return the index of the exclusion that was added
     */
    int addExclusion(int particle1, int particle2);
    /**
     * Get the particles in a pair whose interaction should be excluded.
     *
     * @param index           the index of the exclusion for which to get particle indices
     * @param[out] particle1  the index of the first particle in the pair
     * @param[out] particle2  the index of the second particle in the pair
     */
    void getExclusionParticles(int index, int& particle1, int& particle2) const;
    /**
     * Set the particles in a pair whose interaction should be excluded.
     *
     * @param index      the index of the exclusion for which to set particle indices
     * @param particle1  the index of the first particle in the pair
     * @param particle2  the index of the second particle in the pair
     */
    void setExclusionParticles(int index, int particle1, int particle2);
    /**
     * Identify exclusions based on the molecular topology.  Particles which are separated by up to a specified number of
     * bonds are added as exclusions.
     *
     * @param bonds       the set of bonds based on which to construct exclusions.  Each element specifies the indices of
     *                    two particles that are bonded to each other.
     * @param bondCutoff  pairs of particles that are separated by this many bonds or fewer are added to the list of exclusions
     */
    void createExclusionsFromBonds(const std::vector<std::pair<int, int> >& bonds, int bondCutoff);
    /**
     * Add a tabulated function that may appear in the energy expression.
     *
     * @param name           the name of the function as it appears in expressions
     * @param function       a TabulatedFunction object defining the function.  The TabulatedFunction
     *                       should have been created on the heap with the "new" operator.  The
     *                       Force takes over ownership of it, and deletes it when the Force itself is deleted.
     * @return the index of the function that was added
     */
    int addTabulatedFunction(const std::string& name, OpenMM::TabulatedFunction* function);
    /**
     * Get a const reference to a tabulated function that may appear in the energy expression.
     *
     * @param index     the index of the function to get
     * @return the TabulatedFunction object defining the function
     */
    const OpenMM::TabulatedFunction& getTabulatedFunction(int index) const;
    /**
     * Get a reference to a tabulated function that may appear in the energy expression.
     *
     * @param index     the index of the function to get
     * @return the TabulatedFunction object defining the function
     */
    OpenMM::TabulatedFunction& getTabulatedFunction(int index);
    /**
     * Get the name of a tabulated function that may appear in the energy expression.
     *
     * @param index     the index of the function to get
     * @return the name of the function as it appears in expressions
     */
    const std::string& getTabulatedFunctionName(int index) const;
    /**
     * Add a tabulated function that may appear in the energy expression.
     *
     * @deprecated This method exists only for backward compatibility.  Use addTabulatedFunction() instead.
     */
    int addFunction(const std::string& name, const std::vector<double>& values, double min, double max);
    /**
     * Get the parameters for a tabulated function that may appear in the energy expression.
     *
     * @deprecated This method exists only for backward compatibility.  Use getTabulatedFunctionParameters() instead.
     * If the specified function is not a Continuous1DFunction, this throws an exception.
     */
    void getFunctionParameters(int index, std::string& name, std::vector<double>& values, double& min, double& max) const;
    /**
     * Set the parameters for a tabulated function that may appear in the energy expression.
     *
     * @deprecated This method exists only for backward compatibility.  Use setTabulatedFunctionParameters() instead.
     * If the specified function is not a Continuous1DFunction, this throws an exception.
     */
    void setFunctionParameters(int index, const std::string& name, const std::vector<double>& values, double min, double max);
     /**
     * Add an interaction group.  An interaction will be computed between every particle in set1 and every particle in set2.
     *      
     * @param set1    the first set of particles forming the interaction group
     * @param set2    the second set of particles forming the interaction group
     * @return the index of the interaction group that was added
     */
    int addInteractionGroup(const std::set<int>& set1, const std::set<int>& set2);
    /**
     * Get the parameters for an interaction group.
     *
     * @param index        the index of the interaction group for which to get parameters
     * @param[out] set1    the first set of particles forming the interaction group
     * @param[out] set2    the second set of particles forming the interaction group
     */
    void getInteractionGroupParameters(int index, std::set<int>& set1, std::set<int>& set2) const;
    /**
     * Set the parameters for an interaction group.
     *
     * @param index   the index of the interaction group for which to set parameters
     * @param set1    the first set of particles forming the interaction group
     * @param set2    the second set of particles forming the interaction group
     */
    void setInteractionGroupParameters(int index, const std::set<int>& set1, const std::set<int>& set2);
    /**
     * Update the per-particle parameters in a Context to match those stored in this Force object.  This method provides
     * an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
     * Simply call setParticleParameters() to modify this object's parameters, then call updateParametersInContext()
     * to copy them over to the Context.
     *
     * This method has several limitations.  The only information it updates is the values of per-particle parameters.
     * All other aspects of the Force (the energy function, nonbonded method, cutoff distance, etc.) are unaffected and can
     * only be changed by reinitializing the Context.  Also, this method cannot be used to add new particles, only to change
     * the parameters of existing ones.
     */
    void updateParametersInContext(OpenMM::Context& context);
    /**
     * Returns whether or not this force makes use of periodic boundary
     * conditions.
     *
     * @returns true if force uses PBC and false otherwise
     */
    bool usesPeriodicBoundaryConditions() const {
        return nonbondedMethod == CutoffPeriodic;
    }

/*	int addExpression(const std::string& expression, const char  var);
	void getExpression(int index, std::string& expression, char var);
	void setExpression(int index, const std::string& expression, const char var);
*/
protected:
    OpenMM::ForceImpl* createImpl() const;
private:
    // REMEMBER TO UPDATE THE COPY CONSTRUCTOR IF YOU ADD ANY NEW FIELDS !!
    class ParticleInfo;
    class PerParticleParameterInfo;
    class GlobalParameterInfo;
    class ExpressionInfo;
    class ExclusionInfo;
    class FunctionInfo;
    class InteractionGroupInfo;
    NonbondedMethod nonbondedMethod;
    double cutoffDistance, switchingDistance;
    bool useSwitchingFunction;
    bool useYlm;
    std::string energyExpression;
    char var;
    std::vector<PerParticleParameterInfo> parameters;
    std::vector<GlobalParameterInfo> globalParameters;
    std::vector<ParticleInfo> particles;
    std::vector<ExclusionInfo> exclusions;
    std::vector<FunctionInfo> functions;
    std::vector<InteractionGroupInfo> interactionGroups;
    std::vector<int> energyParameterDerivatives;
};

/**
 * This is an internal class used to record information about a particle.
 * @private
 */
class CustomAnisotropicNonbondedForce::ParticleInfo {
public:
    std::vector<double> parameters;
    int axisType, atomX, atomY, atomZ;
    ParticleInfo() {
	axisType = atomZ = atomX = atomY = -1;
    }
    ParticleInfo(const std::vector<double>& parameters, int axisType, int atomX, int atomY, int atomZ) : parameters(parameters), axisType(axisType), atomX(atomX), atomY(atomY), atomZ(atomZ) {
    }
};

/**
 * This is an internal class used to record information about a per-particle parameter.
 * @private
 */
class CustomAnisotropicNonbondedForce::PerParticleParameterInfo {
public:
    std::string name;
    PerParticleParameterInfo() {
    }
    PerParticleParameterInfo(const std::string& name) : name(name) {
    }
};

/**
 * This is an internal class used to record information about a global parameter.
 * @private
 */
class CustomAnisotropicNonbondedForce::GlobalParameterInfo {
public:
    std::string name;
    double defaultValue;
    GlobalParameterInfo() {
    }
    GlobalParameterInfo(const std::string& name, double defaultValue) : name(name), defaultValue(defaultValue) {
    }
};

/**
 * This is an internal class used to record information about an exclusion.
 * @private
 */
class CustomAnisotropicNonbondedForce::ExclusionInfo {
public:
    int particle1, particle2;
    ExclusionInfo() {
        particle1 = particle2 = -1;
    }
    ExclusionInfo(int particle1, int particle2) :
        particle1(particle1), particle2(particle2) {
    }
};

/**
 * This is an internal class used to record information about a tabulated function.
 * @private
 */
class CustomAnisotropicNonbondedForce::FunctionInfo {
public:
    std::string name;
    OpenMM::TabulatedFunction* function;
    FunctionInfo() {
    }
    FunctionInfo(const std::string& name, OpenMM::TabulatedFunction* function) : name(name), function(function) {
    }
};

/**
 * This is an internal class used to record information about an interaction group.
 * @private
 */
class CustomAnisotropicNonbondedForce::InteractionGroupInfo {
public:
    std::set<int> set1, set2;
    InteractionGroupInfo() {
    }
    InteractionGroupInfo(const std::set<int>& set1, const std::set<int>& set2) :
        set1(set1), set2(set2) {
    }
};
} // namespace OpenMM

#endif /*OPENMM_CUSTOMANISOTROPICNONBONDEDFORCE_H_*/

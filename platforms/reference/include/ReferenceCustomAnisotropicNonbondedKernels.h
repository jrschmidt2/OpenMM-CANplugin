#ifndef REFERENCE_CUSTOMANISOTROPICNONBONDED_KERNELS_H_
#define REFERENCE_CUSTOMANISOTROPICNONBONDED_KERNELS_H_

#include "CustomAnisotropicNonbondedKernels.h"

#include "openmm/reference/SimTKOpenMMRealType.h"
#include "openmm/reference/ReferenceNeighborList.h"
#include "openmm/reference/ReferencePlatform.h"

#include "lepton/CompiledExpression.h"
#include "lepton/CustomFunction.h"

namespace CustomAnisotropicNonbondedPlugin {

class ReferenceCalcCustomAnisotropicNonbondedForceKernel : public CalcCustomAnisotropicNonbondedForceKernel {
public:
	ReferenceCalcCustomAnisotropicNonbondedForceKernel(std::string name, const OpenMM::Platform& platform) : CalcCustomAnisotropicNonbondedForceKernel(name, platform), forceCopy(NULL) {
	}
	~ReferenceCalcCustomAnisotropicNonbondedForceKernel();
	void initialize(const OpenMM::System& system, const CustomAnisotropicNonbondedForce& force);
	double execute(OpenMM::ContextImpl& context, bool includeForces, bool includeEnergy);
	void copyParametersToContext(OpenMM::ContextImpl& context, const CustomAnisotropicNonbondedForce& force);
private:
	int numParticles;
	double **particleParamArray;
	int *atomZs;
	int *atomXs;
	int *atomYs;
	int *axisTypes;
	int *calc;
	int YLM;
	double nonbondedCutoff, switchingDistance, periodicBoxSize[3], longRangeCoefficient;
	bool useSwitchingFunction, hasInitializedLongRangeCorrection, isPeriodic;
	Lepton::CompiledExpression EExpression, forceExpR, forceExpTheta1, forceExpTheta2, forceExpPhi1, forceExpPhi2;
	Lepton::CompiledExpression forceExpY101;
	Lepton::CompiledExpression forceExpY102;
	Lepton::CompiledExpression forceExpY11c1;
	Lepton::CompiledExpression forceExpY11c2;
	Lepton::CompiledExpression forceExpY11s1;
	Lepton::CompiledExpression forceExpY11s2;
	Lepton::CompiledExpression forceExpY201;
	Lepton::CompiledExpression forceExpY202;
	Lepton::CompiledExpression forceExpY21c1;
	Lepton::CompiledExpression forceExpY21c2;
	Lepton::CompiledExpression forceExpY21s1;
	Lepton::CompiledExpression forceExpY21s2;
	Lepton::CompiledExpression forceExpY22c1;
	Lepton::CompiledExpression forceExpY22c2;
	Lepton::CompiledExpression forceExpY22s1;
	Lepton::CompiledExpression forceExpY22s2;
	CustomAnisotropicNonbondedForce* forceCopy;
	std::map<std::string, double> globalParamValues;
	std::vector<std::set<int> > exclusions;
	std::vector<std::string> parameterNames, globalParameterNames, energyParamDerivNames;
	std::vector<std::pair<std::set<int>, std::set<int> > > interactionGroups;
	std::vector<double> longRangeCoefficientDerivs;
	NonbondedMethod nonbondedMethod;
	OpenMM::NeighborList* neighborList;
};
}
#endif

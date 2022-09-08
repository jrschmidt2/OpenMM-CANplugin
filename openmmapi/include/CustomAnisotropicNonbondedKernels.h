#ifndef CUSTOMANISOTROPICNONBONDED_KERNELS_H_
#define CUSTOMANISOTROPICNONBONDED_KERNELS_H_

#include "openmm/KernelImpl.h"
#include "openmm/System.h"
#include "openmm/Platform.h"
#include "CustomAnisotropicNonbondedForce.h"
#include <string>

namespace CustomAnisotropicNonbondedPlugin {

/* 
 * Called by CustomAnisotropicNonbondedForce
 * This kernel is called to calculate forces acting on and energy of the system
 * cr. tjanicki 05/2017
 * ed. tjanicki 06/2017
 *
 */


class CalcCustomAnisotropicNonbondedForceKernel : public OpenMM::KernelImpl {
public:
	enum NonbondedMethod {
		NoCutoff = 0,
		CutoffNonPeriodic = 1,
		CutoffPeriodic = 2
	};

	static std::string Name() {
		return "CalcCustomAnisotropicNonbondedForce";
	}

	CalcCustomAnisotropicNonbondedForceKernel(std::string name, const OpenMM::Platform& platform) : OpenMM::KernelImpl(name, platform) {
	}

	virtual void initialize(const OpenMM::System& system, const CustomAnisotropicNonbondedForce& force)=0;

	virtual double execute(OpenMM::ContextImpl& context, bool includeForces, bool includeEnergy)=0;

	virtual void copyParametersToContext(OpenMM::ContextImpl& context, const CustomAnisotropicNonbondedForce& force)=0;
};
}
#endif

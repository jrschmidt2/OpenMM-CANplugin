#ifndef OPENMM_REFERENCECUSTOMANISOTROPICNONBONDEDKERNELFACTORY_H_
#define OPENMM_REFERENCECUSTOMANISOTROPICNONBONDEDKERNELFACTORY_H_

#include "openmm/KernelFactory.h"

namespace OpenMM {

/*
 * Creates kernels for reference implementation of plugin
 * cr. tjanicki 05/2017
 * ed. tjanicki 06/2017
 *
 */

class ReferenceCustomAnisotropicNonbondedKernelFactory : public KernelFactory {
public:
	KernelImpl* createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const;
};
}
#endif


#include "ReferenceCustomAnisotropicNonbondedKernelFactory.h"
#include "ReferenceCustomAnisotropicNonbondedKernels.h"
#include "openmm/reference/ReferencePlatform.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"

using namespace OpenMM;
using namespace CustomAnisotropicNonbondedPlugin;

extern "C" OPENMM_EXPORT void registerPlatforms() {
}

extern "C" OPENMM_EXPORT void registerKernelFactories() {
   for (int i=0; i<Platform::getNumPlatforms(); i++) {
      Platform& platform = Platform::getPlatform(i);
      if (dynamic_cast<ReferencePlatform*>(&platform) != NULL) {
         ReferenceCustomAnisotropicNonbondedKernelFactory* factory = new ReferenceCustomAnisotropicNonbondedKernelFactory();
         platform.registerKernelFactory(CalcCustomAnisotropicNonbondedForceKernel::Name(),factory);
      }
   }
}

extern "C" OPENMM_EXPORT void registerCustomAnisotropicNonbondedReferenceKernelFactories() {
   registerKernelFactories();
}

KernelImpl* ReferenceCustomAnisotropicNonbondedKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const { 
   ReferencePlatform::PlatformData& data = *static_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
   if (name == CalcCustomAnisotropicNonbondedForceKernel::Name())
      return new ReferenceCalcCustomAnisotropicNonbondedForceKernel(name, platform);
   throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '")+name+"'").c_str());
}


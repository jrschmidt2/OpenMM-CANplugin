#ifdef WIN32
#include <windows.h>
#include <sstream>
#else
#include <dlfcn.h>
#include <dirent.h>
#include <cstdlib>
#endif

#include "CustomAnisotropicNonbondedForce.h"
#include "CustomAnisotropicNonbondedForceProxy.h"
#include "openmm/serialization/SerializationProxy.h"

#if defined(WIN32)
    #include <windows.h>
    extern "C" OPENMM_EXPORT_CAN void registerCustomAnisotropicNonbondedSerializationProxies();
    BOOL WINAPI DllMain(HANDLE hModule, DWORD  ul_reason_for_call, LPVOID lpReserved) {
        if (ul_reason_for_call == DLL_PROCESS_ATTACH)
            registerCustomAnisotropicNonbondedSerializationProxies();
        return TRUE;
    }
#else
    extern "C" void __attribute__((constructor)) registerCustomAnisotropicNonbondedSerializationProxies();
#endif

using namespace CustomAnisotropicNonbondedPlugin;
using namespace OpenMM;

extern "C" OPENMM_EXPORT_CAN void registerCustomAnisotropicNonbondedSerializationProxies() {
    SerializationProxy::registerProxy(typeid(CustomAnisotropicNonbondedForce), new CustomAnisotropicNonbondedForceProxy());
}

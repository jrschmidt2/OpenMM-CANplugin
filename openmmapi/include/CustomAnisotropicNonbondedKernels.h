#ifndef CUSTOMANISOTROPICNONBONDED_KERNELS_H_
#define CUSTOMANISOTROPICNONBONDED_KERNELS_H_
/*--------------------------------------------------------------------*
*                   OpenMM CustomAnisotropicNonbondedPlugin           *
*---------------------------------------------------------------------*
*                                                                     *
* This is part of the OpenMM molecular simulation toolkit originating *
* from Simbios, the NIH National Center for Physics-Based Simulation  *
* of Biological Structures at Stanford, funded under the NIH Roadmap  *
* for Medical Research, grant U54 GM072970. See https://simtk.org.    *
*                                                                     *
* Portions copyright (c) 2014-2021 Stanford University and the        *
* Authors.                                                            *
*                                                                     *
* Authors: Peter Eastman                                              *
*                                                                     *
* Contributors: Tesia D. Janicki, Mary J. Van Vleet                   *
*                                                                     *
* Permission is hereby granted, free of charge, to any person         *
* obtaining a copy of this software and associated documentation      *
* files (the "Software"), to deal in the Software without restriction,*
* including without limitation the rights to use, copy, modify, merge,*
* publish, distribute, sublicense, and/or sell copies of the Software,*
* and to permit persons to whom the Software is furnished to do so,   *
* subject to the following conditions:                                *
*                                                                     * 
* The above copyright notice and this permission notice shall be      *
* included in all copies or substantial portions of the Software.     *
*                                                                     *
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,     *
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF  *
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND               *
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR     *
* COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER         *
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,     *
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE  *
* OR OTHER DEALINGS IN THE SOFTWARE.                                  *
*---------------------------------------------------------------------*/

#include "openmm/KernelImpl.h"
#include "openmm/System.h"
#include "openmm/Platform.h"
#include "CustomAnisotropicNonbondedForce.h"
#include <string>

namespace CustomAnisotropicNonbondedPlugin {
/* 
 * Called by CustomAnisotropicNonbondedForce
 * This kernel is called to calculate forces acting on and energy of the system
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

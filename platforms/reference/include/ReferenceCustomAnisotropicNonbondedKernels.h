#ifndef REFERENCE_CAN_KERNELS_H_
#define REFERENCE_CAN_KERNELS_H_
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
    double nonbondedCutoff, switchingDistance, periodicBoxSize[3], longRangeCoefficient;
    bool useSwitchingFunction, hasInitializedLongRangeCorrection, isPeriodic;
    Lepton::CompiledExpression EExpression, forceExpR, forceExpTheta1, forceExpTheta2, forceExpPhi1, forceExpPhi2;
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

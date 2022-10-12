#ifdef WIN32
#define _USE_MATH_DEFINES
#endif

#include "CudaCustomAnisotropicNonbondedKernels.h"
#include "CudaCustomAnisotropicNonbondedKernelSources.h"
#include "openmm/internal/ContextImpl.h"
#include "internal/CustomAnisotropicNonbondedForceImpl.h"
#include "openmm/cuda/CudaExpressionUtilities.h"
#include "openmm/cuda/CudaForceInfo.h"
#include "openmm/cuda/CudaNonbondedUtilities.h"
#include "lepton/ExpressionTreeNode.h"
#include "lepton/Parser.h"
#include "lepton/Operation.h"
#include <cmath>
#include <tgmath.h>
#include <iterator>
#include <algorithm>
#include <iostream>
#ifdef _MSC_VER
#include <windows.h>
#endif

using namespace CustomAnisotropicNonbondedPlugin;
using namespace OpenMM;
using namespace std;

class CudaCalcCustomAnisotropicNonbondedForceKernel::ForceInfo : public CudaForceInfo {
public:
    ForceInfo(const CustomAnisotropicNonbondedForce& force) : force(force) {
    }
    bool areParticlesIdentical(int particle1, int particle2) {
        int axisType1,axisType2,atomX1,atomX2,atomY1,atomY2,atomZ1,atomZ2;
        vector<double> params1, params2;
        force.getParticleParameters(particle1,params1,axisType1,atomX1,atomY1,atomZ1);
        force.getParticleParameters(particle2,params2,axisType2,atomX2,atomY2,atomZ2);
        for (int i = 0; i < (int) params1.size(); i++)
            if (params1[i] != params2[i]) return false;
        if (atomX1 != atomX2 || atomY1 != atomY2 || atomZ1 != atomZ2 || axisType1 != axisType2) return false;
        return true;
    }
    int getNumParticleGroups() {
        return force.getNumExclusions();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        int particle1,particle2;
        force.getExclusionParticles(index, particle1, particle2);
        particles.resize(2);
        particles[0] = particle1;
        particles[1] = particle2;
    }
    bool areGroupsIdentical(int group1, int group2) {
	return true;
    }

private:
    const CustomAnisotropicNonbondedForce& force;
};

CudaCalcCustomAnisotropicNonbondedForceKernel::~CudaCalcCustomAnisotropicNonbondedForceKernel() {
    cu.setAsCurrent();
//    if (params != NULL) delete params;
//    cu.clearBuffer(kvecs);
}

void CudaCalcCustomAnisotropicNonbondedForceKernel::initialize(const System& system, const CustomAnisotropicNonbondedForce& force) {

    cu.setAsCurrent();
    int forceIndex;
    for (forceIndex = 0; forceIndex < system.getNumForces() && &system.getForce(forceIndex) != &force; ++forceIndex);
    string prefix = (force.getNumInteractionGroups() ==0 ? "custom"+cu.intToString(forceIndex)+"_" : "");
    //cutoffs
    bool useCutoff = (force.getNonbondedMethod() != CustomAnisotropicNonbondedForce::NoCutoff);
    bool usePeriodic = (force.getNonbondedMethod() != CustomAnisotropicNonbondedForce::NoCutoff && force.getNonbondedMethod() != CustomAnisotropicNonbondedForce::CutoffNonPeriodic);
    //setup space
    numParticles = force.getNumParticles();
    numPaddedAtoms = cu.getPaddedNumAtoms();
    axes.initialize<int4>(cu,numPaddedAtoms,"axes"); 
    //axes.initialize<int4>(cu,numParticles,"axes"); 
    kvecs.initialize(cu,3*numPaddedAtoms,(cu.getUseDoublePrecision() ? sizeof(double4) : sizeof(float4)),"kvecs");
    //kvecs.initialize(cu,3*numParticles,(cu.getUseDoublePrecision() ? sizeof(double4) : sizeof(float4)),"kvecs");
    params = new CudaParameterSet(cu, force.getNumPerParticleParameters(), numParticles, "customAnisotropicNonbondedParameters");
    if (force.getNumGlobalParameters() > 0)
        globals.initialize(cu,force.getNumGlobalParameters(),cu.getUseDoublePrecision() ? sizeof(double) : sizeof(float), "customAnisotropicNonbondedGlobals");
    vector<vector<float> > paramVector(numParticles);
    vector<int4> axesVector;
    //vector<int4> axesVector(numParticles);
    //create atomic parameter list, assign axes
    for (int i = 0; i < numParticles; i++) {
        vector<double> parameters;
        int axisType, atomX, atomY, atomZ;
        force.getParticleParameters(i,parameters,axisType,atomX,atomY,atomZ);
        axesVector.push_back(make_int4(atomX,atomY,atomZ,axisType));
        paramVector[i].resize(parameters.size());
        for (int j = 0; j < (int) parameters.size(); j++) {
            paramVector[i][j] = (float) parameters[j];
	}
    }
    for (int i = numParticles; i < numPaddedAtoms; i++) {
        axesVector.push_back(make_int4(0,0,0,0));
    }
    params->setParameterValues(paramVector);
    axes.upload(axesVector);
    //handle exclusions
    if (force.getNumExclusions() > 0) {
        vector<vector<int>> exclusionList(numParticles);
        for (int i = 0; i < force.getNumExclusions(); i++) {
            int particle1, particle2;
            force.getExclusionParticles(i,particle1,particle2);
            exclusionList[particle1].push_back(particle2);
            exclusionList[particle2].push_back(particle1);
        }
        vector<int> exclusionsVec;
        vector<int> exclusionStartIdxVec(numParticles+1);
        for (int i = 0; i < numParticles; i++) {
            sort(exclusionList[i].begin(),exclusionList[i].end());
            exclusionsVec.insert(exclusionsVec.end(), exclusionList[i].begin(),exclusionList[i].end());
            exclusionStartIdxVec[i+1] = exclusionsVec.size();
        }
        exclusions.initialize<int>(cu,exclusionsVec.size(),"customAnisotropicNonbondedExclusions");
        exclusionStartIdx.initialize<int>(cu,exclusionStartIdxVec.size(),"customAnisotropicNonbondedExclusions");
        exclusions.upload(exclusionsVec);
        exclusionStartIdx.upload(exclusionStartIdxVec);
    }
    //set tabulated functions
    map<string, Lepton::CustomFunction*> functions;
    vector<pair<string,string> > functionDefinitions;
    vector<const TabulatedFunction*> functionList;
    stringstream tableArgs;
    tabulatedFunctions.resize(force.getNumTabulatedFunctions());
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        functionList.push_back(&force.getTabulatedFunction(i));
        string name = force.getTabulatedFunctionName(i);
        string arrayName = prefix+"table"+cu.intToString(i);
        functionDefinitions.push_back(make_pair(name,arrayName));
        functions[name] = cu.getExpressionUtilities().getFunctionPlaceholder(force.getTabulatedFunction(i));
        int width;
        vector<float> f = cu.getExpressionUtilities().computeFunctionCoefficients(force.getTabulatedFunction(i),width);
        tabulatedFunctions[i].initialize<float>(cu,f.size(),"TabulatedFunction");
//        cu.addAutoclearBuffer(tabulatedFunctions[i]);
        tabulatedFunctions[i].upload(f);
        cu.getNonbondedUtilities().addArgument(CudaNonbondedUtilities::ParameterInfo(arrayName,"float",width,width*sizeof(float),tabulatedFunctions[i].getDevicePointer()));
        tableArgs << ", const float";
        if (width > 1) tableArgs << width;
        else tableArgs << "* __restrict__ "<< arrayName;
    }
    //create parameter list
    globalParamNames.resize(force.getNumGlobalParameters());
    globalParamValues.resize(force.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParamNames[i] = force.getGlobalParameterName(i);
        globalParamValues[i] = (float) force.getGlobalParameterDefaultValue(i);
    }
    if (globals.isInitialized()) globals.upload(globalParamValues);

    //parse force expressions
    Lepton::ParsedExpression expression = Lepton::Parser::parse(force.getEnergyFunction(),functions).optimize();
    map<string, Lepton::ParsedExpression> forceExpressions;
    forceExpressions["dEdR = "] = expression.differentiate("r").optimize();
    forceExpressions["dEdTheta1 = "] = expression.differentiate("theta1").optimize();
    forceExpressions["dEdTheta2 = "] = expression.differentiate("theta2").optimize();
    forceExpressions["dEdPhi1 = "] = expression.differentiate("phi1").optimize();
    forceExpressions["dEdPhi2 = "] = expression.differentiate("phi2").optimize();
    forceExpressions["energy += "] = expression;

    stringstream compute;
    //set variables
    map<string,string> variables;
    variables["r"] = "r"; 
    variables["theta1"] = "computeAzim(-kvec1z,-rij)";
    variables["theta2"] = "computeAzim(-kvec2z,rij)";
    variables["phi1"] = "computePolar(kvec1x,kvec1z,kvec1z-rij)";
    variables["phi2"] = "computePolar(kvec2x,kvec2z,kvec2z+rij)";
    for (int i = 0; i < force.getNumPerParticleParameters(); i++) {
        const string& name = force.getPerParticleParameterName(i);
        variables[name+"1"] = prefix+"params"+params->getParameterSuffix(i,"1");
        variables[name+"2"] = prefix+"params"+params->getParameterSuffix(i,"2");
    }
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        const string& name = force.getGlobalParameterName(i);
        variables[name] = "globals["+cu.intToString(i)+"]";
    }

    stringstream extraArgs;
    if(force.getNumGlobalParameters() > 0) extraArgs << ", const float* __restrict__ globals";
    for (int i = 0; i < (int) params->getBuffers().size(); i++) {
        CudaNonbondedUtilities::ParameterInfo& buffer = params->getBuffers()[i];
        extraArgs << ", const " + buffer.getType()+"* __restrict__ "+buffer.getName();
        compute << buffer.getType()+" "+prefix+"params"+cu.intToString(i+1)+"1 = "+buffer.getName()+"[ii];\n";
        compute << buffer.getType()+" "+prefix+"params"+cu.intToString(i+1)+"2 = "+buffer.getName()+"[jj];\n";
    }

    //make kernels
    compute << cu.getExpressionUtilities().createExpressions(forceExpressions, variables, functionList, functionDefinitions,prefix+"temp");
    map<string,string> replacements;
    replacements["COMPUTE_FORCE"] = compute.str();
    replacements["PARAMETER_ARGUMENTS"] = extraArgs.str()+tableArgs.str();
    map<string,string> defines;
    defines["M_PI"] = cu.doubleToString(M_PI);
    defines["NUM_ATOMS"] = cu.intToString(numParticles);
    defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
    if (useCutoff) {
        defines["USE_CUTOFF"] = "1";
        defines["CUTOFF_SQUARED"] = cu.doubleToString(force.getCutoffDistance()*force.getCutoffDistance());
    }
    if (usePeriodic) defines["USE_PERIODIC"] = "1";
    if (force.getNumExclusions() > 0) defines["USE_EXCLUSIONS"] = "1";
    string source = cu.replaceStrings(CudaCustomAnisotropicNonbondedKernelSources::vectorOps+CudaCustomAnisotropicNonbondedKernelSources::customAnisotropicNonbonded, replacements);
    CUmodule module = cu.createModule(source,defines);
    computeAxesKernel = cu.getKernel(module,"accessAxisParameter");
    mykernel = cu.getKernel(module, "computeCAN");
    cu.getNonbondedUtilities().setUsePadding(false);
    cu.addForce(new ForceInfo(force));

}

double CudaCalcCustomAnisotropicNonbondedForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    if (numParticles == 0) return 0.0;
    if (globals.isInitialized()) {
        bool changed = false;
        for (int i = 0; i < (int) globalParamNames.size(); i++) {
            float value = (float) context.getParameter(globalParamNames[i]);
            if (value != globalParamValues[i]) changed = true;
            globalParamValues[i] = value;
        }
        if (changed) {
            globals.upload(globalParamValues);
        } 
    }
    if (!hasInitializedKernel) {
        hasInitializedKernel=true;
        //add args for axes kernel
        axesargs.push_back(&cu.getPosq().getDevicePointer());
        axesargs.push_back(&axes.getDevicePointer());
        axesargs.push_back(&kvecs.getDevicePointer());
        //add all args for kernel
        args.push_back(&cu.getForce().getDevicePointer());
        args.push_back(&cu.getEnergyBuffer().getDevicePointer());
        args.push_back(&cu.getPosq().getDevicePointer());
        args.push_back(&axes.getDevicePointer());
        args.push_back(&kvecs.getDevicePointer());
        args.push_back(cu.getPeriodicBoxSizePointer());
        args.push_back(cu.getInvPeriodicBoxSizePointer());
        args.push_back(cu.getPeriodicBoxVecXPointer());
        args.push_back(cu.getPeriodicBoxVecYPointer());
        args.push_back(cu.getPeriodicBoxVecZPointer());
        if (exclusions.isInitialized()) {
            args.push_back(&exclusions.getDevicePointer());
            args.push_back(&exclusionStartIdx.getDevicePointer());
        }
        if (globals.isInitialized()) args.push_back(&globals.getDevicePointer());
        for (auto& buffer : params->getBuffers()) args.push_back(&buffer.getMemory());
        for (auto& function : tabulatedFunctions) args.push_back(&function.getDevicePointer());
    }
    cu.executeKernel(computeAxesKernel,&axesargs[0],numParticles);
    //int sharedMemorySize = 3*CudaContext::ThreadBlockSize*sizeof(float4);
    //cu.executeKernel(mykernel,&args[0],numParticles,CudaContext::ThreadBlockSize, sharedMemorySize);
    int blockSize = 256;
    int sharedMemorySize = 3*blockSize*sizeof(float4);
    cu.executeKernel(mykernel,&args[0],numParticles,blockSize, sharedMemorySize);
    return 0.0;
}

void CudaCalcCustomAnisotropicNonbondedForceKernel::copyParametersToContext(ContextImpl& context, const CustomAnisotropicNonbondedForce& force) {

    cu.setAsCurrent();
    int numContext = cu.getPlatformData().contexts.size();
    int startIndex = cu.getContextIndex()*force.getNumParticles()/numContext;
    int endIndex = (cu.getContextIndex()+1)*force.getNumParticles()/numContext;
    if (numParticles != endIndex-startIndex) throw OpenMMException("updateParametersInContext: The number of particles has changed");
    vector<int4> axesVector(numParticles);
    if (numParticles > 0) {
        vector<vector<float> > paramVector(numParticles);
        vector<double> parameters;
        for (int i = 0; i < numParticles; i++) {
            force.getParticleParameters(startIndex+i,parameters,axesVector[i].w,axesVector[i].x,axesVector[i].y,axesVector[i].z);
            paramVector[i].resize(parameters.size());
            for (int j = 0; j < (int) parameters.size(); j++) {
                paramVector[i][j] = (float) parameters[j];
	    }
        }
        axes.upload(axesVector);
        params->setParameterValues(paramVector);
    }
    cu.invalidateMolecules();

}

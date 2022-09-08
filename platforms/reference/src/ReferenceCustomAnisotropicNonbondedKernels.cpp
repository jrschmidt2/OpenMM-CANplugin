#include "ReferenceCustomAnisotropicNonbondedKernels.h"
#include "internal/CustomAnisotropicNonbondedForceImpl.h"
#include "ReferenceCustomAnisotropicNonbondedIxn.h"

#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/reference/ReferenceConstraints.h"
#include "openmm/reference/ReferenceTabulatedFunction.h"

#include "lepton/Operation.h"
#include "lepton/ParsedExpression.h"
#include "lepton/CustomFunction.h"
#include "lepton/Parser.h"


using OpenMM::Vec3;
using OpenMM::ContextImpl;
using OpenMM::System;
using OpenMM::NeighborList;
using OpenMM::ReferencePlatform;
using OpenMM::OpenMMException;

using namespace CustomAnisotropicNonbondedPlugin;
using namespace std;
using std::string;

//ALLOCATION
static vector<Vec3>& extractPositions(ContextImpl& context) {
	ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
	return *((vector<Vec3>*) data->positions);
}

static vector<Vec3>& extractForces(ContextImpl& context) {
	ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
	return *((vector<Vec3>*) data->forces);
}

static Vec3* extractBoxVectors(ContextImpl& context) {
	ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
	return (Vec3*) data->periodicBoxVectors;
}
static double** allocateRealArray(int length, int width) {
	double** array= new double*[length];
	for (int i = 0; i < length; ++i)
		array[i] = new double[width];
	return array;
}

static void disposeRealArray(double** array, int size) {
	if (array) {
		for (int i = 0; i < size; ++i)
			delete[] array[i];
		delete[] array;
	}
}

static void validateVariables(const Lepton::ExpressionTreeNode& node, const set<string>& variables) {
	const Lepton::Operation& op = node.getOperation();
	if (op.getId() == Lepton::Operation::VARIABLE && variables.find(op.getName()) == variables.end())
		throw OpenMMException("Unknown variable in expression: "+ op.getName());
	for (auto& child : node.getChildren())
		validateVariables(child, variables);
}

//DESTRUCTOR
ReferenceCalcCustomAnisotropicNonbondedForceKernel::~ReferenceCalcCustomAnisotropicNonbondedForceKernel() {
	disposeRealArray(particleParamArray, numParticles);
	delete[] axisTypes;
	delete[] atomXs;
	delete[] atomYs;
	delete[] atomZs;
	delete[] calc;

	if(neighborList != NULL)
		delete neighborList;
	if(forceCopy != NULL)
		delete forceCopy;
}

//INITIALIZE
void ReferenceCalcCustomAnisotropicNonbondedForceKernel::initialize(const System& system, const CustomAnisotropicNonbondedForce& force) {

	numParticles = force.getNumParticles();
	int numParameters = force.getNumPerParticleParameters();
	
	//Exclusions
	exclusions.resize(numParticles);
	for (int i = 0; i < force.getNumExclusions(); i++) {
		int particle1, particle2;
		force.getExclusionParticles(i, particle1, particle2);
		exclusions[particle1].insert(particle2);
		exclusions[particle2].insert(particle1);
	}
	//Parameters & Methods
	atomXs = new int[numParticles];
	atomYs = new int[numParticles];
	atomZs = new int[numParticles];
	axisTypes = new int[numParticles];	
	particleParamArray = allocateRealArray(numParticles, numParameters);

	for (int i = 0; i < numParticles; ++i) {
		vector<double> parameters;
		int axisType,atomZ,atomX,atomY;
		force.getParticleParameters(i, parameters,axisType,atomX,atomY,atomZ);
		for (int j = 0; j < numParameters; j++) 
			particleParamArray[i][j] = parameters[j];
		axisTypes[i] = axisType;
		atomXs[i] = atomX;
		atomYs[i] = atomY;
		atomZs[i] = atomZ;
	}
	
	nonbondedMethod = CalcCustomAnisotropicNonbondedForceKernel::NonbondedMethod(force.getNonbondedMethod());
	nonbondedCutoff = force.getCutoffDistance();
	if (nonbondedMethod == NoCutoff) {
		neighborList = NULL;
		useSwitchingFunction = false;
	}
	else {
		neighborList = new NeighborList();
		useSwitchingFunction = force.getUseSwitchingFunction();
		switchingDistance = force.getSwitchingDistance();
	}

	//Function Preparation
	string inputExp = force.getEnergyFunction();

	map<string, Lepton::CustomFunction*> functions;
	for (int i = 0; i < force.getNumFunctions(); i++)
		functions[force.getTabulatedFunctionName(i)] = createReferenceTabulatedFunction(force.getTabulatedFunction(i));

	Lepton::ParsedExpression tempR,tempTheta1,tempTheta2,tempPhi1,tempPhi2,expression;
	Lepton::ParsedExpression tempY101,tempY102,tempY11c1,tempY11c2,tempY201,tempY202,tempY21c1,tempY21c2,tempY22c1,tempY22c2;
	Lepton::ParsedExpression tempY11s1,tempY11s2,tempY21s1,tempY21s2,tempY22s1,tempY22s2;
	expression = Lepton::Parser::parse(inputExp,functions).optimize();
	calc = new int[4];
	int counter = 0; 
	if (inputExp.find("theta1") != -1) {
		calc[0] = 0;
		tempTheta1 = expression.differentiate("theta1");
		forceExpTheta1 = tempTheta1.createCompiledExpression();
		counter += 1;
	}
	else calc[0] = 1;
	if (inputExp.find("theta2") != -1) {
		calc[1] = 0;
		tempTheta2 = expression.differentiate("theta2");
		forceExpTheta2 = tempTheta2.createCompiledExpression();
		counter += 1;
	}
	else calc[1] = 1;
 	if (inputExp.find("phi1") != -1) {
		calc[2] = 0;
		tempPhi1 = expression.differentiate("phi1");
		forceExpPhi1 = tempPhi1.createCompiledExpression();
		counter += 1;
	}
	else calc[2] = 1;
	if (inputExp.find("phi2") != -1) {
		calc[3] = 0;
		tempPhi2 = expression.differentiate("phi2");
		forceExpPhi2 = tempPhi2.createCompiledExpression();
		counter += 1;
	}
	else calc[3] = 1;
	if (counter == 0) {
		tempY101 = expression.differentiate("y10_1");
		tempY102 = expression.differentiate("y10_2");
		tempY11c1 = expression.differentiate("y11c_1");
		tempY11c2 = expression.differentiate("y11c_2");
		tempY11s1 = expression.differentiate("y11s_1");
		tempY11s2 = expression.differentiate("y11s_2");
		tempY201 = expression.differentiate("y20_1");
		tempY202 = expression.differentiate("y20_2");
		tempY21c1 = expression.differentiate("y21c_1");
		tempY21c2 = expression.differentiate("y21c_2");
		tempY21s1 = expression.differentiate("y21s_1");
		tempY21s2 = expression.differentiate("y21s_2");
		tempY22c1 = expression.differentiate("y22c_1");
		tempY22c2 = expression.differentiate("y22c_2");
		tempY22s1 = expression.differentiate("y22s_1");
		tempY22s2 = expression.differentiate("y22s_2");

		forceExpY101 = tempY101.createCompiledExpression();
		forceExpY102 = tempY102.createCompiledExpression();
		forceExpY11c1 = tempY11c1.createCompiledExpression();
		forceExpY11c2 = tempY11c2.createCompiledExpression();
		forceExpY11s1 = tempY11s1.createCompiledExpression();
		forceExpY11s2 = tempY11s2.createCompiledExpression();
		forceExpY201 = tempY201.createCompiledExpression();
		forceExpY202 = tempY202.createCompiledExpression();
		forceExpY21c1 = tempY21c1.createCompiledExpression();
		forceExpY21c2 = tempY21c2.createCompiledExpression();
		forceExpY21s1 = tempY21s1.createCompiledExpression();
		forceExpY21s2 = tempY21s2.createCompiledExpression();
		forceExpY22c1 = tempY22c1.createCompiledExpression();
		forceExpY22c2 = tempY22c2.createCompiledExpression();
		forceExpY22s1 = tempY22s1.createCompiledExpression();
		forceExpY22s2 = tempY22s2.createCompiledExpression();
	}
	YLM = counter;

	//differentiate
	tempR = expression.differentiate("r");

	//compile expressions
	EExpression = expression.createCompiledExpression();
	forceExpR = tempR.createCompiledExpression();

	for (int i = 0; i < numParameters; i ++)
		parameterNames.push_back(force.getPerParticleParameterName(i));
	for (int i = 0; i < force.getNumGlobalParameters(); i++) {
		globalParameterNames.push_back(force.getGlobalParameterName(i));
		globalParamValues[force.getGlobalParameterName(i)] = force.getGlobalParameterDefaultValue(i);
	}
	
	set<string> variables;
	variables.insert("r");
	variables.insert("theta1");
	variables.insert("theta2");
	variables.insert("phi1");
	variables.insert("phi2");
	variables.insert("y00");
	variables.insert("y10_1");
	variables.insert("y10_2");
	variables.insert("y11c_1");
	variables.insert("y11c_2");
	variables.insert("y11s_1");
	variables.insert("y11s_2");
	variables.insert("y20_1");
	variables.insert("y20_2");
	variables.insert("y21c_1");
	variables.insert("y21c_2");
	variables.insert("y21s_1");
	variables.insert("y21s_2");
	variables.insert("y22c_1");
	variables.insert("y22c_2");
	variables.insert("y22s_1");
	variables.insert("y22s_2");
	for (int i = 0; i < numParameters; i++) {
		variables.insert(parameterNames[i]+"1");
		variables.insert(parameterNames[i]+"2");
	}
	variables.insert(globalParameterNames.begin(),globalParameterNames.end());
	validateVariables(expression.getRootNode(),variables);
	for (auto& function: functions)
		delete function.second;

	isPeriodic = (nonbondedMethod == CutoffPeriodic);

/*
	for (int i = 0; i < force.getNumInteractionGroups(); i++) {
		set<int> set1,set2;
		force.getInteractionGroupParameters(i, set1, set2);
		interactionGroups.push_back(make_pair(set1, set2));
	}
*/
}
//EXECUTE
double ReferenceCalcCustomAnisotropicNonbondedForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
	vector<Vec3>& posData = extractPositions(context);
	vector<Vec3>& forceData = extractForces(context);
	Vec3* boxVectors = extractBoxVectors(context);

	ReferenceCustomAnisotropicNonbondedIxn ixn(EExpression,forceExpR,forceExpTheta1,forceExpTheta2,forceExpPhi1,forceExpPhi2,forceExpY101,forceExpY102,forceExpY11c1,forceExpY11c2,forceExpY11s1,forceExpY11s2,forceExpY201,forceExpY202,forceExpY21c1,forceExpY21c2,forceExpY21s1,forceExpY21s2,forceExpY22c1,forceExpY22c2,forceExpY22s1,forceExpY22s2,parameterNames);

	if (nonbondedMethod != NoCutoff) {
		computeNeighborListVoxelHash(*neighborList, numParticles,posData,exclusions,extractBoxVectors(context),isPeriodic, nonbondedCutoff, 0.0);
		ixn.setUseCutoff(nonbondedCutoff,*neighborList);
	}

	//check periodicity
	if(isPeriodic) {
		double minAllowedSize = 2*nonbondedCutoff;
		if (boxVectors[0][0] < minAllowedSize || boxVectors[1][1] < minAllowedSize || boxVectors[2][2] < minAllowedSize)
			throw OpenMMException("The periodic box size has decreased to less than twice the nonbonded cutoff.");
		ixn.setPeriodic(boxVectors);
	}
	bool globalParamsChanged = false;
	for (auto& name : globalParameterNames) {
		double value = context.getParameter(name);
		if (globalParamValues[name] != value)
			globalParamsChanged = true;
		globalParamValues[name] = value;
	}
	double energy = 0;
	//calculate interaction
	ixn.calculatePairIxn(numParticles, posData, axisTypes, atomXs, atomYs, atomZs, calc, particleParamArray, exclusions, 0, globalParamValues, forceData, includeEnergy ? &energy: NULL,YLM);
	return energy;
}

//COPY PARAMETERS TO CONTEXT
void ReferenceCalcCustomAnisotropicNonbondedForceKernel::copyParametersToContext(ContextImpl& context, const CustomAnisotropicNonbondedForce& force) {
	if (numParticles != force.getNumParticles())
		throw OpenMMException("updateParametersInContext: The number of particles has changed");

	int numParameters = force.getNumPerParticleParameters();
	vector<double> params;
	for (int i = 0; i < numParticles; ++i) {
		vector<double> parameters;
		int axisType, atomZ, atomX, atomY;
		force.getParticleParameters(i, parameters,axisType,atomX,atomY,atomZ);
		for (int j = 0 ; j < numParameters; j++)
			particleParamArray[i][j] = parameters[j];
		axisTypes[i] = axisType;
		atomZs[i] = atomZ;
		atomXs[i] = atomX;
		atomYs[i] = atomY;
	}

	if (forceCopy!= NULL) {
		*forceCopy = force;
	}
}

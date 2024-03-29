
#include "openmm/internal/AssertionUtilities.h"
#include "CustomAnisotropicNonbondedForce.h"
#include "openmm/serialization/XmlSerializer.h"
#include "openmm/Platform.h"
#include <iostream>
#include <sstream>

using namespace CustomAnisotropicNonbondedPlugin;
using namespace OpenMM;
using namespace std;

void testSerialization() {
    // Create a Force.

    CustomAnisotropicNonbondedForce force("5*sin(x)^2+y*z");
    force.setForceGroup(3);
    force.setNonbondedMethod(CustomAnisotropicNonbondedForce::CutoffPeriodic);
    force.setUseSwitchingFunction(true);
    force.setSwitchingDistance(2.0);
    force.setCutoffDistance(2.1);
    force.addGlobalParameter("x", 1.3);
    force.addGlobalParameter("y", 2.221);
    force.addPerParticleParameter("z");
    force.addEnergyParameterDerivative("y");
    vector<double> params(1);
    params[0] = 1.0;
    force.addParticle(params,0,1,0,0);
    params[0] = -3.3;
    force.addParticle(params,0,0,0,0);
    params[0] = 2.1;
    force.addParticle(params,0,0,0,0);
    force.addExclusion(0, 1);
    force.addExclusion(1, 2);
    vector<double> values(10);
    for (int i = 0; i < 10; i++)
        values[i] = sin((double) i);
    force.addFunction("f", values, 0.5, 1.5);
    std::set<int> set1, set2;
    set1.insert(0);
    set2.insert(1);
    set2.insert(2);
    force.addInteractionGroup(set1, set2);

    // Serialize and then deserialize it.

    stringstream buffer;
    XmlSerializer::serialize<CustomAnisotropicNonbondedForce>(&force, "Force", buffer);
    CustomAnisotropicNonbondedForce* copy = XmlSerializer::deserialize<CustomAnisotropicNonbondedForce>(buffer);

    // Compare the two forces to see if they are identical.

    CustomAnisotropicNonbondedForce& force2 = *copy;
    ASSERT_EQUAL(force.getForceGroup(), force2.getForceGroup());
    ASSERT_EQUAL(force.getEnergyFunction(), force2.getEnergyFunction());
    ASSERT_EQUAL(force.getNonbondedMethod(), force2.getNonbondedMethod());
    ASSERT_EQUAL(force.getCutoffDistance(), force2.getCutoffDistance());
    ASSERT_EQUAL(force.getSwitchingDistance(), force2.getSwitchingDistance());
    ASSERT_EQUAL(force.getUseSwitchingFunction(), force2.getUseSwitchingFunction());
    ASSERT_EQUAL(force.getNumPerParticleParameters(), force2.getNumPerParticleParameters());
    for (int i = 0; i < force.getNumPerParticleParameters(); i++)
        ASSERT_EQUAL(force.getPerParticleParameterName(i), force2.getPerParticleParameterName(i));
    ASSERT_EQUAL(force.getNumGlobalParameters(), force2.getNumGlobalParameters());
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        ASSERT_EQUAL(force.getGlobalParameterName(i), force2.getGlobalParameterName(i));
        ASSERT_EQUAL(force.getGlobalParameterDefaultValue(i), force2.getGlobalParameterDefaultValue(i));
    }
    ASSERT_EQUAL(force.getNumEnergyParameterDerivatives(), force2.getNumEnergyParameterDerivatives());
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++)
        ASSERT_EQUAL(force.getEnergyParameterDerivativeName(i), force2.getEnergyParameterDerivativeName(i));
    ASSERT_EQUAL(force.getNumParticles(), force2.getNumParticles());
    for (int i = 0; i < force.getNumParticles(); i++) {
        vector<double> params1, params2;
        int axisType1,axisType2,atomX1,atomX2,atomY1,atomY2,atomZ1,atomZ2;
        force.getParticleParameters(i, params1,axisType1,atomX1,atomY1,atomZ1);
        force2.getParticleParameters(i, params2,axisType2,atomX2,atomY2,atomZ2);
        ASSERT_EQUAL(params1.size(), params2.size());
        ASSERT_EQUAL(axisType1,axisType2);
        ASSERT_EQUAL(atomX1,atomX2);
        ASSERT_EQUAL(atomY1,atomY2);
        ASSERT_EQUAL(atomZ1,atomZ2);
        for (int j = 0; j < (int) params1.size(); j++) {
            ASSERT_EQUAL(params1[j], params2[j]);
		}
    }
    ASSERT_EQUAL(force.getNumExclusions(), force2.getNumExclusions());
    for (int i = 0; i < force.getNumExclusions(); i++) {
        int a1, a2, b1, b2;
        force.getExclusionParticles(i, a1, b1);
        force2.getExclusionParticles(i, a2, b2);
        ASSERT_EQUAL(a1, a2);
        ASSERT_EQUAL(b1, b2);
    }
    ASSERT_EQUAL(force.getNumTabulatedFunctions(), force2.getNumTabulatedFunctions());
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++) {
        double min1, min2, max1, max2;
        vector<double> val1, val2;
        dynamic_cast<Continuous1DFunction&>(force.getTabulatedFunction(i)).getFunctionParameters(val1, min1, max1);
        dynamic_cast<Continuous1DFunction&>(force2.getTabulatedFunction(i)).getFunctionParameters(val2, min2, max2);
        ASSERT_EQUAL(force.getTabulatedFunctionName(i), force2.getTabulatedFunctionName(i));
        ASSERT_EQUAL(min1, min2);
        ASSERT_EQUAL(max1, max2);
        ASSERT_EQUAL(val1.size(), val2.size());
        for (int j = 0; j < (int) val1.size(); j++)
            ASSERT_EQUAL(val1[j], val2[j]);
    }
    ASSERT_EQUAL(force.getNumInteractionGroups(), force2.getNumInteractionGroups());
    std::set<int> set1c, set2c;
    force2.getInteractionGroupParameters(0, set1c, set2c);
    ASSERT_EQUAL_CONTAINERS(set1, set1c);
    ASSERT_EQUAL_CONTAINERS(set2, set2c);
}

int main() {
    try {
        testSerialization();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}

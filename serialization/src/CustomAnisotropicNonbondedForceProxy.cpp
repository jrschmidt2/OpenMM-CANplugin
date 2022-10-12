
#include "CustomAnisotropicNonbondedForceProxy.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/Force.h"
#include "CustomAnisotropicNonbondedForce.h"
#include <sstream>

using namespace CustomAnisotropicNonbondedPlugin;
using namespace OpenMM;
using namespace std;

CustomAnisotropicNonbondedForceProxy::CustomAnisotropicNonbondedForceProxy() : SerializationProxy("CustomAnisotropicNonbondedForce") {
}

void CustomAnisotropicNonbondedForceProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 2);
    const CustomAnisotropicNonbondedForce& force = *reinterpret_cast<const CustomAnisotropicNonbondedForce*>(object);
    node.setIntProperty("forceGroup", force.getForceGroup());
    node.setStringProperty("energy", force.getEnergyFunction());
    node.setIntProperty("method", (int) force.getNonbondedMethod());
    node.setDoubleProperty("cutoff", force.getCutoffDistance());
    node.setBoolProperty("useSwitchingFunction", force.getUseSwitchingFunction());
    node.setDoubleProperty("switchingDistance", force.getSwitchingDistance());
    SerializationNode& perParticleParams = node.createChildNode("PerParticleParameters");
    for (int i = 0; i < force.getNumPerParticleParameters(); i++) {
        perParticleParams.createChildNode("Parameter").setStringProperty("name", force.getPerParticleParameterName(i));
    }
    SerializationNode& globalParams = node.createChildNode("GlobalParameters");
    for (int i = 0; i < force.getNumGlobalParameters(); i++) {
        globalParams.createChildNode("Parameter").setStringProperty("name", force.getGlobalParameterName(i)).setDoubleProperty("default", force.getGlobalParameterDefaultValue(i));
    }
    SerializationNode& energyDerivs = node.createChildNode("EnergyParameterDerivatives");
    for (int i = 0; i < force.getNumEnergyParameterDerivatives(); i++) {
        energyDerivs.createChildNode("Parameter").setStringProperty("name", force.getEnergyParameterDerivativeName(i));
    }
    SerializationNode& particles = node.createChildNode("Particles");
    for (int i = 0; i < force.getNumParticles(); i++) {
        int axisType, atomZ, atomX, atomY;
        vector<double> params;
        force.getParticleParameters(i, params,axisType,atomX,atomY,atomZ);
        SerializationNode& particle = particles.createChildNode("Particle");
        for (int j = 0; j < (int) params.size(); j++) {
            stringstream key;
            key << "param";
            key << j+1;
            particle.setDoubleProperty(key.str(), params[j]);
        }
        particle.setIntProperty("axisType", axisType);
        particle.setIntProperty("atomX", atomX);
        particle.setIntProperty("atomY", atomY);
        particle.setIntProperty("atomZ",atomZ);
    }
    SerializationNode& exclusions = node.createChildNode("Exclusions");
    for (int i = 0; i < force.getNumExclusions(); i++) {
        int particle1, particle2;
        force.getExclusionParticles(i, particle1, particle2);
        exclusions.createChildNode("Exclusion").setIntProperty("p1", particle1).setIntProperty("p2", particle2);
    }
    SerializationNode& functions = node.createChildNode("Functions");
    for (int i = 0; i < force.getNumTabulatedFunctions(); i++)
        functions.createChildNode("Function", &force.getTabulatedFunction(i)).setStringProperty("name", force.getTabulatedFunctionName(i));
        SerializationNode& interactionGroups = node.createChildNode("InteractionGroups");
            for (int i = 0; i < force.getNumInteractionGroups(); i++) {
                SerializationNode& interactionGroup = interactionGroups.createChildNode("InteractionGroup");
                std::set<int> set1;
                std::set<int> set2;
                force.getInteractionGroupParameters(i, set1, set2);
                SerializationNode& set1node = interactionGroup.createChildNode("Set1");
                for (int p : set1)
                    set1node.createChildNode("Particle").setIntProperty("index", p);
                SerializationNode& set2node = interactionGroup.createChildNode("Set2");
                for (int p : set2)
                    set2node.createChildNode("Particle").setIntProperty("index", p);
            }
}

void* CustomAnisotropicNonbondedForceProxy::deserialize(const SerializationNode& node) const {
    int version = node.getIntProperty("version");
    if (version < 1 || version > 2)
        throw OpenMMException("Unsupported version number");
    CustomAnisotropicNonbondedForce* force = NULL;
    try {
        CustomAnisotropicNonbondedForce* force = new CustomAnisotropicNonbondedForce(node.getStringProperty("energy"));
        force->setForceGroup(node.getIntProperty("forceGroup", 0));
        force->setNonbondedMethod((CustomAnisotropicNonbondedForce::NonbondedMethod) node.getIntProperty("method"));
        force->setCutoffDistance(node.getDoubleProperty("cutoff"));
        force->setUseSwitchingFunction(node.getBoolProperty("useSwitchingFunction", false));
        force->setSwitchingDistance(node.getDoubleProperty("switchingDistance", -1.0));
        const SerializationNode& perParticleParams = node.getChildNode("PerParticleParameters");
        for (auto& parameter : perParticleParams.getChildren())
            force->addPerParticleParameter(parameter.getStringProperty("name"));
        const SerializationNode& globalParams = node.getChildNode("GlobalParameters");
        for (auto& parameter : globalParams.getChildren())
            force->addGlobalParameter(parameter.getStringProperty("name"), parameter.getDoubleProperty("default"));
        if (version > 1) {
            const SerializationNode& energyDerivs = node.getChildNode("EnergyParameterDerivatives");
            for (auto& parameter : energyDerivs.getChildren())
                force->addEnergyParameterDerivative(parameter.getStringProperty("name"));
        }
        const SerializationNode& particles = node.getChildNode("Particles");
        vector<double> params(force->getNumPerParticleParameters());
        for (auto& particle : particles.getChildren()) {
            for (int j = 0; j < (int) params.size(); j++) {
                stringstream key;
                key << "param";
                key << j+1;
                params[j] = particle.getDoubleProperty(key.str());
            }
            int axisType = particle.getIntProperty("axisType");
            int atomX = particle.getIntProperty("atomX");
            int atomY = particle.getIntProperty("atomY");
            int atomZ = particle.getIntProperty("atomZ");
            force->addParticle(params,axisType,atomX,atomY,atomZ);
        }
        const SerializationNode& exclusions = node.getChildNode("Exclusions");
        for (auto& exclusion : exclusions.getChildren())
            force->addExclusion(exclusion.getIntProperty("p1"), exclusion.getIntProperty("p2"));
        const SerializationNode& functions = node.getChildNode("Functions");
        for (auto& function : functions.getChildren()) {
            if (function.hasProperty("type")) {
                force->addTabulatedFunction(function.getStringProperty("name"), function.decodeObject<TabulatedFunction>());
            }
            else {
                // This is an old file created before TabulatedFunction existed.

                const SerializationNode& valuesNode = function.getChildNode("Values");
                vector<double> values;
                for (auto& child : valuesNode.getChildren())
                    values.push_back(child.getDoubleProperty("v"));
                force->addTabulatedFunction(function.getStringProperty("name"), new Continuous1DFunction(values, function.getDoubleProperty("min"), function.getDoubleProperty("max")));
            }
        }
     bool hasInteractionGroups = false; // Older files will be missing this block.
        for (auto& child : node.getChildren())
            if (child.getName() == "InteractionGroups")
                hasInteractionGroups = true;
        if (hasInteractionGroups) {
            const SerializationNode& interactionGroups = node.getChildNode("InteractionGroups");
            for (auto& interactionGroup : interactionGroups.getChildren()) {
                // Get set 1.
                const SerializationNode& set1node = interactionGroup.getChildNode("Set1");
                std::set<int> set1;
                for (auto& child : set1node.getChildren())
                set1.insert(child.getIntProperty("index"));
                // Get set 2.
                const SerializationNode& set2node = interactionGroup.getChildNode("Set2");
                std::set<int> set2;
                for (auto& child : set2node.getChildren())
                set2.insert(child.getIntProperty("index"));
                force->addInteractionGroup(set1, set2);
                }
        }
        return force;
    }
    catch (...) {
        if (force != NULL)
            delete force;
        throw;
    }
}

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
%module canplugin

%import(module="simtk.openmm") "swig/OpenMMSwigHeaders.i"
%include "swig/typemaps.i"

%include "std_vector.i"
%include "std_string.i"
%include "std_pair.i"
%include "std_set.i"
namespace std {
  %template(pairii) pair<int,int>;
  %template(vectord) vector<double>;
  %template(vectori) vector<int>;
  %template(vectorpairii) vector< pair<int,int> >;
  %template(seti) set<int>;
};

%include "typemaps.i"

%{
#include "CustomAnisotropicNonbondedForce.h"
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "OpenMMDrude.h"
#include "openmm/RPMDIntegrator.h"
#include "openmm/RPMDMonteCarloBarostat.h"
%}

%pythoncode %{
import simtk.openmm as mm
import simtk.unit as unit
%}

/*
 * Add units
*/
%pythonappend CustomAnisotropicNonbondedPlugin::CustomAnisotropicNonbondedForce::getCutoffDistance() const %{
    val = unit.Quantity(val, unit.nanometer)
%}
%pythonappend CustomAnisotropicNonbondedPlugin::CustomAnisotropicNonbondedForce::getSwitchingDistance() const %{
    val = unit.Quantity(val, unit.nanometer)
%}
/*
* Extra things for tabulated functions
*/
%pythonprepend CustomAnisotropicNonbondedPlugin::CustomAnisotropicNonbondedForce::addTabulatedFunction(const std::string &name, OpenMM::TabulatedFunction *function) %{
    if not function.thisown:
        s = ("the %s object does not own its corresponding OpenMM object"
             % self.__class__.__name__)
        raise Exception(s)
%}
%pythonappend CustomAnisotropicNonbondedPlugin::CustomAnisotropicNonbondedForce::addTabulatedFunction(const std::string & name, OpenMM::TabulatedFunction * function) %{
   function.thisown=0
%}

/*
* Convert exceptions
*/
%exception {
    try {
        $action
    } catch (std::exception &e) {
        PyErr_SetString(PyExc_Exception, const_cast<char*>(e.what()));
        return NULL;
    }
}


/*
*  Classify force
*/

namespace CustomAnisotropicNonbondedPlugin {
class CustomAnisotropicNonbondedForce : public OpenMM::Force {
public:
    /*
    * Methods in CAN
    */
    CustomAnisotropicNonbondedForce(const std::string& energy);

    enum NonbondedMethod{
        NoCutoff = 0,
        CutoffNonPeriodic = 1,
        CutoffPeriodic = 2
    };
    enum AxisTypes{
        ZThenX = 0,
        Bisector = 1,
        ZBisect = 2,
        ThreeFold = 3,
        ZOnly = 4,
        NoAxisType = 5,
        LastAxisTypeIndex = 6
    };
    int getNumParticles() const;
    int getNumExclusions() const;
    int getNumPerParticleParameters() const;
    int getNumGlobalParameters() const;
    int getNumTabulatedFunctions() const;
    int getNumFunctions() const;
    int getNumInteractionGroups() const;
    int getNumEnergyParameterDerivatives() const;

    const std::string getEnergyFunction() const;
    void setEnergyFunction(const std::string& energy);
   
    CustomAnisotropicNonbondedForce::NonbondedMethod getNonbondedMethod() const;
    void setNonbondedMethod(CustomAnisotropicNonbondedForce::NonbondedMethod method);

    double getCutoffDistance() const;
    void setCutoffDistance(double distance);

    bool getUseSwitchingFunction() const;
    void setUseSwitchingFunction(bool use);
    double getSwitchingDistance() const;
    void setSwitchingDistance(double distance);
   
    int checkAngles() const;

    int addPerParticleParameter(const std::string& name);
    const std::string& getPerParticleParameterName(int index) const;
    void setPerParticleParameterName(int index, const std::string& name);
    int addGlobalParameter(const std::string& name, double defaultValue);
    const std::string& getGlobalParameterName(int index) const;
    void setGlobalParameterName(int index, const std::string& name);
    double getGlobalParameterDefaultValue(int index) const;
    void setGlobalParameterDefaultValue(int index, double defaultValue);

    void addEnergyParameterDerivative(const std::string& name);
    const std::string& getEnergyParameterDerivativeName(int index) const;

    int addParticle(const std::vector<double>& parameters=std::vector<double>(), int axisType=5, int atomX=-1, int atomY=-1, int atomZ=-1);

    %apply std::vector<double>& OUTPUT {std::vector<double>& parameters};
    %apply int& OUTPUT {int& axisType};
    %apply int& OUTPUT {int& atomX};
    %apply int& OUTPUT {int& atomY};
    %apply int& OUTPUT {int& atomZ};
    void getParticleParameters(int index, std::vector<double>& parameters, int& axisType, int& atomX, int& atomY, int& atomZ) const;
    %clear std::vector<double>& parameters;
    %clear int& axisType;
    %clear int& atomX;
    %clear int& atomY;
    %clear int& atomZ;

    void setParticleParameters(int index, const std::vector<double>& parameters, int axisType, int atomX, int atomY, int atomZ);

    int addExclusion(int particle1, int particle2);

    %apply int& OUTPUT {int& particle1};
    %apply int& OUTPUT {int& particle2};
    void getExclusionParticles(int index, int& particle1, int& particle2) const;
    %clear int& particle1;
    %clear int& particle2;

    void setExclusionParticles(int index, int particle1, int particle2);
    void createExclusionsFromBonds(const std::vector<std::pair<int,int> >& bonds, int bondCutoff); 
    int addTabulatedFunction(const std::string& name, OpenMM::TabulatedFunction* function);
    const OpenMM::TabulatedFunction& getTabulatedFunction(int index) const;
    OpenMM::TabulatedFunction& getTabulatedFunction(int index);
    const std::string& getTabulatedFunctionName(int index) const;

    int addFunction(const std::string& name, const std::vector<double>& values, double min, double max);
   
    %apply std::string& OUTPUT {std::string& name}
    %apply std::vector<double>& OUTPUT {std::vector<double>& values}
    %apply double& OUTPUT {double& min}
    %apply double& OUTPUT {double& max}
    void getFunctionParameters(int index, std::string& name, std::vector<double>& values, double& min, double& max) const;
    %clear std::string& name;
    %clear std::vector<double>& values;
    %clear double& min;
    %clear double& max;

    void setFunctionParameters(int index, const std::string& name, const std::vector<double>& values, double min, double max);
    int addInteractionGroup(const std::set<int>& set1, const std::set<int>& set2);
   
    %apply std::set<int>& OUTPUT {std::set<int>& set1};
    %apply std::set<int>& OUTPUT {std::set<int>& set2};
    void getInteractionGroupParameters(int index, std::set<int>& set1, std::set<int>& set2) const;
    %clear std::set<int>& set1;
    %clear std::set<int>& set2;

    void setInteractionGroupParameters(int index, const std::set<int>& set1, const std::set<int>& set2);
    %apply OpenMM::Context & OUTPUT { Context & context }; 
    void updateParametersInContext(OpenMM::Context& context);
    %clear OpenMM::Context & context;
    bool usesPeriodicBoundaryConditions() const;
    /*
    * Cast Force to CustomAnisotropicNonbondedForce
    */
    %extend {
        static CustomAnisotropicNonbondedPlugin::CustomAnisotropicNonbondedForce& cast(OpenMM::Force& force) {
            return dynamic_cast<CustomAnisotropicNonbondedPlugin::CustomAnisotropicNonbondedForce&>(force);
        }
        static bool isinstance(OpenMM::Force& force) {
            return (dynamic_cast<CustomAnisotropicNonbondedPlugin::CustomAnisotropicNonbondedForce*>(&force) != NULL);
        }
    }


};
}

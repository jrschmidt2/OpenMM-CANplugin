#ifndef __ReferenceCustomAnisotropicNonbondedIxn_H__
#define __ReferenceCustomAnisotropicNonbondedIxn_H__
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

#include "openmm/Vec3.h"
#include "openmm/reference/ReferencePairIxn.h"
#include "openmm/reference/ReferenceNeighborList.h"
#include "openmm/internal/CompiledExpressionSet.h"

#include <map>
#include <set>
#include <utility>
#include <vector>

namespace CustomAnisotropicNonbondedPlugin {

class ReferenceCustomAnisotropicNonbondedIxn {

    private:

        bool cutoff;
        bool useSwitch;
        bool periodic;
        const OpenMM::NeighborList* neighborList;
        OpenMM::Vec3 periodicBoxVectors[3];
        double cutoffDistance, switchingDistance;
        std::vector<std::string> paramNames;
//      std::vector<std::pair<std::set<int>, std::set<int> > > interactionGroups;
        Lepton::CompiledExpression EExp,fExpR,fExpTheta1,fExpTheta2,fExpPhi1,fExpPhi2;
        OpenMM::CompiledExpressionSet expressionSet;
        std::vector<int> particleParamIndex;
        int rIndex, theta1Index,theta2Index,phi1Index,phi2Index;

        /**---------------------------------------------------------------------------------------
         Calculate custom pair ixn between two atoms

         @param atom1            the index of the first atom
         @param atom2            the index of the second atom
         @param atomCoordinates  atom coordinates
         @param forces           force array (forces added)
         @param axisTypes        numeric identifier of local axis type
         @param atomXs           indices of atoms along x local coordinate
         @param atomYs           indices of atoms along y local coordinate
         @param atomZs           indices of atoms along z local coordinate
         @param calc             array indicating which angles to compute
         @param kxvector[1,2]    local x coordinate axis for atom[1,2]
         @param kyvector[1,2]    local x coordinate axis for atom[1,2]
         @param kzvector[1,2]    local z coordinate axis for atom[1,2]
         @param variables        parameters for given function
         @param totalEnergy      total energy
         --------------------------------------------------------------------------------------- */
        void calculateOneIxn(int atom1, int atom2, std::vector<OpenMM::Vec3>& atomCoordinates, std::vector<OpenMM::Vec3>& forces, int* axisTypes,int* atomXs, int* atomYs,int* atomZs,int* calc,
                           OpenMM::Vec3& kxvector1, OpenMM::Vec3& kxvector2, OpenMM::Vec3& kyvector1,OpenMM::Vec3& kyvector2,OpenMM::Vec3& kzvector1, OpenMM::Vec3& kzvector2,
                           double* totalEnergy);

        /*----------------------------------------------------------------------------------------
        Calculate angle between two vectors without fixed orientation [0,PI]
        @param vec1		first vector in angle
        @param vec2		second vector in angle
        -----------------------------------------------------------------------------------------*/
        static double computeAzim(const OpenMM::Vec3& vec1, const OpenMM::Vec3& vec2);

        /*----------------------------------------------------------------------------------------
        Calculates angle between two vectors in a plane with fixed orientation [0,2PI]
        @param rij_in		vector to be projected onto plane before angle calculation
        @param pvec		vector fixed in plane prior to angle calculation
        @param ovec		orthonormal vector to plane
        -----------------------------------------------------------------------------------------*/
        static double computePolar(const OpenMM::Vec3& rij_in, const OpenMM::Vec3& pvec, const OpenMM::Vec3& ovec);

    public:

        /**---------------------------------------------------------------------------------------
         Constructor
         --------------------------------------------------------------------------------------- */
        ReferenceCustomAnisotropicNonbondedIxn(const Lepton::CompiledExpression& EExpression, const Lepton::CompiledExpression& forceExpR,
                                               const Lepton::CompiledExpression& forceExpTheta1, const Lepton::CompiledExpression& forceExpTheta2,
                                               const Lepton::CompiledExpression& forceExpPhi1, const Lepton::CompiledExpression& forceExpPhi2,
                                               const std::vector<std::string>& parameterNames);

        /**---------------------------------------------------------------------------------------
         Destructor
         --------------------------------------------------------------------------------------- */
        ~ReferenceCustomAnisotropicNonbondedIxn();

        /**--------------------------------------------------------------------------------------
 	 Normalizes vector
	 @param vectorToNormalize	vector to normalize
	 */
        double normalizeVec3(OpenMM::Vec3& vectorToNormalize) const;

        /**---------------------------------------------------------------------------------------
         Set the force to use a cutoff.
         @param distance            the cutoff distance
         @param neighbors           the neighbor list to use
         --------------------------------------------------------------------------------------- */
        void setUseCutoff(double distance, const OpenMM::NeighborList& neighbors);

        /**---------------------------------------------------------------------------------------
         Restrict the force to a list of interaction groups.
         @param distance            the cutoff distance
         @param neighbors           the neighbor list to use
         --------------------------------------------------------------------------------------- */
        void setInteractionGroups(const std::vector<std::pair<std::set<int>, std::set<int> > >& groups);

        /**---------------------------------------------------------------------------------------
         Set the force to use a switching function.
         @param distance            the switching distance
         --------------------------------------------------------------------------------------- */
        void setUseSwitchingFunction(double distance);

        /**---------------------------------------------------------------------------------------
         Set the force to use periodic boundary conditions.  This requires that a cutoff has
         already been set, and the smallest side of the periodic box is at least twice the cutoff
         distance.
         @param vectors    the vectors defining the periodic box
         --------------------------------------------------------------------------------------- */
        void setPeriodic(OpenMM::Vec3* vectors);

	/*--------------------------------------------------------------------------------------
         Calculates local coordinate axes for given atom
         @param particleI    particle for which axes should be calculated
         @param particleX    particle in the kx direction with particleI as the origin
         @param particleY    particle in the ky direction with particleI as the origin
         @param particleZ    particle in the kz direction with particleI as the origin
         @param axisType    define type of axis to calculate
         --------------------------------------------------------------------------------------- */
//        double * accessAxisParameter(std::vector<OpenMM::Vec3>& atomCoordinates, int particleI, int particleX, int particleY, int particleZ, int axisType) const;

        /**---------------------------------------------------------------------------------------
         Calculate custom pair ixn
         @param numberOfAtoms    number of atoms
         @param atomCoordinates  atom coordinates
         @param axisTypes        integer vector of axis types needed for all atoms in accessAxisParameter
         @param atomXs           integer vector of X indices needed for all atoms in accessAxisParameter
         @param atomYs           integer vector of Y indices needed for all atoms in accessAxisParameter
         @param atomZs           integer vector of Z indices needed for all atoms in accessAxisParameter
         @param calc             array indicating which angles to compute
         @param atomParameters   atom parameters (charges, c6, c12, ...)     atomParameters[atomIndex][paramterIndex]
         @param exclusions       atom exclusion indices
                                 exclusions[atomIndex] contains the list of exclusions for that atom
         @param fixedParameters  non atom parameters (not currently used)
         @param globalParameters the values of global parameters
         @param forces           force array (forces added)
         @param energyByAtom     atom energy
         @param totalEnergy      total energy
         --------------------------------------------------------------------------------------- */

        void calculatePairIxn(int numberOfAtoms, std::vector<OpenMM::Vec3>& atomCoordinates, int* axisTypes, 
                              int* atomXs, int* atomYs, int* atomZs, int* calc,
                              double** atomParameters, std::vector<std::set<int> >& exclusions, double* fixedParameters, 
                              const std::map<std::string, double>& globalParameters, std::vector<OpenMM::Vec3>& forces, 
                              double*  totalEnergy);
// ---------------------------------------------------------------------------------------

};


} // namespace CustomAnisotropicNonbondedPlugin

#endif // __ReferenceCustomAnisotropicNonbondedIxn_H__

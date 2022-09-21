Custom Anisotropic Nonbonded Force Plugin
Plugin for OpenMM (v7.1)  
Tesia Janicki  
cr. 5.12.2017 ed. 9.21.2022 
components borrowed from existing OpenMM forces to ensure smooth integration
---
**Overview** <br />
This repository houses a plugin to incorporate user-defined force expressions which include atomic anisotropy. Atomic anisotropy is defined with respect to the local coordinate system of each atom within a pair interaction (e.g. between atoms i and j, shown in Figure 1). Where isotropic forces are defined as a function of only interatomic distance, forces incorporating atomic anisotropy require positional definitions of all atoms. For example, in addition to interatomic distance, atom j must be defined by a polar phi1 and azimuthal theta1 of atom i's coordinate system, and atom i must be defined by polar phi2 and azimuthal theta2 of atom j's coordinate system. Details on how each atom's local coordinate system is defined by the user as a function of local bonding is described in further detail below. <br />
Figure 1: Depiction of local coordinate systems for atomic anisotropy definitions.
<img src="https://github.com/tesiadj16/OpenMM-CANplugin/blob/master/axOrient.PNG" width = "500"> 
---
**Plugin Directory Framework** <br />
I. openmm_can  	<br />									
 A. openmmcustomanisotropicnonbonded  	<br />
   1. openmmapi									
   2. platforms									
       * cuda
       * reference				
   3. python									
   4. serialization									
---
**Necessary Software for Installation**
 1. OpenMM 7.1 (preferred); this plugin has not been tested on earlier versions	
 2. CMAKE 3.1 (minimum) for current cmake files					
 3. C++ 11
 4. SWIG 3.0.8 (minimum)									
---
**Installation Notes**
 1. While in 'openmm_can', run cmake on the plugin subdirectory (i.e. 'ccmake openmmcustomanisotropicnonbonded')															
    Press 'c' for compile										
    Make any necessary modifications to the given fields											
    Compile ('c') if any changes were made 													
    Press 'g' to generate															
 2. While in 'openmm_can', enter 'make' to compile												
 3. Enter 'make install' to install the plugin													
 4. Enter 'make test' to run tests on the installed plugin											
 5. If installing python wrappers:
	* Enter 'make PythonInstall'
	* Identify forcefield.py, the python script used to interpret the user's force field implementation. This file exists in wrappers/python/simtk/openmm/app within the openmm source code.
	* Permit use of plugin libraries by adding the following line to the header:
```
import canplugin as can
```
	* Append the contents of appendPython.txt to forcefield.py so the python layer can interpret user input.
	* Re-make and install OpenMM with python wrappers.
---
**Execution Notes: General**
 1. The plugin assumes inclusion of anisotropy in force expression. These computations are not 'turned off' in the absence of anisotropic components.  Absence of angles or spherical harmonics will produce a warning which directs the user to CustomNonbondedForce for more efficient custom isotropic nonbonded expressions.
 2. Forces may be defined in terms of the pair distance (r) and the anisotropic components (listed below).  As such, these variables should not be redefined as custom variables in the force expression (e.g. theta1 only ever means the polar angle in the local coordinate system of atom1).
 3. Anisotropy may be incorporated either as a function of angles.
	* Angles and spherical harmonic expressions *may not* be used in the same force expression to prevent confusion. An error will be returned in this case.
	* Angles are defined in terms of polar (phi1, phi2) and azimuthal (theta1, theta2) components for the local coordinate systems of atom1 and atom2 in a given pair interaction.
	
 4. A given atom i is defined in terms of its index, and the indices of atoms defining the local symmetry about atom i. This document refers to the (at most) three atoms definig this symmetry as atomX, atomY, atomZ, as described in the following table: <br />
<br />
Table 2: Axis Definitions <br />
<img src="https://github.com/tesiadj16/OpenMM-CANplugin/blob/master/axDef.PNG" width = "500"> 
<br />
** Integers in parentheses should be input for axisType when using c++ directly; in python, a script will interpret this identification automatically.
<br />
*** Unused Atom(X/Y/Z) are indicated as (-1) when using c++ directly; in python, simply omit the unused value.
 5. Definitions of local coordinate systems do not change during the simulation, but the numerical values are updated for each timestep.  For example, if the jth atom is defined as atomZ for the ith atom for a ZOnly axis type, rij will always define the axis of most symmetry.

**Execution Notes: C++**
 1. The plugin may be accessed in c++ submission scripts with similar syntax to other OpenMM forces.
 2. When compiling, be sure to include the relevant library: -lCustomAnisotropicNonbondedPlugin
 3. Define the force as CustomAnisotropicNonbondedForce("energy expression")
 4. Global and local (per atom) parameters must be explicitly defined, as in other OpenMM custom forces.
 5. Particles are added to the system using the addParticle function, with any unused axis identifier set to (-1). The axisType must be defined explicitly, per integer identifier in Table 2. Per-particle parameters are defined as a vector of their values.
```
 addParticle(peratomparametervector,axisType,atomX,atomY,atomZ)
```
**Execution Notes: Python**
 1. If installed properly, the plugin may be accessed in a similar way as other OpenMM forces.
 2. The force plugin should be imported as:
```
from canplugin import CustomAnisotropicNonbondedForce
```
 3. After defining the force expression, each atom in a given molecule should be identified in terms of its atom type and axis parameters, as described in **Execution Notes: General**.  Axis parameters should be input in Z, X, Y order; if the given axis type does not require all identifiers, they should be omitted from the descriptor (as for CO2 in the xml example below). Any other per-particle parameters should follow the axis parameters.
 4. A sample xml input for CO2 for a fictitious force is given below:
```
<ForceField>
 <AtomTypes>
  <Type name="1" class="C_CO2" element="C" mass="12.01060"/>
  <Type name="2" class="O1_CO2" element="O" mass="15.99943"/>
  <Type name="3" class="O2_CO2" element="O" mass="15.99943"/>
 </AtomTypes>
 <Residues>
  <Residue name="co2">
   <Atom name="C" type="1"/>
   <Atom name="O1" type="2"/>
   <Atom name="O2" type="3"/>
   <Bond from="0" to="1"/>
   <Bond from="0" to="2"/>
  </Residue>
 </Residues>
 <CustomAnisotropicNonbondedForce bondCutoff="5.0"
	energy = "A*r*(sin(theta1)+sin(theta2)) + B1*B2" >
	<GlobalParameter name="A" defaultValue="1.0"/>
	<PerParticleParameter name="B" />
	<Atom type="1" AtomZ= "2" B= "2.0" />
	<Atom type="2" AtomZ= "1" B= "1.0" />
	<Atom type="3" AtomZ= "1" B= "1.5" />
 </CustomAnisotropicNonbondedForce>
</ForceField>

```

---

from __future__ import print_function
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from sys import stdout, exit
import numpy as np
import canplugin as can
import mdtraj

# Specify thermodynamic conditions
temperature = 298*unit.kelvin
pressure = 1*unit.atmospheres
# Specify cutoffs
two_body_cutoff = 1.2*unit.nanometers
electrostatic_two_body_cutoff = 1.2*unit.nanometers
three_body_cutoff = 0.7*unit.nanometers


# Load xml definition
app.topology.Topology.loadBondDefinitions('residues.xml')
print('Loaded bond definitions ')
pdb = app.PDBFile('liquid_benzene.pdb')
print('Loaded PDB file')
forcefield = app.ForceField('benzene_2b_3b_lrc_noljpme_07162021_no15.xml')
print('Loaded force field')

##System
system = forcefield.createSystem(model_topology, 
                nonbondedMethod=app.NoCutoff,
                nonbondedCutoff=electrostatic_two_body_cutoff,
                constraints=None,
                rigidWater=True,
                polarization='mutual',
                ewaldErrorTolerance=0.0005)
print('Created system')

# Set distance cutoffs, constraints, and other force-specific options
for force in system.getForces():
    if isinstance(force, mm.CustomHbondForce):
        force.setNonbondedMethod(mm.CustomHbondForce.CutoffPeriodic)
        force.setCutoffDistance(two_body_cutoff)
    elif can.CustomAnisotropicNonbondedForce.isinstance(force):
        force1 = can.CustomAnisotropicNonbondedForce.cast(force)
        force1.setNonbondedMethod(can.CustomAnisotropicNonbondedForce.CutoffPeriodic)
        force1.setCutoffDistance(two_body_cutoff)
        force1.setUseSwitchingFunction(False)
    elif isinstance(force, mm.CustomNonbondedForce):
        force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        force.setCutoffDistance(two_body_cutoff)
        force.setUseLongRangeCorrection(True)
        force.setUseSwitchingFunction(False)
    elif isinstance(force, mm.AmoebaMultipoleForce):
        force.setNonbondedMethod(mm.AmoebaMultipoleForce.PME)
        force.setCutoffDistance(electrostatic_two_body_cutoff)
    elif isinstance(force, mm.CustomManyParticleForce):
        force.setNonbondedMethod(mm.CustomManyParticleForce.CutoffPeriodic)
        force.setCutoffDistance(three_body_cutoff)
    elif isinstance(force, mm.HarmonicBondForce):
        force.setUsesPeriodicBoundaryConditions(True)
    elif isinstance(force, mm.HarmonicAngleForce):
        force.setUsesPeriodicBoundaryConditions(True)
    elif isinstance(force, mm.CustomTorsionForce):
        force.setUsesPeriodicBoundaryConditions(True)
    else:
        pass

##Integrator
integrator = mm.LangevinIntegrator(temperature, 2.0/unit.picoseconds, 1.0*unit.femtoseconds)

##Barostat
barostat = mm.MonteCarloBarostat(pressure, temperature, 25)
system.addForce(barostat)

##Platform & Simulation
platform = mm.Platform.getPlatformByName('CUDA')
properties={}
properties["CudaPrecision"] = "mixed"
simulation = app.Simulation(model_topology, system, integrator, platform,properties)
simulation.context.setPositions(model_positions)

##Initial minimization
state = simulation.context.getState(getEnergy=True)
print(state.getPotentialEnergy())

simulation.context.setVelocitiesToTemperature(temperature)
simulation.minimizeEnergy()
simulation.step(2000000)

##Running & Recording
simulation.reporters.append(app.StateDataReporter('output10.out', 1000, step=True,
        potentialEnergy=True, kineticEnergy=True, temperature=True,volume=True,
        density=True, speed=True, time=True, separator='\t'))
simulation.reporters.append(app.PDBReporter('output10.pdb',1000))

print('Running Production...')
for i in range(10000):
    simulation.step(1000)
print('Done!')

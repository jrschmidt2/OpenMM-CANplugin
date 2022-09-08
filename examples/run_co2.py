from __future__ import print_function
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
from sys import stdout
from canplugin import CustomAnisotropicNonbondedForce
from time import time

#Need a force reporter, which is not internal to OpenMM

class ForceReporter(object):
	def __init__(self, file, reportInterval):
		self._out = open(file,'w')
		self._out.write('Force Components (x,y,z) by Atom in kJ/mol \n')
		self._reportInterval = reportInterval
	def __del__(self):
		self._out.close()
	def describeNextReport(self, simulation):
		steps = self._reportInterval - simulation.currentStep%self._reportInterval
		return (steps, False, False, True, False)
	def report(self, simulation, state):
		forces = state.getForces().value_in_unit(unit.kilojoule/unit.mole/unit.nanometer)
		self._out.write('Step %8.0f \n' % simulation.currentStep)
		for f in forces:
			self._out.write('%20.5f %20.5f %20.5f \n' % (f[0], f[1], f[2]))
####

#constants
tempset = 298*unit.kelvin
presset=1*unit.atmospheres

#set model
app.topology.Topology.loadBondDefinitions('residues.xml')
pdb=app.PDBFile('co2_nox.pdb')
model = app.modeller.Modeller(pdb.topology, pdb.positions)
model_topology = model.getTopology()
model_positions = model.getPositions()

#set force field and system
#forcefield = app.ForceField('CO2_example.xml')
forcefield = app.ForceField('CO2_example_YLM.xml')
system = forcefield.createSystem(model_topology, nonbondedMethod=app.CutoffPeriodic,nonbondedCutoff=10.0*unit.angstrom,constraints=app.AllBonds)

#simulation details
integrator = mm.VerletIntegrator(0.01*unit.femtoseconds)
platform = mm.Platform.getPlatformByName('Reference')
simulation = app.Simulation(model_topology,system,integrator,platform)
simulation.context.setPositions(model_positions)

#setup reporters
simulation.reporters.append(app.PDBReporter('trajectory.pdb',1))
simulation.reporters.append(ForceReporter('forces.dat',1))
simulation.reporters.append(app.StateDataReporter('output.dat',1,step=True,potentialEnergy=True,kineticEnergy=True,totalEnergy=True,density=True,volume=True,temperature=True,progress=True,remainingTime=True,totalSteps=1,separator=','))

#run simulation
print('Running Production...')
start = time()
simulation.step(10)
end = time()
print((end-start)/10)
print('Done!')


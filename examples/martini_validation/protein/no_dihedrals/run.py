import msprot
import sqlite3
from   openmm.app import *
from   openmm import *
from   openmm.unit import *
from   sys import stdout

martini = msprot.ReadMartini()
msprot.MartinizeDMS('protein_cg.pdb', martini=martini, output='input.martini.dms')
msprot.dumpsql('input.martini.dms')

system, dms = msprot.DMS2openmm('input.martini.dms').make()
#forces = { force.__class__.__name__ : force for force in system.getForces() }

### PRESSURE
#system.addForce(MonteCarloMembraneBarostat(1*bar, 0*bar*nanometer, 310*kelvin, 
#    MonteCarloMembraneBarostat.XYIsotropic, MonteCarloMembraneBarostat.ZFree))
system.addForce(MonteCarloBarostat(1*bar, 310*kelvin))

### RUN SIMS
integrator = LangevinMiddleIntegrator(310*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(dms.topology, system, integrator)
simulation.context.setPositions(dms.positions)
print(simulation.context.getState(getEnergy=True).getPotentialEnergy())
#simulation.minimizeEnergy()
#simulation.reporters.append(DCDReporter('input.dcd', 1000))
#simulation.reporters.append(StateDataReporter('input.csv', 1000, step=True, potentialEnergy=True, temperature=True))
#simulation.step(1000000)

msprot.gmx_energy(c='protein_cg.pdb', p='topol.top')



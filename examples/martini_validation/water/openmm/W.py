import msprot
import sqlite3
from   openmm.app import *
from   openmm import *
from   openmm.unit import *

martini = msprot.ReadMartini()
msprot.MartinizeDMS('W.pdb', martini=martini)
msprot.dumpsql('W.martini.dms')

system, dms = msprot.DMS2openmm('W.martini.dms').make()
forces = { force.__class__.__name__ : force for force in system.getForces() }

### PRESSURE
system.addForce(MonteCarloBarostat(1*bar, 310*kelvin, 250))

### RUN SIMS
integrator = LangevinMiddleIntegrator(310*kelvin, 1/picosecond, 0.02*picoseconds)
simulation = Simulation(dms.topology, system, integrator)
simulation.context.setPositions(dms.positions)
simulation.minimizeEnergy()
simulation.reporters.append(DCDReporter('W.dcd', 1000))
simulation.step(1000000)


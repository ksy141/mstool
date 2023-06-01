import msprot
import sqlite3
from   openmm.app import *
from   openmm import *
from   openmm.unit import *

martini = msprot.ReadMartini()
msprot.MartinizeDMS('../input.pdb', martini=martini, output='input.martini.dms')
msprot.dumpsql('input.martini.dms')

system, dms = msprot.DMS2openmm('input.martini.dms').make()
forces = { force.__class__.__name__ : force for force in system.getForces() }

### PRESSURE
system.addForce(MonteCarloMembraneBarostat(1*bar, 0*bar*nanometer, 310*kelvin, 
    MonteCarloMembraneBarostat.XYIsotropic, MonteCarloMembraneBarostat.ZFree))
#system.addForce(MonteCarloBarostat(1*bar, 310*kelvin))

### RUN SIMS
integrator = LangevinMiddleIntegrator(310*kelvin, 1/picosecond, 0.02*picoseconds)
simulation = Simulation(dms.topology, system, integrator)
simulation.context.setPositions(dms.positions)
simulation.minimizeEnergy()
simulation.reporters.append(DCDReporter('input.dcd', 1000))
simulation.step(500000)


from openmm.unit import *

def getEnergy(simulation):
    return simulation.context.getState(getEnergy=True).getPotentialEnergy()._value

def Pressure(system, P=1.0, T=310.0):
    system.addForce(MonteCarloBarostat(P*bar, T*kelvin))
    return system

def MembranePressure(system, P=1.0, T=310.0, r=0.0):
    system.addForce(MonteCarloMembraneBarostat(P*bar, r*bar*nanometer, T*kelvin,
            MonteCarloMembraneBarostat.XYIsotropic, MonteCarloMembraneBarostat.ZFree))
    return system


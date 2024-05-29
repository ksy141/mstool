import numpy as np

def SDMinimizer(simulation, h0=0.01, F_max_tol = 1e3, nsteps=1000):
    '''
    https://manual.gromacs.org/current/reference-manual/algorithms/energy-minimization.html
    openMM and gromacs use kJ/mol and nm;
    h0: initial maximum displacement
    '''

    for i in range(1, nsteps+1):
        state0 = simulation.context.getState(getEnergy=True, getForces=True, getPositions=True)
        E0     = state0.getPotentialEnergy()._value
        F0     = state0.getForces(asNumpy=True)._value
        p0     = state0.getPositions(asNumpy=True)._value
        F0max  = np.linalg.norm(F0, axis=1).max()
        
        if F0max < F_max_tol:
            print(f'Fmax: {F0max:.1e} < F_max_tol: {F_max_tol:.1e}')
            break
        
        p1     = p0 + F0/F0max * h0
        simulation.context.setPositions(p1)
        state1 = simulation.context.getState(getEnergy=True, getForces=True, getPositions=True)
        E1     = state1.getPotentialEnergy()._value
        
        if E1 < E0:
            report = 'accepted'
            h0 = 1.2 * h0
            
        else:
            report = 'rejected'
            h0 = 0.2 * h0
            simulation.context.setPositions(p0)
            
        print(f'step{i:5d}: {report:s}. E0: {E0:.1e} kJ/mol. E1: {E1:.1e} kJ/mol. Fmax: {F0max:.1e}')
            
    return simulation

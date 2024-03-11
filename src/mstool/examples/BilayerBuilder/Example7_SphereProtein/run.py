import mstool

mstool.SphereProtein(radius=100,
                     protein={mstool.Universe(mstool.MPAA2): 5,
                              mstool.Universe(mstool.GPCR):  5},
                     out='protein.dms')

mstool.SphereBuilder(radius=100,
                     protein='protein.dms',
                     upper={'POPC': 1600, 'DOPS': 990, 'CHL1': 10},
                     lower={'POPC': 1000, 'DOPS': 400, 'CHL1': 10},
                     aa_nsteps=aa_nsteps,
                     remove_solvent=True,
                     rockENM=False)

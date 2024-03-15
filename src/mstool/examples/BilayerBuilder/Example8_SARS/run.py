import mstool

mstool.SphereBuilder(radius=500, 
                     protein={mstool.Universe('S3.pdb'): 30,
                              mstool.Universe('M2.pdb'): 200,
                              mstool.Universe('E5.pdb'): 10},
                     upper={'POPC': 17100, 'DOPS': 17100, 'CHL1': 17100},
                     lower={'POPC': 15100, 'DOPS': 15100, 'CHL1': 15100},
                     solvate=False,
                     water=200,
                     cg_nsteps=10)


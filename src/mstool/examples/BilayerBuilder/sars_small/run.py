import mstool

mstool.SphereBuilder(radius=100, 
                     protein={mstool.Universe('S3.pdb'): 5, 
                              mstool.Universe('M2.pdb'): 2, 
                              mstool.Universe('E5.pdb'): 2},
                     upper={'POPC': 2600},
                     lower={'POPC': 1400},
                     solvate=False,
                     water=500)


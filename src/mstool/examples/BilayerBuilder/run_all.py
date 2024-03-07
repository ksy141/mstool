# copy this script to somewhere else
# and then execute this script
# to save time, I reduced the number of NVT steps (nsteps=1000)
# The default is nsteps=10000, which is 2 ps.

import mstool
import shutil
import os

aa_nsteps = 1000

### Example 1
os.makedirs('Example1_POPC')
os.chdir('Example1_POPC')
mstool.BilayerBuilder(upper={'POPC': 40}, lower={'POPC': 40}, aa_nsteps=aa_nsteps)
os.chdir('../')

### Example 2
os.makedirs('Example2_Lipids')
os.chdir('Example2_Lipids')
lipids = {'DPPC': 5, 'DOPC': 5, 'DMPC': 5, 'DSPC': 5, 'POPC': 5,
          'DOPS': 5, 'POPS': 5, 'POPG': 5, 'DOPG': 5, 'CHL1': 5,
          'POPA': 5, 'DOPA': 5, 'POPE': 5, 'DOPE': 5}
mstool.BilayerBuilder(upper=lipids, lower=lipids, 
                      aa_nsteps=aa_nsteps, remove_solvent=True)
os.chdir('../')

### Example 3
os.makedirs('Example3_Triolein')
os.chdir('Example3_Triolein')
mstool.BilayerBuilder(upper={'POPC':100, 'TRIO':10},
                      lower={'POPC':100, 'TRIO':10},
                      mapping_add=mstool.TRIOMAPPING,
                      martini_add=mstool.TRIOMARTINI,
                      ff_add=mstool.TRIOFF,
                      aa_nsteps=aa_nsteps,
                      remove_solvent=True)
os.chdir('../')

### Example 4
os.makedirs('Example4_MembraneProtein')
os.chdir('Example4_MembraneProtein')
mstool.BilayerBuilder(protein=mstool.MPAA2,
                      upper={'POPC':100},
                      lower={'POPC':100},
                      aa_nsteps=aa_nsteps,
                      remove_solvent=True)
os.chdir('../')


### Example 5
os.makedirs('Example5_GPCR')
os.chdir('Example5_GPCR')
mstool.BilayerBuilder(protein=mstool.GPCR,
                      upper={'POPC':100},
                      lower={'POPC':100},
                      aa_nsteps=aa_nsteps,
                      remove_solvent=True)
os.chdir('../')


### Example 6
os.makedirs('Example6_Sphere')
os.chdir('Example6_Sphere')
mstool.SphereBuilder(radius=60, 
                     upper={'POPC': 1090, 'CHL1': 10},
                     lower={'POPC':  390, 'CHL1': 10},
                     aa_nsteps=aa_nsteps,
                     remove_solvent=True)
os.chdir('../')


### Example 7
os.makedirs('Example7_SphereProtein')
os.chdir('Example7_SphereProtein')
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

os.chdir('../')


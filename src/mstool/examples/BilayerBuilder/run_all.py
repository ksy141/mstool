# copy this script to somewhere else
# and then execute this script
# to save time, I reduced the number of NVT steps (nsteps=1000)
# The default is nsteps=10000, which is 2 ps.

import mstool

aa_nsteps = 1000

### Example 1
mstool.BilayerBuilder(workdir='Example1_POPC', 
                      upper={'POPC': 40}, lower={'POPC': 40}, 
                      aa_nsteps=aa_nsteps)

### Example 2
lipids = {'DPPC': 5, 'DOPC': 5, 'DMPC': 5, 'DSPC': 5, 'POPC': 5,
          'DOPS': 5, 'POPS': 5, 'POPG': 5, 'DOPG': 5, 'CHL1': 5,
          'POPA': 5, 'DOPA': 5, 'POPE': 5, 'DOPE': 5}
mstool.BilayerBuilder(workdir='Example2_Lipids',
                      upper=lipids, lower=lipids, 
                      aa_nsteps=aa_nsteps, remove_solvent=True)

### Example 3
mstool.BilayerBuilder(workdir='Example3_Triolein',
                      upper={'POPC':100, 'TRIO':10},
                      lower={'POPC':100, 'TRIO':10},
                      mapping_add=mstool.TRIOMAPPING,
                      martini_add=mstool.TRIOMARTINI,
                      ff_add=mstool.TRIOFF,
                      aa_nsteps=aa_nsteps,
                      remove_solvent=True)

### Example 4
mstool.BilayerBuilder(workdir='Example4_MembraneProtein',
                      protein=mstool.MPAA2,
                      upper={'POPC':100},
                      lower={'POPC':100},
                      aa_nsteps=aa_nsteps,
                      remove_solvent=True)


### Example 5
mstool.BilayerBuilder(workdir='Example5_GPCR',
                      protein=mstool.GPCR,
                      upper={'POPC':100},
                      lower={'POPC':100},
                      aa_nsteps=aa_nsteps,
                      remove_solvent=True)


### Example 6
upper = {'DPPC': 80, 'DOPC': 80, 'DMPC': 80, 'DSPC': 80, 'POPC': 80,
         'DOPS': 80, 'POPS': 80, 'POPG': 80, 'DOPG': 80, 'CHL1': 80,
         'POPA': 80, 'DOPA': 80, 'POPE': 80, 'DOPE': 80}

lower = {'DPPC': 30, 'DOPC': 30, 'DMPC': 30, 'DSPC': 30, 'POPC': 30,
         'DOPS': 30, 'POPS': 30, 'POPG': 30, 'DOPG': 30, 'CHL1': 30,
         'POPA': 30, 'DOPA': 30, 'POPE': 30, 'DOPE': 30}

mstool.SphereBuilder(workdir='Example6_Sphere',
                     radius=60, 
                     upper=upper,
                     lower=lower,
                     aa_nsteps=aa_nsteps,
                     remove_solvent=True)


### Example 7
mstool.SphereBuilder(workdir='Example7_SphereProtein',
                     radius=100, 
                     protein={mstool.Universe(mstool.MPAA2): 5, 
                              mstool.Universe(mstool.GPCR):  5},
                     upper={'POPC': 1600, 'DOPS': 990, 'CHL1': 10},
                     lower={'POPC': 1000, 'DOPS': 400, 'CHL1': 10},
                     aa_nsteps=aa_nsteps,
                     remove_solvent=True)


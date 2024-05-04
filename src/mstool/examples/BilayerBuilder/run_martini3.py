import mstool

# Although it's feasible to execute a membrane builder using Martini 3.0, I advise against it. 
# Equilibration appears to take longer compared to Martini 2, 
# though this observation likely pertains to my simulation setup and parameters 
# rather than being inherent to Martini 3.

martini = ['/Users/siyoungkim/mstool/src/mstool/FF/martini3.0.0/martini.itp',
           '/Users/siyoungkim/mstool/src/mstool/FF/martini3.0.0/phospholipids.itp',
           '/Users/siyoungkim/mstool/src/mstool/FF/martini3.0.0/protein.itp']

mapping = ['/Users/siyoungkim/mstool/src/mstool/mapping/martini3.protein.c36m.dat',
           '/Users/siyoungkim/mstool/src/mstool/mapping/martini.lipid.c36.dat']

### Example 1
mstool.BilayerBuilder(workdir='Example1_POPC', 
                      upper={'POPC': 100}, lower={'POPC': 100},
                      martini=martini, mapping=mapping, cg_nsteps=500000)

### Example 2
lipids = {'DPPC': 5, 'DOPC': 5, 'POPC': 5,
          'DOPS': 5, 'POPS': 5, 'POPG': 5, 'DOPG': 5,
          'POPA': 5, 'DOPA': 5, 'POPE': 5, 'DOPE': 5}
mstool.BilayerBuilder(workdir='Example2_Lipids',
                      upper=lipids, lower=lipids, 
                      martini=martini, mapping=mapping, cg_nsteps=500000)


### Example 4
mstool.BilayerBuilder(workdir='Example4_MembraneProtein',
                      protein=mstool.MPAA2,
                      upper={'POPC':100, 'POPS': 10, 'POPE': 10},
                      lower={'POPC':100, 'POPS': 10, 'POPE': 10},
                      martini=martini, mapping=mapping, cg_nsteps=500000)



### Example 5
mstool.BilayerBuilder(workdir='Example5_GPCR',
                      protein=mstool.GPCR,
                      upper={'DOPC':100, 'DOPS': 10, 'DOPE': 10},
                      lower={'DOPC':100, 'DOPS': 10, 'DOPE': 10},
                      martini=martini, mapping=mapping, cg_nsteps=500000)



import mstool

### Example 2A
mstool.BilayerBuilder(workdir='Example2A_POPC', 
                      upper={'POPC': 40}, lower={'POPC': 40})

### Example 2B
lipids = {'DPPC': 5, 'DOPC': 5, 'DMPC': 5, 'DSPC': 5, 'POPC': 5,
          'DOPS': 5, 'POPS': 5, 'POPG': 5, 'DOPG': 5, 'CHL1': 5,
          'POPA': 5, 'DOPA': 5, 'POPE': 5, 'DOPE': 5}
mstool.BilayerBuilder(workdir='Example2B_Lipids',
                      upper=lipids, lower=lipids)

### Example 2C
mstool.BilayerBuilder(workdir='Example2C_Triolein',
                      upper={'POPC':100, 'TRIO':6},
                      lower={'POPC':100, 'TRIO':6},
                      mapping_add=mstool.TRIOMAPPING,
                      martini_add=mstool.TRIOMARTINI,
                      ff_add=mstool.TRIOFF)

### Example 2D
mstool.BilayerBuilder(workdir='Example2D_ompF',
                      protein=mstool.MPAA2,
                      upper={'POPC':100, 'POPS': 10, 'POPE': 10, 'CHL1': 10},
                      lower={'POPC':100, 'POPS': 10, 'POPE': 10, 'CHL1': 10})


### Example 2E
mstool.BilayerBuilder(workdir='Example2E_3SN6',
                      protein=mstool.GPCR,
                      upper={'DOPC':100, 'DOPS': 10, 'DOPE': 10, 'CHL1': 10},
                      lower={'DOPC':100, 'DOPS': 10, 'DOPE': 10, 'CHL1': 10})

### Example 2F
mstool.BilayerBuilder(workdir='Example2F_4XNV',
                      protein=mstool.GPCR2,
                      upper={'DOPC':100, 'DOPS': 10, 'DOPE': 10, 'CHL1': 10},
                      lower={'DOPC':100, 'DOPS': 10, 'DOPE': 10, 'CHL1': 10})


### Example 3A
mstool.BilayerBuilder(workdir='Example3A_crowded', 
                      protein={mstool.MPAA2: 4, mstool.GPCR: 5}, 
                      upper={'POPC':1200}, lower={'POPC': 1200}, 
                      cg_nsteps=2000000)


### Example 3B
mstool.BilayerBuilder(workdir='Example3B_LD',
                      upper={'POPC':135, 'CHL1':5},
                      lower={'POPC':135, 'CHL1':5},
                      between={'TRIO': 200, 'CHYO': 200},
                      sep=80.0, ff_add=mstool.TRIOFF,
                      martini_add=mstool.TRIOMARTINI,
                      mapping_add=mstool.TRIOMAPPING)


### Example 3B LD2
# To build a lipid droplet whose radius is 50A,
# radius = 50/2 = 25 and sep = 50
mstool.SphereBuilder(workdir='Example3B_LD2',
                     radius=25, sep=50,
                     upper={'POPC': 700, 'TRIO': 40},
                     between={'TRIO': 300},
                     martini_add=mstool.TRIOMARTINI,
                     mapping_add=mstool.TRIOMAPPING,
                     ff_add=mstool.TRIOFF)

### Example 3C
mstool.BilayerBuilder(workdir='Example3C_Seipin',
                      protein=mstool.SEIPIN,
                      upper={'POPC':437, 'TRIO': 27},
                      lower={'POPC':360, 'TRIO': 21},
                      martini_add=mstool.TRIOMARTINI,
                      mapping_add=mstool.TRIOMAPPING,
                      ff_add=mstool.TRIOFF,
                      cg_nsteps=500000,
                      dx=7.2, removedr=3.0)


### Example 3C (is it possible to get rid of all lipids inside protein?)
mstool.BilayerBuilder(workdir='Example3C_Seipin2',
                      protein=mstool.SEIPIN,
                      upper={'POPC':437, 'TRIO': 27},
                      lower={'POPC':360, 'TRIO': 21},
                      martini_add=mstool.TRIOMARTINI,
                      mapping_add=mstool.TRIOMAPPING,
                      ff_add=mstool.TRIOFF,
                      cg_nsteps=500000,
                      dx=7.2, removedr=3.0,
                      rcut_use_enclosed_protein=True)

### Example 4D
mstool.SphereBuilder(workdir='Example4D_SphereProtein',
                     radius=100, 
                     protein={mstool.Universe(mstool.MPAA2): 5, 
                              mstool.Universe(mstool.GPCR):  5},
                     upper={'POPC': 1600, 'DOPS': 990, 'CHL1': 10},
                     lower={'POPC': 1000, 'DOPS': 400, 'CHL1': 10})



import mstool
mstool.BilayerBuilder(protein=mstool.SEIPIN,
                      upper={'POPC':437, 'TRIO': 27},
                      lower={'POPC':360, 'TRIO': 21},
                      martini_add=mstool.TRIOMARTINI,
                      mapping_add=mstool.TRIOMAPPING,
                      ff_add=mstool.TRIOFF,
                      cg_nsteps=500000,
                      dx=7.2, removedr=3.0)

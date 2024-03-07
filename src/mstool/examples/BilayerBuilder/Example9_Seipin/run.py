import mstool
import shutil
shutil.copyfile(mstool.SEIPIN, 'protein.pdb')
mstool.BilayerBuilder(protein='protein.pdb',
                      upper={'POPC':437, 'TRIO': 27},
                      lower={'POPC':360, 'TRIO': 21},
                      martini_add=mstool.TRIOMARTINI,
                      mapping_add=mstool.TRIOMAPPING,
                      ff_add=mstool.TRIOFF,
                      cg_nsteps=500000)


import mstool
import shutil
shutil.copyfile(mstool.MPAA2, 'protein.pdb')
mstool.BilayerBuilder(protein='protein.pdb',
                      upper={'POPC':100},
                      lower={'POPC':100})

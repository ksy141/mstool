import mstool
from mstool.core.gmx2dms import GMX2DMS
d = mstool.ReadMartini('topol.top')
GMX2DMS(structure='../relaxed.pdb', martini=d, out='output.dms')

mstool.getDMSEnergy(dms_in='output.dms')
mstool.getGMXEnergy('../relaxed.pdb', add='-nt 1')



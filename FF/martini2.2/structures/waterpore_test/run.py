# python ../makepore.py 8

import mstool

d = mstool.ReadMartini('topol.top', define={'POSRES_Z': 'True', 'POSRES_XYZ': 'True'})
mstool.GMX2DMS('P008.pdb', d, 'P008.dms')
mstool.getGMXEnergy('P008.pdb', p='topol.top', add='-nt 1')
mstool.getDMSEnergy('P008.dms', nonbondedMethod='CutoffNonPeriodic')
mstool.runMartiniEM('P008.dms', 'P008.EM.dms', nonbondedMethod='CutoffNonPeriodic')


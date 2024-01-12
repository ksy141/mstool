import mstool
m = mstool.ReadMartini('topol.top', define={'POSRES_XYZ': '1000.0', 'FLEXIBLE': 'True'})
mstool.GMX2DMS('../CHOL.pdb', martini=m, out='chol.martini.dms')
mstool.runMartiniEM( dms_in='chol.martini.dms', out='chol.EM.pdb', soft=False, nonbondedMethod='CutoffNonPeriodic')

m = mstool.ReadMartini('topol.top', define={'POSRES_XYZ': '1000.0'})
mstool.GMX2DMS('chol.EM.pdb', martini=m, out='chol.martini2.dms')
mstool.runMartiniEM( dms_in='chol.martini2.dms', out='chol.EM2.pdb', soft=False, nonbondedMethod='CutoffNonPeriodic')




import mstool
mstool.Map('1J46.pdb', '1J46_CG.pdb')
mstool.Backmap('1J46_CG.pdb', pbc=False, nsteps=0)

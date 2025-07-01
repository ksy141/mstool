import mstool
mstool.Map('1RNA.pdb', '1RNA_CG.pdb')
mstool.Backmap('1RNA_CG.pdb', pbc=False, nsteps=0)

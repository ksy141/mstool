import mstool
mstool.Map('1BNA.pdb', '1BNA_CG.pdb')
mstool.Backmap('1BNA_CG.pdb', pbc=False, nsteps=0)

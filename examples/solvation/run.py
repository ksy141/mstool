import mstool
mstool.SolvateMartini(out='waterbox.pdb', dimensions=[42, 42, 42])
mstool.SolvateMartini(structure='cg.pdb', out='solvated.pdb', center=True)
mstool.SolvateMartini(structure='cg.pdb', out='solvatedt.pdb', t=20)


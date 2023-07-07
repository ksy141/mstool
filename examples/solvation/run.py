import mstool
mstool.SolvateMartini(out='waterbox.pdb', dimensions=[42, 42, 42])
mstool.SolvateMartini(structure='cg.pdb', out='solvated.pdb', center=True)
mstool.SolvateMartini(structure='cg.pdb', out='solvatedt.pdb', t=20)


waterbox1 = mstool.SolvateMartini(dimensions=[50, 50, 50], pos='POT', neg='CLA')
waterbox2 = mstool.SolvateMartini(dimensions=[50, 50, 50], pos='MG',  neg='CLA')

waterbox3 = mstool.SolvateMartini(dimensions=[50, 50, 50], pos='POT', neg='CLA')
mstool.ionize(waterbox3, pos='MG', neg='CLA', out='waterbox3.pdb')


import mstool

#mstool.Backmap('lipid/input_cg_lipid.dms', 'lipid/input_aa_lipid.dms', mapping='mapping.dat')
mstool.REM(rock='protein/workdir/step9_final.dms', 
           structure='lipid/input_aa_lipid.dms', 
           out='final.dms', 
           mapping='mapping.dat',
           ff_add='trio.xml',
           rockrcut=0.5)



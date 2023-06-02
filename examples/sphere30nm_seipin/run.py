import mstool
import shutil
import os

workdir = 'workdir'
#if os.path.exists(workdir):
#    shutil.rmtree(workdir)
#os.mkdir(workdir)
#
#mstool.Backmap('ext_306_lipid.dms', workdir + '/input_aa_lipid.dms', mapping='mapping.dat')

#u = mstool.Universe(workdir + '/input_aa_lipid.dms')
#u.cell = [[900, 0, 0], [0, 900, 0], [0, 0, 900]]
#u.write(workdir + '/input_aa_lipid2.dms')

mstool.REM(protein='protein/workdir/step9_final.pdb',
           structure=workdir + '/input_aa_lipid2.dms', 
           out=workdir + '/final.dms', 
           mapping_add='mapping.dat',
           ff_add='trio.xml', 
           pbc=False)

mstool.CheckStructure(workdir + '/final.dms', log=workdir + '/log.txt')
mstool.Universe(workdir + '/final.dms').write(workdir + '/final.pdb')



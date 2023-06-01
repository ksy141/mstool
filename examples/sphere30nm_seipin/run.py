import mstool
import shutil
import os

workdir = 'workdir'
#if os.path.exists(workdir):
#    shutil.rmtree(workdir)
#os.mkdir(workdir)
#
#mstool.Backmap('ext_306_lipid.dms', workdir + '/input_aa_lipid.dms', mapping='mapping.dat')
mstool.REM(protein='protein/workdir/step9_final.pdb',
           structure=workdir + '/input_aa_lipid.dms', 
           out=workdir + '/final.dms', 
           mapping='mapping.dat',
           ff_add='trio.xml', 
           pbc=False)

msprot.CheckStructure(workdir + '/final.dms', mapping='mapping.dat', log=workdir + '/log.txt')
msprot.Universe(workdir + '/final.dms').write(workdir + '/final.pdb')



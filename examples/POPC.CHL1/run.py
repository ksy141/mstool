import msprot
import shutil
import os

workdir = 'workdir'
if os.path.exists(workdir):
    shutil.rmtree(workdir)
os.mkdir(workdir)

msprot.Backmap('cg.pdb', workdir + '/aa.pdb')
msprot.REM(workdir + '/aa.pdb', workdir + '/aa_final.pdb')
msprot.CheckStructure(workdir + '/aa_final.pdb')


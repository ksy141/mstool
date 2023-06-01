import msprot
import shutil
import os

workdir = 'workdir'
if os.path.exists(workdir):
    shutil.rmtree(workdir)
os.mkdir(workdir)

msprot.Map('filtered.pdb', workdir + '/cg.pdb')
msprot.Backmap(workdir + '/cg.pdb',  workdir + '/aa.pdb')
msprot.REM(workdir + '/aa.pdb', workdir + '/aa_final.pdb')
msprot.CheckStructure(workdir + '/aa_final.pdb')


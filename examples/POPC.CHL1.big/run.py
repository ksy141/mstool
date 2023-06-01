import mstool
import shutil
import os

workdir = 'workdir'
if os.path.exists(workdir):
    shutil.rmtree(workdir)
os.mkdir(workdir)

mstool.Backmap('cg.pdb', workdir + '/aa.pdb')
mstool.REM(workdir + '/aa.pdb', workdir + '/aa_final.pdb')
mstool.CheckStructure(workdir + '/aa_final.pdb')


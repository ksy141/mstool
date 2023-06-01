import mstool
import shutil
import os

workdir = 'workdir'
if os.path.exists(workdir):
    shutil.rmtree(workdir)
os.mkdir(workdir)

mstool.Map('filtered.pdb', workdir + '/cg.pdb')
mstool.Backmap(workdir + '/cg.pdb',  workdir + '/aa.pdb')
mstool.REM(workdir + '/aa.pdb', workdir + '/aa_final.pdb')
mstool.CheckStructure(workdir + '/aa_final.pdb')


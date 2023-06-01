import msprot
import shutil
import os

workdir = 'workdir'
if os.path.exists(workdir):
    shutil.rmtree(workdir)
os.mkdir(workdir)

msprot.Backmap('ext_449.pdb', workdir + '/aa.pdb', mapping_add='mapping.dat')
msprot.REM(workdir + '/aa.pdb', workdir + '/aa_final.pdb', mapping_add='mapping.dat', ff_add='trio.xml')
msprot.CheckStructure(workdir + '/aa_final.pdb')


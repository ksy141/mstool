import mstool
import shutil
import os

workdir = 'workdir'
if os.path.exists(workdir):
    shutil.rmtree(workdir)
os.mkdir(workdir)

mstool.REM(rock='protein.pdb', structure='aa_lipid.pdb', out=workdir + '/aa_final.pdb')
mstool.CheckStructure(workdir + '/aa_final.pdb')

shutil.move('ROCK.xml',     workdir + '/')
shutil.move('ROCK.dms',     workdir + '/')
shutil.move('ROCK_rem.pdb', workdir + '/')


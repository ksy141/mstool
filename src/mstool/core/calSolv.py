import subprocess
import os
import numpy as np
import time
from  .martinizedms import MartinizeDMS
from  .readmartini import ReadMartini
from  .universe     import Universe
from  .solvate_martini import OctanolMartini
from  .dmsfile import DMSFile 

def calSolv(resname=None, workdir='workdir', 
        octfrac=0.7, t=50, martini=None, T=310, nsteps=10000):

    if martini is None:
        martini = ReadMartini(martini3=True)
    os.makedirs(workdir)


    ### Add the molecule of interest at the center
    if resname:
        mol   = martini.martini['molecules'][resname]
        names = mol['atoms']['name']
        
        data = {'name': [], 'chain': 'A', 'resname': resname, 
                'x':    [], 'y': [], 'z': []}
        
        for name in names:
            data['name'].append(name)
            x, y, z = np.random.rand(3)
            data['x'].append(x)
            data['y'].append(y)
            data['z'].append(z)
        u = Universe(data)

    else:
        u = Universe()

    u.dimensions = np.array([t, t, t, 90., 90., 90.])
    u.cell       = [[t, 0, 0], [0, t, 0], [0, 0, t]]
    u.write(f'{workdir}/step1_mol.dms')
    

    ### Octaonlize & Save typed DMS
    newu = OctanolMartini(structure=u, molfrac=octfrac)
    newu.atoms[['x','y','z']] += newu.dimensions[0:3]/2
    newu.write(f'{workdir}/step2_solv.dms')
    MartinizeDMS(f'{workdir}/step2_solv.dms', martini=martini, out=f'{workdir}/step3_ff.dms')


    ### Run REM
    dms = DMSFile(f'{workdir}/step3_ff.dms')
    dms.createSystem(REM=True, tapering='shift', martini=True, 
                     nonbondedCutoff=1.1, nonbondedMethod='CutoffPeriodic',
                     improper_prefactor=0.99, removeCMMotion=True)
    dms.runEMNPT(dt=0.02, nsteps=nsteps, frictionCoeff=10.0, 
                 barfreq=10, dcdfreq=100, csvfreq=100, tension=0, P=1.0, T=T, 
                 semiisotropic=False, out=f'{workdir}/step4_rem.dms')


    ### Run EM
    dms.createSystem(REM=False, tapering='shift', martini=True, 
                     nonbondedCutoff=1.1, nonbondedMethod='CutoffPeriodic', 
                     improper_prefactor=0.99, removeCMMotion=True)
    dms.runEMNPT(dt=0.02, nsteps=nsteps, frictionCoeff=10.0, 
                 barfreq=10, dcdfreq=100, csvfreq=100, tension=0, P=1.0, T=T, 
                 semiisotropic=False, out=f'{workdir}/step5_em.dms')
   

    ### Make Top
    u = Universe(f'{workdir}/step5_em.dms')
    u.write(f'{workdir}/step5_em.gro')
    u.makeTopology(top=f'{workdir}/topol.top', martini=martini)




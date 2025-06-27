from .universe import Universe, Merge, Wrap
from .backmap import Backmap
from .rem import REM
from ..lib.distance import distance_matrix
import os
import numpy as np
import time
import pickle
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

class DivideConquer:
    def __init__(self, structure, Nx, Ny, Nz, workdir='workdir', AA=None, pbc=True, wrap=False, std_threshold=80):
        '''As it relies on cog, lipids should be whole/unwrapped. However, proteins should be unwrapped. systems should be [0, L]'''
        self.structure = structure
        self.workdir   = workdir
        self.wrap      = wrap
        self.AA = AA
        self.Nx = Nx
        self.Ny = Ny
        self.Nz = Nz
        
        self.translate = np.zeros(3)
        u = Universe(structure)
        if wrap:
            u = Wrap(u, std_threshold=std_threshold)

        if u.dimensions is None or u.dimensions[0] * u.dimensions[1] * u.dimensions[2] == 0 or not pbc:
            dimensions     = u.atoms[['x','y','z']].max(axis=0).to_numpy() \
                           - u.atoms[['x','y','z']].min(axis=0).to_numpy() \
                           + np.array([30, 30, 30])
            center         = u.atoms[['x','y','z']].max(axis=0).to_numpy() \
                           + u.atoms[['x','y','z']].min(axis=0).to_numpy()
            center        /= 2
            u.dimensions   = np.array([dimensions[0], dimensions[1], dimensions[2], 90, 90, 90])
            u.cell         = np.array([[dimensions[0], 0, 0], [0, dimensions[1], 0], [0, 0, dimensions[2]]])
            self.translate = -center + u.dimensions[0:3] / 2
            u.atoms[['x','y','z']] += self.translate
        
        self.dimensions = u.dimensions
        self.cell = u.cell
        self.u = u

        # Divide
        self.dx = u.dimensions[0] / self.Nx
        self.dy = u.dimensions[1] / self.Ny
        self.dz = u.dimensions[2] / self.Nz


    def create(self, shift=False, protein_distance_cutoff=25, **kwargs):
        '''use shift=True if center is origin [0, 0, 0]'''

        time1 = time.time()

        # Make a directory
        os.makedirs(self.workdir, exist_ok=True)
       
        # Should backmap protein first
        u = self.u
        if shift:
            u.atoms[['x','y','z']] += u.dimensions[0:3] / 2
            self.translate += u.dimensions[0:3] / 2
        
        if np.all(self.translate == np.zeros(3)) and not self.wrap:
            # writing takes some time. Do not write unless necessary
            pass
        else:
            u.write(f'{self.workdir}/input.pdb')
            u.write(f'{self.workdir}/input.dms')

        protein = Universe(u.select('protein'))
        protein.dimensions = u.dimensions
        if len(protein.atoms) > 0:
            os.makedirs(self.workdir + '/protein', exist_ok=True)
            protein.write(self.workdir + '/protein/protein.pdb')
            protein.write(self.workdir + '/protein/protein.dms')

            protein_pos = {}
            for chain, group in protein.select('protein & @BB').groupby('chain'):
                protein_pos[chain] = group[['x','y','z']].to_numpy()
            
            protein_dm = np.zeros((len(protein_pos.keys()), len(protein_pos.keys())))
            for i, (chaini, posi) in enumerate(protein_pos.items()):
                for j, (chainj, posj) in enumerate(protein_pos.items()):
                    if i >= j: continue
                    protein_dm[i, j] = np.min(distance_matrix(posi, posj, u.dimensions))
                    protein_dm[j, i] = protein_dm[i, j]
            
            Z = linkage(squareform(protein_dm))
            labels = fcluster(Z, t=protein_distance_cutoff, criterion='distance')
            unique_labels = np.unique(labels)

            for label in unique_labels:
                os.makedirs(f'{self.workdir}/protein/{label}', exist_ok=True)
                label_chains  = np.array(list(protein_pos.keys()))[label == labels]
                label_protein = Universe(protein.atoms[protein.atoms['chain'].isin(label_chains)])
                label_protein.dimensions = u.dimensions
                label_protein.cell = u.cell
                label_protein.write(f'{self.workdir}/protein/{label}/protein.dms')

                flipped = 100
                for i in range(3):
                    try:
                        # sometimes NaN
                        bm = Backmap(structure = f'{self.workdir}/protein/{label}/protein.dms', 
                                     workdir   = f'{self.workdir}/protein/{label}/workdir',
                                     turn_off_EMNVT=True, use_existing_workdir=True, **kwargs)
                        flipped = bm.check.output['peptide'] + bm.check.output['chiral'] 
                        if flipped == 0:
                            break

                    except:
                        pass
            
            proteinAA = Universe()
            for label in unique_labels:
                proteinAA = Merge(proteinAA.atoms, Universe(f'{self.workdir}/protein/{label}/workdir/step4_final.dms').atoms)
            proteinAA.dimensions = u.dimensions
            proteinAA.cell = u.cell


        else:
            proteinAA = Universe()


        # Combine protein + AA
        if self.AA is not None:
            AA = Universe(self.AA)
            AA.atoms[['x','y','z']] += self.translate
            newAA = Merge(proteinAA.atoms, AA.atoms)
        else:
            newAA = proteinAA

        newAA.dimensions = u.dimensions
        newAA.cell = u.cell
        if len(newAA.atoms) > 0:
            newAA.write(f'{self.workdir}/AA.dms')
            newAA.write(f'{self.workdir}/AA.pdb')

        # protein + AA atoms must be within [0, L]
        newAA.atoms[['x','y','z']] -= u.dimensions[0:3] * np.floor(newAA.atoms[['x','y','z']].to_numpy() / u.dimensions[0:3])
        #center = u.atoms[['x','y','z']].mean(axis=0)
        #if np.linalg.norm(center) > np.linalg.norm(center - u.dimensions[0:3] / 2):
        #    # wrap between [0, L]
        #    newAA.atoms[['x','y','z']] -= u.dimensions[0:3] * np.floor(newAA.atoms[['x','y','z']].to_numpy() / u.dimensions[0:3])
        #else:
        #    # wrap between [-L/2, L/2]
        #    newAA.atoms[['x','y','z']] -= u.dimensions[0:3] * np.round(newAA.atoms[['x','y','z']].to_numpy() / u.dimensions[0:3])

        # Non-protein
        lipids  = Universe(u.atoms[~u.select('protein', returnbA=True)])
        lipids.dimensions = u.dimensions
        cog = lipids.atoms.groupby('resn').agg({'x': 'mean', 'y': 'mean', 'z': 'mean'})

        def partition(i, dx, Nx, pos):
            bAx = None
            if i == 0:
                bAx = (pos < (i+1) * dx)
            elif i == Nx-1:
                bAx = (i * dx <= pos)
            else:
                bAx = (i * dx <= pos) & (pos < (i+1) * dx)
            return bAx

        for i in range(self.Nx):
            for j in range(self.Ny):
                for k in range(self.Nz):
                    subworkdir = f'{self.workdir}/{i}_{j}_{k}'
                    os.makedirs(subworkdir, exist_ok=True)
                    
                    # lipids
                    bAx = partition(i, self.dx, self.Nx, cog['x'])
                    bAy = partition(j, self.dy, self.Ny, cog['y'])
                    bAz = partition(k, self.dz, self.Nz, cog['z'])
                    resn = cog[bAx & bAy & bAz]

                    sublipids = Universe(lipids.atoms[lipids.atoms['resn'].isin(resn.index)])
                    if len(sublipids.atoms) > 0:
                        sublipids.write(subworkdir + '/input.dms')
                    
                    # ROCK
                    bAx = partition(i, self.dx, self.Nx, newAA.atoms['x'])
                    bAy = partition(j, self.dy, self.Ny, newAA.atoms['y'])
                    bAz = partition(k, self.dz, self.Nz, newAA.atoms['z'])
                    subAA = Universe(newAA.atoms[newAA.atoms['resn'].isin(newAA.atoms[bAx & bAy & bAz]['resn'])])
                    if len(subAA.atoms) > 0:
                        subAA.write(subworkdir + '/AA.dms')

        time2 = time.time()
        with open(f'{self.workdir}/time_create.txt', 'w') as W:
            W.write(f'{(time2 - time1)/60:.3f} min\n')

    def _run_subregion(self, i, j, k, Kwall=1000, delta=0.1, rerun_unsuccesful=False, **kwargs):
        """
        Work on one (i,j,k) sub-box.
        All *extra* keyword arguments are forwarded unchanged to Backmap.
        """
        subdir = f'{self.workdir}/{i}_{j}_{k}'

        wall = {
            "x0": i * self.dx * 0.1, "x1": (i + 1) * self.dx * 0.1,
            "y0": j * self.dy * 0.1, "y1": (j + 1) * self.dy * 0.1,
            "z0": k * self.dz * 0.1, "z1": (k + 1) * self.dz * 0.1,
            "Kwall": Kwall, "delta": delta,
        }
        
        # Save dictionary to a file
        with open(f'{subdir}/wall.pkl', 'wb') as f:
            pickle.dump(wall, f)

        ## Load dictionary from a file
        #with open('{subdir}/wall.pkl', 'rb') as f:
        #    wall = pickle.load(f)

        if rerun_unsuccesful and os.path.exists(f'{subdir}/workdir/step4_final.dms'):
            return
        
        if os.path.exists(f'{subdir}/input.dms'):
            flipped = 100
            for idx in range(3):
                bm = Backmap(
                        structure=f'{subdir}/input.dms',
                        workdir=f'{subdir}/workdir',
                        AA=f'{subdir}/AA.dms' if os.path.exists(f'{subdir}/AA.dms') else None,
                        pbc=False, wall=wall, turn_off_EMNVT=True, use_existing_workdir=True, **kwargs
                )

                flipped = bm.check.output['cistrans'] + bm.check.output['chiral']
                if flipped == 0:
                    break


    def run(self, Kwall=1000, delta=0.1, rerun_unsuccesful=False, **kwargs):
        for i in range(self.Nx):
            for j in range(self.Ny):
                for k in range(self.Nz):
                    self._run_subregion(i, j, k, Kwall, delta, rerun_unsuccesful=rerun_unsuccesful, **kwargs)


    def runParallel(self, Kwall=1000, delta=0.1, num_workers=None, rerun_unsuccesful=False, **kwargs):
        """
        Parallelised version of the original `run`.  All keyword arguments
        that are *not* explicitly consumed (`Kwall`, `delta`, `num_workers`)
        are cascaded to every `Backmap` call via **kwargs.
        """
        from functools import partial
        from itertools import product
        from multiprocessing import Pool, cpu_count

        if num_workers is None:
            num_workers = min(cpu_count(), self.Nx * self.Ny * self.Nz)

        # Pre-bind the constant parameters (Kwall, delta, and the user’s **kwargs)
        task = partial(self._run_subregion, Kwall=Kwall, delta=delta, rerun_unsuccesful=rerun_unsuccesful, **kwargs)

        with Pool(processes=num_workers) as pool:
            pool.starmap(task, product(range(self.Nx), range(self.Ny), range(self.Nz)))


    def combine(self, **kwargs):
        final = Universe()
        for i in range(self.Nx):
            for j in range(self.Ny):
                for k in range(self.Nz):
                    if os.path.exists(f'{self.workdir}/{i}_{j}_{k}/workdir/step4_nonprotein.dms'):
                        final = Merge(final.atoms, Universe(f'{self.workdir}/{i}_{j}_{k}/workdir/step4_nonprotein.dms').atoms)
                    elif os.path.exists(f'{self.workdir}/{i}_{j}_{k}/workdir/step4_final.dms'):
                        final = Merge(final.atoms, Universe(f'{self.workdir}/{i}_{j}_{k}/workdir/step4_final.dms').atoms)
                    else:
                        continue

        final.dimensions = self.dimensions
        final.cell = self.cell
        final.write(f'{self.workdir}/step1_ungroup.dms')

        Backmap(structure=f'{self.workdir}/step1_ungroup.dms',
                workdir=self.workdir,
                AA=f'{self.workdir}/AA.pdb' if os.path.exists(f'{self.workdir}/AA.pdb') else None,
                use_existing_workdir=True,
                skip_ungroup=True, **kwargs)


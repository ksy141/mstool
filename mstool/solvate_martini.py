from   .utils.universe import Universe
from   .utils.distance import distance_matrix
import numpy as np
import pandas as pd

class SolvateMartini:
    def __init__(self, out, dimensions=None, structure=None, t=None, 
    solventdr=4.93, removedr=5.0, waterslab=0.8, center=True):
        
        # make a water box
        if dimensions:
            u = Universe()
            u.cell = np.array([[dimensions[0], 0, 0],
                               [0, dimensions[1], 0], 
                               [0, 0, dimensions[2]]])

            u.dimensions = np.array([dimensions[0], dimensions[1], dimensions[2], 
                                     90, 90, 90])
            wateru = self.solvate(u, solventdr, removedr, waterslab, center)
            wateru.write(out)
            return None

        u = Universe(structure)

        if t:
            # make new dimensions / cell based on solute particles
            # exisitng dimensions / cell do not matter
            maxx = u.atoms['x'].max() - u.atoms['x'].min()
            maxy = u.atoms['y'].max() - u.atoms['y'].min()
            maxz = u.atoms['z'].max() - u.atoms['z'].min()
            maxd = max(maxx, maxy, maxz)
            dim  = maxd + 2 * t
            u.dimensions = np.array([dim] * 3 + [90] * 3)
            u.cell = np.array([[dim, 0, 0], [0, dim, 0], [0, 0, dim]])
        
        # make a water box    
        wateru = self.solvate(u, solventdr, removedr, waterslab, center)

        pbc = u.dimensions[0:3]
        zero_dims = pbc[0] * pbc[1] * pbc[2]
        assert zero_dims != 0, 'check your dimensions'

        # calculate a distance matrix between water atoms and solute atoms
        dm = distance_matrix(wateru.atoms[['x','y','z']].to_numpy(), 
                             u.atoms[['x','y','z']].to_numpy(), 
                             dimensions=pbc)

        # bA is a selection for the water atoms that overlap with solute atoms
        bA = np.any(dm < removedr, axis=1)
        
        # re-assign residue id
        wateru.atoms.loc[~bA, 'resid'] = np.arange(1, len(wateru.atoms.loc[~bA]) + 1)
        
        # combine and save
        u.atoms = pd.concat([u.atoms, wateru.atoms[~bA]], ignore_index=True)
        u.write(out)


        
    def solvate(self, u, solventdr, removedr, waterslab, center):
        Nx  = u.dimensions[0]               // solventdr
        Nx2 = (u.dimensions[0] - waterslab) // solventdr
        if Nx != Nx2: Nx = Nx2

        Ny  = u.dimensions[1]               // solventdr
        Ny2 = (u.dimensions[1] - waterslab) // solventdr
        if Ny != Ny2: Ny = Ny2

        Nz  = u.dimensions[2]               // solventdr
        Nz2 = (u.dimensions[2] - waterslab) // solventdr
        if Nz != Nz2: Nz = Nz2

        xx = np.arange(Nx) * solventdr
        yy = np.arange(Ny) * solventdr
        zz = np.arange(Nz) * solventdr
        
        if center:
            xx -= u.dimensions[0] / 2
            yy -= u.dimensions[1] / 2
            zz -= u.dimensions[2] / 2
        
        # increment in every direction (very useful)
        xyz = np.array(np.meshgrid(xx, yy, zz)).T.reshape(-1, 3)

        data = {'resname': 'W', 'name': 'W', 'resid': np.arange(1, len(xyz) + 1), 'chain': 'W', 
                'x': xyz[:,0], 'y': xyz[:,1], 'z': xyz[:,2]}

        wateru = Universe(data=data)
        wateru.dimensions = u.dimensions
        wateru.cell = u.cell

        return wateru


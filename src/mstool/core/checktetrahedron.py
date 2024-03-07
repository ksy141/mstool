import numpy as np
from   openmm.app import *
from   .universe import Universe
from   ..utils.openmmutils import getTopPDB

def CheckTetrahedron(structure, ff=[], ff_add=[], tol=1.0):
    pos = Universe(structure).positions
    pdb = getTopPDB(structure, ff=ff, ff_add=ff_add)
    pdbatoms = [atom  for atom in pdb.topology.atoms()]
    bonded   = [set() for atom in pdb.topology.atoms()]

    for bond in pdb.topology.bonds():
        i1 = bond.atom1.index
        i2 = bond.atom2.index
        bonded[i1].add(i2)
        bonded[i2].add(i1)

    bondedlist = [list(b) for b in bonded]
    
    print("####################################################################")
    print("Incorrect tetrahedron geometry is harmless")
    print("because it will be fixed within 0.1 ns of AA MD simulations.")
    print("Checking anyway...")
    print("Tetrahedron checking - started")
    for center, bb in enumerate(bondedlist):
        if len(bb) != 4: continue
        posc = pos[center]
        pos0 = pos[bb[0]]
        pos1 = pos[bb[1]]
        pos2 = pos[bb[2]]
        pos3 = pos[bb[3]]
        
        dr0 = pos0 - posc
        dr1 = pos1 - posc
        dr2 = pos2 - posc
        dr3 = pos3 - posc

        dr0 /= np.linalg.norm(dr0)
        dr1 /= np.linalg.norm(dr1)
        dr2 /= np.linalg.norm(dr2)
        dr3 /= np.linalg.norm(dr3)

        dr = dr0 + dr1 + dr2 + dr3 
        if np.linalg.norm(dr) > tol:
            print(pdbatoms[center], np.linalg.norm(dr))
    print("Tetrahedron checking - finished")
    print("####################################################################")

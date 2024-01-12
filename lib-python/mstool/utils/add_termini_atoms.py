from ..core.universe import Universe
from .protein_sel import *
import pandas as pd

def addTerminiAtoms(structure, out=None):
    u = Universe(structure)
    resnames = three2one.keys()
    bA = u.atoms.resname.isin(resnames)
    proteinu = Universe(data=u.atoms[bA])
    otheru   = Universe(data=u.atoms[~bA])

    chains = list(dict.fromkeys(proteinu.atoms['chain']))
    data = {'name':    [],
            'resname': [],
            'resid':   [],
            'chain':   [],
            'segname': [],
            'x':       [],
            'y':       [],
            'z':       [],
            'bfactor': 0.0}


    for chain in chains:
        bA1 = proteinu.atoms.chain == chain
        Ntermresid = int(proteinu.atoms[bA1]['resid'].min())
        Ctermresid = int(proteinu.atoms[bA1]['resid'].max())
        
        # Ntermini
        bA2 = proteinu.atoms.resid == Ntermresid
        bA3 = proteinu.atoms.name.isin(['HT2', 'HT3', 'H2', 'H3'])
        if len(proteinu.atoms[bA1 & bA2 & bA3]) == 0:
            bA4  = proteinu.atoms.name == 'N'
            ref  = proteinu.atoms[bA1 & bA2 & bA4]
            if len(ref) == 0: assert 0 == 1, 'no N atom in the Ntermini'

            resn, segn, x, y, z = proteinu.atoms[bA1 & bA2 & bA4][['resname', 'segname', 'x', 'y', 'z']].values[0]
            #if resn == 'PRO':
            #    data['name'].extend(['HT2'])
            #    data['resname'].extend([resn])
            #    data['segname'].extend([segn])
            #    data['resid'].extend([Ntermresid])
            #    data['chain'].extend([chain, chain])
            #    data['x'].extend([x + 0.1])
            #    data['y'].extend([y + 0.1])
            #    data['z'].extend([z + 0.1])

            if resn == 'PRO':
                data['name'].extend(['HT1', 'HT2'])
                data['resname'].extend([resn, resn])
                data['segname'].extend([segn, segn])
                data['resid'].extend([Ntermresid, Ntermresid])
                data['chain'].extend([chain, chain])
                data['x'].extend([x + 0.1, x + 0.2])
                data['y'].extend([y + 0.1, y + 0.2])
                data['z'].extend([z + 0.1, z + 0.2])

            else:
                data['name'].extend(['HT2', 'HT3'])
                data['resname'].extend([resn, resn])
                data['segname'].extend([segn, segn])
                data['resid'].extend([Ntermresid, Ntermresid])
                data['chain'].extend([chain, chain])
                data['x'].extend([x + 0.1, x + 0.2])
                data['y'].extend([y + 0.1, y + 0.2])
                data['z'].extend([z + 0.1, z + 0.2])


        # Ctermini
        bA2 = proteinu.atoms.resid == Ctermresid
        bA3 = proteinu.atoms.name.isin(['OT2'])
        if len(proteinu.atoms[bA1 & bA2 & bA3]) == 0:
            bA4 = proteinu.atoms.name == 'C'
            resn, segn, x, y, z = proteinu.atoms[bA1 & bA2 & bA4][['resname', 'segname', 'x','y','z']].values[0]

            data['name'].extend(['OT2'])
            data['resname'].extend([resn])
            data['segname'].extend([segn])
            data['resid'].extend([Ctermresid])
            data['chain'].extend([chain])
            data['x'].extend([x + 0.1])
            data['y'].extend([y + 0.1])
            data['z'].extend([z + 0.1])
    

    if len(data['resid']) == 0:
        if out: u.write(out)
        return u

    protein = pd.concat([proteinu.atoms, pd.DataFrame(data)])
    protein.sort_values(by=['chain', 'resid', 'name'], inplace=True)
    u.atoms = pd.concat([protein, otheru.atoms], ignore_index=True)

    if out: u.write(out)
    return u


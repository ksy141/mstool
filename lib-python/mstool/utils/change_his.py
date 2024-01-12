from ..core.universe import Universe
from .amberselection import amberSelection
from .protein_sel    import three2one
import pandas as pd
import numpy as np

def changeHIS(structure, out=None, mutate=[], sort_by=['chain','resid','name']):
    u = Universe(structure)
    resnames = three2one.keys()
    bA = u.atoms.resname.isin(resnames)

    proteinu = Universe(data=u.atoms[bA])
    otheru   = Universe(data=u.atoms[~bA])

    data = {'name':    [],
            'resname': [],
            'resid':   [],
            'chain':   [],
            'segname': [],
            'x':       [],
            'y':       [],
            'z':       [],
            'bfactor': 0.0}
    removeH = []

    for m in mutate:
        s = amberSelection(m[0])
        chain = s['chain']
        resid = s['resid']
        resn  = m[1]

        bA1 = chain == proteinu.atoms.chain
        bA2 = resid == proteinu.atoms.resid

        # remove HD1 / HE2 of the residue
        bA3 = proteinu.atoms['name'].isin(['HD1', 'HE2'])
        removeH.extend(list(proteinu.atoms[bA1 & bA2 & bA3].index))
            
        # check whether the residue is histidine
        # get segn and resid too
        bACA  = 'CA'  == proteinu.atoms['name']
        resnCA, segn, resid = \
            proteinu.atoms[bA1 & bA2 & bACA][['resname', 'segname', 'resid']].values[0]
        
        if resnCA not in ['HIS', 'HID', 'HIE', 'HIP', 'HSD', 'HSE', 'HSP']:
            assert 0 == 1, '%s is not histidine but %s' %(m[0], resnCA)

        # change a residue name
        proteinu.atoms.loc[bA1 & bA2, 'resname'] = resn
        
        bACE1 = 'CE1' == proteinu.atoms['name']
        posCE1 = proteinu.atoms[bA1 & bA2 & bACE1][['x','y','z']].values[0]

        if resn == 'HIE' or resn == 'HSE' or resn == 'HIP' or resn == 'HSP':
            bACD2 = 'CD2' == proteinu.atoms['name']
            bANE2 = 'NE2' == proteinu.atoms['name']

            posCD2 = proteinu.atoms[bA1 & bA2 & bACD2][['x','y','z']].values[0]
            posNE2 = proteinu.atoms[bA1 & bA2 & bANE2][['x','y','z']].values[0]
            dr = (posNE2-posCE1)/np.linalg.norm(posNE2-posCE1) + (posNE2-posCD2)/np.linalg.norm(posNE2-posCD2)
            dr /= np.linalg.norm(dr)

            posHE2 = posNE2 + dr

            data['name'].append('HE2')
            data['resname'].append(resn)
            data['segname'].append(segn)
            data['resid'].append(int(resid))
            data['chain'].append(chain)
            data['x'].append(posHE2[0])
            data['y'].append(posHE2[1])
            data['z'].append(posHE2[2])

        if resn == 'HID' or resn == 'HSD' or resn == 'HIP' or resn == 'HSP':
            bACG  = 'CG'  == proteinu.atoms['name']
            bAND1 = 'ND1' == proteinu.atoms['name']

            posCG  = proteinu.atoms[bA1 & bA2 & bACG][['x','y','z']].values[0]
            posND1 = proteinu.atoms[bA1 & bA2 & bAND1][['x','y','z']].values[0]
            dr = (posND1-posCE1)/np.linalg.norm(posND1-posCE1) + (posND1-posCG)/np.linalg.norm(posND1-posCG)
            dr /= np.linalg.norm(dr)

            posHD1 = posND1 + dr
            data['name'].append('HD1')
            data['resname'].append(resn)
            data['segname'].append(segn)
            data['resid'].append(int(resid))
            data['chain'].append(chain)
            data['x'].append(posHD1[0])
            data['y'].append(posHD1[1])
            data['z'].append(posHD1[2])
    
    proteinu.atoms.drop(removeH, inplace=True)
    protein = pd.concat([proteinu.atoms, pd.DataFrame(data)], ignore_index=True)
    protein.sort_values(by=sort_by, inplace=True)
    u.atoms = pd.concat([protein, otheru.atoms], ignore_index=True)

    if out: u.write(out)
    return u



from ..core.universe import Universe
from .protein_sel import *
from .util import completeTetra
import pandas as pd
import numpy as np

def capTerminiResidues(structure, out=None, cap_termini={}):
    u = Universe(structure)
    resnames = three2one.keys()
    bA = u.atoms.resname.isin(resnames)

    # exclude H2, H3, HT2, HT3, OT2, OXT2
    termatoms = ['H2', 'H3', 'HT2', 'HT3', 'OT2', 'OXT2']
    proteinu  = Universe(data=u.atoms[bA])
    bA_prot   = proteinu.atoms['name'].isin(termatoms)
    proteinu  = Universe(data=proteinu.atoms[~bA_prot])

    # exclude ACE, NMA
    otheru  = Universe(data=u.atoms[~bA])
    bA_term = otheru.atoms.resname.isin(['ACE', 'NMA', 'NME'])
    otheru  = Universe(data=otheru.atoms[~bA_term])

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

    
    removeProHT1 = []
    for chainnum, chain in enumerate(chains):
        try:
            Ncap, Ccap = cap_termini[chainnum]
        except:
            Ncap = 'NH3'
            Ccap = 'COO'

        if Ncap != 'NH3' and Ncap != 'ACE':
            print('Available N-terminus: NH3 or ACE... Using NH3')
            Ncap = 'NH3'

        if Ccap != 'COO' and Ccap != 'NMA':
            print('Available C-terminus: COO or NMA... Using COO')
            Ccap = 'COO'
        
        print(f'N-terminus of chain {chainnum}: {Ncap}')
        print(f'C-terminus of chain {chainnum}: {Ccap}')


        bA1 = proteinu.atoms.chain == chain
        Ntermresid = proteinu.atoms[bA1]['resid'].min()
        Ctermresid = proteinu.atoms[bA1]['resid'].max()

        # Ntermini
        bA2  = proteinu.atoms.resid == Ntermresid
        bA3  = proteinu.atoms['name'] == 'N'
        bA4  = proteinu.atoms['name'].isin(['HT1', 'HN', 'H', 'H1'])
        bA5  = proteinu.atoms['name'] == 'CA'

        resn, segn, Nx, Ny, Nz = proteinu.atoms[bA1 & bA2 & bA3][['resname', 'segname', 'x', 'y', 'z']].values[0]
        Cx, Cy, Cz = proteinu.atoms[bA1 & bA2 & bA5][['x','y','z']].values[0]
        dr1 = np.array([Cx, Cy, Cz]) - np.array([Nx, Ny, Nz])
        dr1 = dr1 / np.linalg.norm(dr1)

        # For proline, it is important to know whether there is HN.
        # if there is HN: you should add only HT2
        # if there is no HN (e.g., ACE): you should add HT1 and HT2 (adding HT2 and HT3 results in error...)
        
        Hatoms = proteinu.atoms[bA1 & bA2 & bA4]
        if len(Hatoms) != 0:
            Hx, Hy, Hz = Hatoms[['x','y','z']].values[0]
            dr2 = np.array([Hx, Hy, Hz]) - np.array([Nx, Ny, Nz])
            dr2 = dr2 / np.linalg.norm(dr2)
        
        if Ncap == 'NH3':
            if resn != 'PRO':
                dr3, dr4 = completeTetra(dr1, dr2)
                dr3 *= 1.1
                dr4 *= 1.1
                data['name'].extend(['HT2', 'HT3'])
                data['resname'].extend([resn, resn])
                data['segname'].extend([segn, segn])
                data['resid'].extend([Ntermresid, Ntermresid])
                data['chain'].extend([chain, chain])
                data['x'].extend([Nx + dr3[0], Nx + dr4[0]])
                data['y'].extend([Ny + dr3[1], Ny + dr4[1]])
                data['z'].extend([Nz + dr3[2], Nz + dr4[2]])

            elif resn == 'PRO' and len(Hatoms) > 0:
                # proline has HN
                bA6 = proteinu.atoms['name'] == 'CD'
                CDx, CDy, CDz = proteinu.atoms[bA1 & bA2 & bA6][['x','y', 'z']].values[0]
                dr2_CD = np.array([CDx, CDy, CDz]) - np.array([Nx, Ny, Nz])
                dr3, dr4 = completeTetra(dr1, dr2_CD)
                dr3 *= 1.1
                dr4 *= 1.1
                data['name'].extend(['HT2'])
                data['resname'].extend([resn])
                data['segname'].extend([segn])
                data['resid'].extend([Ntermresid])
                data['chain'].extend([chain])
                data['x'].extend([Nx + dr3[0]])
                data['y'].extend([Ny + dr3[1]])
                data['z'].extend([Nz + dr3[2]])
                
                # change HN location
                proteinu.atoms.loc[bA1 & bA2 & bA4, 'x'] = Nx + dr4[0]
                proteinu.atoms.loc[bA1 & bA2 & bA4, 'y'] = Ny + dr4[1]
                proteinu.atoms.loc[bA1 & bA2 & bA4, 'z'] = Nz + dr4[2]


            elif resn == 'PRO' and len(Hatoms) == 0:
                # proline does not have HN
                bA6 = proteinu.atoms['name'] == 'CD'
                CDx, CDy, CDz = proteinu.atoms[bA1 & bA2 & bA6][['x','y', 'z']].values[0]
                dr2_CD = np.array([CDx, CDy, CDz]) - np.array([Nx, Ny, Nz])
                dr3, dr4 = completeTetra(dr1, dr2_CD)
                dr3 *= 1.1
                dr4 *= 1.1
                data['name'].extend(['HT1', 'HT2'])
                data['resname'].extend([resn, resn])
                data['segname'].extend([segn, segn])
                data['resid'].extend([Ntermresid, Ntermresid])
                data['chain'].extend([chain, chain])
                data['x'].extend([Nx + dr3[0], Nx + dr4[0]])
                data['y'].extend([Ny + dr3[1], Ny + dr4[1]])
                data['z'].extend([Nz + dr3[2], Nz + dr4[2]])


        elif Ncap == 'ACE':
            if resn == 'PRO':
                bA6 = proteinu.atoms['name'] == 'CD'
                CDx, CDy, CDz = proteinu.atoms[bA1 & bA2 & bA6][['x','y', 'z']].values[0]
                dr2 = np.array([CDx, CDy, CDz]) - np.array([Nx, Ny, Nz])
                dr2 /= np.linalg.norm(dr2)
                
                try:
                    bA7 = proteinu.atoms['name'].isin(['HT1', 'HN', 'H'])
                    removeProHT1.append(proteinu.atoms[bA1 & bA2 & bA7].index[0])
                except:
                    pass
 
            posCY  = np.array([Nx, Ny, Nz]) - (dr1 + dr2) / np.linalg.norm(dr1 + dr2) * 1.34
            posOY  = posCY - 1.23 * dr2
            posCAY = posCY - dr1 * 1.53
            posHY1 = posCAY - (dr1 + dr2) / np.linalg.norm(dr1 + dr2) * 1.1
            posHY2, posHY3 = completeTetra(posCY - posCAY, posHY1 - posCAY)
            posHY2 = posCAY + posHY2 * 1.1
            posHY3 = posCAY + posHY3 * 1.1
            
            data['name'].extend(['CY', 'OY', 'CAY', 'HY1', 'HY2', 'HY3'])
            data['resname'].extend(['ACE'] * 6)
            data['segname'].extend([segn] * 6)
            data['resid'].extend([Ntermresid - 1] * 6)
            data['chain'].extend([chain] * 6)
            data['x'].extend([posCY[0], posOY[0], posCAY[0], posHY1[0], posHY2[0], posHY3[0]])
            data['y'].extend([posCY[1], posOY[1], posCAY[1], posHY1[1], posHY2[1], posHY3[1]])
            data['z'].extend([posCY[2], posOY[2], posCAY[2], posHY1[2], posHY2[2], posHY3[2]])
            

        # Ctermini
        bA2  = proteinu.atoms.resid == Ctermresid
        bA3  = proteinu.atoms['name'] == 'C'
        bA4  = proteinu.atoms['name'].isin(['OT1', 'O', 'O1'])
        bA5  = proteinu.atoms['name'] == 'CA'
        bA6  = proteinu.atoms['name'] == 'N'

        resn, segn, Cx, Cy, Cz = proteinu.atoms[bA1 & bA2 & bA3][['resname', 'segname', 'x', 'y', 'z']].values[0]
        Ox, Oy, Oz = proteinu.atoms[bA1 & bA2 & bA4][['x','y','z']].values[0]
        CAx, CAy, CAz = proteinu.atoms[bA1 & bA2 & bA5][['x','y','z']].values[0]
        Nx, Ny, Nz = proteinu.atoms[bA1 & bA2 & bA6][['x','y','z']].values[0]

        dr1 = np.array([CAx, CAy, CAz]) - np.array([Nx, Ny, Nz])
        dr2 = np.array([Ox, Oy, Oz]) - np.array([Cx, Cy, Cz])
        dr3 = np.array([Cx, Cy, Cz]) - np.array([CAx, CAy, CAz])

        dr1 /= np.linalg.norm(dr1)
        dr2 /= np.linalg.norm(dr2)
        dr3 /= np.linalg.norm(dr3)

        if Ccap == 'COO':
            posOT2 = np.array([Cx, Cy, Cz]) + (dr3 - dr2) / np.linalg.norm(dr3 - dr2) * 1.25
            data['name'].extend(['OT2'])
            data['resname'].extend([resn])
            data['segname'].extend([segn])
            data['resid'].extend([Ctermresid])
            data['chain'].extend([chain])
            data['x'].extend([posOT2[0]])
            data['y'].extend([posOT2[1]])
            data['z'].extend([posOT2[2]])

        elif Ccap == 'NMA':
            posN   = np.array([Cx, Cy, Cz]) + (dr3 - dr2) / np.linalg.norm(dr3 - dr2) * 1.34
            posH   = posN - dr2
            posCA  = posN + dr3 * 1.46
            posHA1 = posCA + dr2 * 1.1
            posHA2, posHA3 = completeTetra(posHA1-posCA, posN-posCA)
            posHA2 *= 1.1
            posHA3 *= 1.1
            posHA2 += posCA
            posHA3 += posCA

            data['name'].extend(['N', 'H', 'CH3', 'HH31', 'HH32', 'HH33'])
            data['resname'].extend(['NMA'] * 6)
            data['segname'].extend([segn] * 6)
            data['resid'].extend([Ctermresid + 1] * 6)
            data['chain'].extend([chain] * 6)
            data['x'].extend([posN[0], posH[0], posCA[0], posHA1[0], posHA2[0], posHA3[0]])
            data['y'].extend([posN[1], posH[1], posCA[1], posHA1[1], posHA2[1], posHA3[1]])
            data['z'].extend([posN[2], posH[2], posCA[2], posHA1[2], posHA2[2], posHA3[2]])

    
    if len(data['resid']) == 0:
        if out: u.write(out)
        return u

    proteinu.atoms.drop(removeProHT1, inplace=True)
    protein = pd.concat([proteinu.atoms, pd.DataFrame(data)], ignore_index=True)
    protein.sort_values(by=['chain', 'resid', 'name'], inplace=True)
    u.atoms = pd.concat([protein, otheru.atoms], ignore_index=True)

    if out: u.write(out)
    return u


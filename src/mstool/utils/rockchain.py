import numpy as np
import pandas as pd
from   ..core.universe import Universe

class RockChain:
    def __init__(self, structure, out='ROCK', resname='ROCK',
                 rockHtype='HAL3', rockCtype='CTL3'):
        
        self.dms = out + '.dms'
        u = Universe(structure)
        u.atoms.resid   = 1
        u.atoms.mass    = 100.0
        u.atoms.bfactor = 1.0
        u.atoms.charge  = 0.0
        bAH = u.atoms['name'].str.startswith('H')
        u.atoms.loc[bAH,  'type'] = rockHtype
        u.atoms.loc[~bAH, 'type'] = rockCtype
        chains = pd.unique(u.atoms.chain)

        ### Save it to XML
        self.xml = out + '.xml'
        W = open(self.xml, 'w')
        W.write('<ForceField>\n')
        W.write('  <Residues>\n')

        ### iter by chain
        for i, chain in enumerate(chains):
            bA = u.atoms.chain == chain
            resn = resname + str(i)
            u.atoms.loc[bA, 'resname'] = resn
            # u.atoms.loc[bA, 'chain'] = resn

            names = []
            for j, (idx, atom) in enumerate(u.atoms[bA].iterrows()):
                names.append(atom['type'][0] + str(j))
            u.atoms.loc[bA, 'name'] = names

            W.write(f'    <Residue name="{resn}">\n')
            for nn, tt in zip(u.atoms[bA].name.tolist(), u.atoms[bA].type.tolist()):
                W.write(f'      <Atom charge="0.00" name="{nn}" type="{tt}"/>\n')
            W.write(f'    </Residue>\n')

        W.write('  </Residues>\n')
        W.write('</ForceField>\n')
        W.close()

        u.write(self.dms)
        self.bonds = []


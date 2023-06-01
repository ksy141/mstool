from .utils.universe import Universe
import numpy as np
import pandas as pd

class FillSide:
    def __init__(self, mapping, universe):
        '''
        Place SC* atoms near BB atoms if they do not exist.
        '''

        dg = {'name': [], 'resid': [], 'resname': [],'chain':[],
              'x': [],'y': [],'z': [], 'bfactor': 0.0}

        records = []
        for index, atom in universe.atoms.iterrows():
            record = atom.chain + str(atom.resid) + atom.resname
            if record in records: continue

            ### Mapping was all good for this residue
            mm  = mapping.RESI[atom.resname]['CGAtoms']
            bA1 = universe.atoms.chain   == atom.chain
            bA2 = universe.atoms.resid   == atom.resid
            bA3 = universe.atoms.resname == atom.resname
            df  = universe.atoms[bA1 & bA2 & bA3]
            if len(df) == len(mm) and set(df.name) == set(mm.keys()):
                records.append(record)
                continue

            ### Some atoms are not present
            ref = df.iloc[-1][['x','y','z']].to_numpy()
            for cgname in mm.keys():
                if len(df[ df.name == cgname ]) == 1:
                    continue
                elif len(df[ df.name == cgname ]) == 0:
                    dg['name'].append(cgname)
                    dg['resid'].append(atom.resid)
                    dg['resname'].append(atom.resname)
                    dg['chain'].append(atom.chain)
                    pos = ref + np.random.rand(3) - 0.5
                    dg['x'].append(pos[0])
                    dg['y'].append(pos[1])
                    dg['z'].append(pos[2])
                    records.append(record)

                else:
                    assert 0 == 1, '# of cg atom > 1'

        newu = Universe(data=dg)
        self.universe = universe
        self.universe.atoms = pd.concat([self.universe.atoms, newu.atoms])
        self.universe.sort()

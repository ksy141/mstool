import pandas as pd
import numpy  as np
from   .utils.universe import Universe
from   .utils.seq      import Seq
#from   .protein_sel    import *

class Fill(Universe):
    def __init__(self, mapping, universe, sequence):

        self.universe = universe
        self.mapping  = mapping
        self.sequence = sequence

        # python 3.7+
        self.chains = list(dict.fromkeys(self.universe.atoms.chain))
        assert isinstance(self.universe, Universe), "universe is not Universe instance"
        assert isinstance(self.sequence, Seq),      "sequence is not Seq instance"

        self.seqfromU = Seq(universe=self.universe)
        self.seqfromF = self.sequence

        assert len(self.seqfromU.seq.keys()) == len(self.seqfromF.seq.keys()), 'len(chains)'
        assert len(self.chains) == len(self.seqfromF.seq.keys()), 'len(chains)'
        self.nchains = len(self.seqfromF.seq.keys())

        ### Add missing loops
        self.fillloops()

        ### Add / delete atoms
        self.fillside()

        super().__init__(data=self.universe.atoms)
        self.sort()
        self.cell = universe.cell
        self.dimensions = universe.dimensions


    def fillloops(self):
        '''
        find a missing loop region, fill BB atoms
        '''
        for chainnum in range(self.nchains):
            s = self.seqfromU.seq[chainnum]['one']
            missing = np.where(np.array([*s]) == '-')[0]
            grouped = self.consecutive(missing)

            # disregard the missing N-termini residues
            if grouped[0][0] == 0 and len(grouped) != 1:
                grouped = grouped[1:]

            elif grouped[0][0] == 0 and len(grouped) == 1:
                print('there is nothing to fill except for the N-termini')

            for group in grouped:
                self.fillloop(chainnum, group)


    def consecutive(self, data, stepsize=1):
        return np.split(data, np.where(np.diff(data) != stepsize)[0]+1)


    def fillloop(self, chainnum, group):
        dg = {'name': [], 'resid': [], 'resname': [],'chain':[],
              'x': [],'y': [],'z': [], 'bfactor': 0.0}

        N = len(group)
        chain = self.chains[chainnum]

        bA1 = self.universe.atoms.chain  == chain
        bA2 = (self.universe.atoms['name']  == 'CA') | (self.universe.atoms['name'] == 'BB')
        bA3 = self.universe.atoms.resid == group[0] - 1 + 1
        bA4 = self.universe.atoms.resid == group[-1] + 1 + 1

        start = self.universe.atoms[bA1 & bA2 & bA3][['x', 'y', 'z']].to_numpy()[0]
        end   = self.universe.atoms[bA1 & bA2 & bA4][['x', 'y', 'z']].to_numpy()[0]
        dr    = end - start
        step  = np.linspace(0, 1, N+2)

        for i in range(len(group)):
            resid = group[i] + 1
            posv  = start + step[i+1] * dr
            resname = self.sequence.resid(resid=resid, chainnum=chainnum)[1]

            # names = self.mapping.RESI[resname]['CGAtoms'].keys()
            # for name in names:
            #     dg['name'].append(name)
            #     dg['resid'].append(resid)
            #     dg['resname'].append(resname)
            #     dg['chain'].append(chain)
            #     pos = posv + np.random.rand(3) - 0.5
            #     dg['x'].append(pos[0])
            #     dg['y'].append(pos[1])
            #     dg['z'].append(pos[2])

            # fill BB first -> fillside will fix it
            dg['name'].append('BB')
            dg['resid'].append(resid)
            dg['resname'].append(resname)
            dg['chain'].append(chain)
            pos = posv + np.random.rand(3) - 0.5
            dg['x'].append(pos[0])
            dg['y'].append(pos[1])
            dg['z'].append(pos[2])

        newu = Universe(data=dg)
        self.universe.atoms = pd.concat([self.universe.atoms, newu.atoms], ignore_index=True)




    def fillside(self):
        '''
        Place SC* atoms near BB atoms if they do not exist
        '''
        dg = {'name': [], 'resid': [], 'resname': [], 'chain': [],
              'x': [], 'y': [], 'z': [], 'bfactor': 0.0}

        records = []
        for index, atom in self.universe.atoms.iterrows():
            record = atom.chain + str(atom.resid) + atom.resname
            if record in records: continue

            ### Read mapping scheme of the residue
            mm  = self.mapping.RESI[atom.resname]['CGAtoms']

            ### Remove extra atoms (from mutations)
            if atom['name'] not in mm.keys():
                print('dropping ' + atom['name'] + ' ' + record)
                self.universe.atoms.drop(index, inplace=True)

            ### Select all atoms of the residue            
            bA1 = self.universe.atoms.chain   == atom.chain
            bA2 = self.universe.atoms.resid   == atom.resid
            bA3 = self.universe.atoms.resname == atom.resname
            df  = self.universe.atoms[bA1 & bA2 & bA3]

            ### Mapping was all good for this residue
            if len(df) == len(mm) and set(df['name']) == set(mm.keys()):
                records.append(record)
                continue

            ### Some atoms are not present
            ref = df.iloc[-1][['x','y','z']].to_numpy()
            for cgname in mm.keys():
                if len(df[ df['name'] == cgname ]) == 1:
                    continue
                elif len(df[ df['name'] == cgname ]) == 0:
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
        self.universe.atoms = pd.concat([self.universe.atoms, newu.atoms], ignore_index=True)


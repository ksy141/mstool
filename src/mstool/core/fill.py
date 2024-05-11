import pandas as pd
import numpy  as np
from   .universe     import Universe
from   .readmappings import ReadMappings
from   .seq          import Seq

class Fill(Universe):
    def __init__(self, structure, sequence, out=None,
                 extend_termini=[], 
                 mapping=[], mapping_add={}):

        self.u              = Universe(structure)
        self.mapping        = ReadMappings(mapping, mapping_add)
        self.seqfromU       = Seq(structure=structure)
        self.seqfromF       = Seq(fasta=sequence)
        self.extend_termini = extend_termini

        # make a bfactor = 1 for BB atoms that exist in the structure
        # the rest of them; set to 0
        self.u.atoms.bfactor = 0.0
        bA = ((self.u.atoms.name == 'BB') | (self.u.atoms.name == 'SC1'))
        self.u.atoms.loc[bA, 'bfactor'] = 1.0

        # python 3.7+
        self.chains = list(dict.fromkeys(self.u.atoms.chain))

        assert len(self.seqfromU.seq.keys()) == len(self.seqfromF.seq.keys()), 'len(chains)'
        assert len(self.chains) == len(self.seqfromF.seq.keys()), 'len(chains)'
        self.nchains = len(self.seqfromF.seq.keys())

        ### Add missing loops
        self.FillLoops()

        ### Extend termini
        if len(extend_termini.keys()) != 0:
            self.ExtendTermini()


        super().__init__(data=self.u.atoms)
        self.sort()
        self.cell       = self.u.cell
        self.dimensions = self.u.dimensions

        ### save
        if out: self.write(out)


    def ExtendTermini(self):
        for chainnum in range(self.nchains):
            if chainnum not in self.extend_termini.keys(): continue
            sU = self.seqfromU.seq[chainnum]['one']
            sF = self.seqfromF.seq[chainnum]['one']

            # ---ASEF---
            # first_resid = 4
            # last_resid  = 7
            for i, oneletter in enumerate([*sU]):
                if oneletter != '-':
                    first_resid = i + 1
                    break
            last_resid  = len(sU)

            Nterm, Cterm = self.extend_termini[chainnum]

            ### N-termini
            # first_resid - 1 is the index of the first resid
            if Nterm == 'all':
                grouped = np.arange(0, first_resid - 1)

            elif Nterm > first_resid - 1.5:
                grouped = np.arange(0, first_resid - 1)

            else:
                grouped = np.arange(first_resid - 1 - Nterm, first_resid - 1)

            # Nterm = 0 -> grouped = []
            if len(grouped) != 0:
                self.fillloop(chainnum, grouped)


            ### C-termini
            # last_resid - 1 is the index of the last resid
            # last_resid - 1 + 1 is the index of the first missing resid
            if Cterm == 'all':
                grouped = np.arange(last_resid, len(sF))

            elif last_resid + Cterm >= len(sF):
                grouped = np.arange(last_resid, len(sF))

            else:
                grouped = np.arange(last_resid, last_resid + Cterm)

            # Cterm = 0 -> grouped = []
            if len(grouped) != 0:
                self.fillloop(chainnum, grouped)



    def FillLoops(self):
        '''
        find a missing loop region, fill BB atoms
        '''
        for chainnum in range(self.nchains):
            s = self.seqfromU.seq[chainnum]['one']
            missing = np.where(np.array([*s]) == '-')[0]
            grouped = self.consecutive(missing)

            # complete protein structure - skip
            try:
                grouped[0][0]
            except:
                continue

            # disregard the missing N-termini residues
            if grouped[0][0] == 0 and len(grouped) != 1:
                grouped = grouped[1:]

            elif grouped[0][0] == 0 and len(grouped) == 1:
                print(f'chain{chainnum}: there is nothing to fill')
                continue

            for group in grouped:
                self.fillloop(chainnum, group)


    def consecutive(self, data, stepsize=1):
        '''
        split an array so that the difference between the values of each array is stepsize.
        consecutive([136, 137, 138, 142, 143, 145, 146]) -> 
        [array([136, 137, 138]), array([142, 143]), array([145, 146])]
        '''
        return np.split(data, np.where(np.diff(data) != stepsize)[0] + 1)


    def fillloop(self, chainnum, group):
        # s = '---AB---DEFG'
        # missing = np.where(np.array([*s]) == '-')[0]
        # grouped = np.split(missing, np.where(np.diff(missing) != 1)[0] + 1)
        # grouped -> [array([0, 1, 2]), array([5, 6, 7])]
        # grouped is the indices (0-indexed) of the "missing"

        dg = {'name': [], 'resid': [], 'resname': [],'chain':[],
              'x': [],'y': [],'z': [], 'bfactor': 0.0}

        N = len(group)
        chain = self.chains[chainnum]

        bA1 = self.u.atoms.chain  == chain
        bA2 = (self.u.atoms['name']  == 'CA') | (self.u.atoms['name'] == 'BB')
        bA3 = self.u.atoms.resid == (group[0]  - 1) + 1
        bA4 = self.u.atoms.resid == (group[-1] + 1) + 1

        start = self.u.atoms[bA1 & bA2 & bA3][['x', 'y', 'z']].to_numpy()
        end   = self.u.atoms[bA1 & bA2 & bA4][['x', 'y', 'z']].to_numpy()

        step = np.linspace(0, 1, N+2)
        
        if len(start) == 0 and len(end) == 0:
            assert 0 == 1, 'review grouped'

        elif len(start) == 1 and len(end) == 0:
            # C-termini
            start = start[0]
            dr    = 0.0

        elif len(start) == 0 and len(end) == 1:
            # N-termini
            start = end[0]
            dr    = 0.0

        elif len(start) == 1 and len(end) == 1:
            # loop has both ends
            start = start[0]
            end   = end[0]
            dr    = end - start

        else:
            assert 0 == 1, 'review grouped'


        for i in range(len(group)):
            resid = group[i] + 1
            posv  = start + step[i+1] * dr
            resname = self.seqfromF.resid(resid=resid, chainnum=chainnum)[1]

            names = self.mapping.RESI[resname]['CGAtoms'].keys()
            for name in names:
                dg['name'].append(name)
                dg['resid'].append(resid)
                dg['resname'].append(resname)
                dg['chain'].append(chain)
                pos = posv + np.random.rand(3) - 0.5
                dg['x'].append(pos[0])
                dg['y'].append(pos[1])
                dg['z'].append(pos[2])

        newu = Universe(data=dg)
        self.u.atoms = pd.concat([self.u.atoms, newu.atoms], ignore_index=True)



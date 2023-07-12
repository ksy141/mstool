from .universe import Universe
from .seq      import Seq

class Mutate(Universe):
    def __init__(self, structure, fasta, out=None):

        self.universe = Universe(structure)
        self.seqfromU = Seq(structure=structure)
        self.seqfromF = Seq(fasta)
        self.atoms    = self.universe.atoms

        self.fixmutation()
        df = self.universe.atoms[ self.atoms.bfactor > 0.5 ]

        super().__init__(data = df)
        self.sort()
        self.cell       = self.universe.cell
        self.dimensions = self.universe.dimensions

        ### save
        if out is not None:
            self.write(out)

    def fixmutation(self):
        '''
        fix mutation (at atomistic level)
        more like save backbone + CB atoms
        '''
        
        for key in self.seqfromU.seq.keys():

            if key not in self.seqfromF.seq.keys():
                continue

            U = self.seqfromU.seq[key]                
            F = self.seqfromF.seq[key]

            Fneeded = F['three'][:len(U['three'])]
            Uneeded = U['three']
            Uchain  = U['chain']

            for i, (Faa, Uaa) in enumerate(zip(Fneeded, Uneeded), 1):

                if Uaa == '---':
                    continue

                if Uaa == Faa:
                    continue

                if Faa in ['HIS', 'HIE', 'HIP', 'HID', 'HSE', 'HSD', 'HSP']:
                    Faa = 'HIS'

                if Faa == 'HIS' and (Uaa in ['HIS', 'HIE', 'HIP', 'HID', 'HSE', 'HSD', 'HSP']):
                    continue
              
                print(f'Mutating /{Uchain}:{i}:{Uaa} to /{Uchain}:{i}:{Faa}')
                bA1 = self.atoms.chain == self.seqfromU.seq[key]['chain']
                bA2 = self.atoms.resid == i
                self.atoms.loc[bA1 & bA2, 'resname'] = Faa

                if Faa == 'GLY':
                    bA3 = self.atoms['name'].isin(['N', 'CA', 'C', 'O', 'HA1', 'HA'])

                else:
                    bA3 = self.atoms['name'].isin(['N', 'CA', 'C', 'O', 'HA1', 'HA', 'CB'])
                
                # backbone
                self.atoms.loc[bA1 & bA2 & bA3, 'bfactor'] = 1.0
                
                # side chain
                self.atoms.loc[bA1 & bA2 & (~bA3), 'bfactor'] = 0.0




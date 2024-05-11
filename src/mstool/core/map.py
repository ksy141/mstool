import pandas as pd
import numpy  as np
from   .universe     import Universe
from   .readmappings import ReadMappings

class Map(Universe):
    def __init__(self, structure, out=None, mapping=[], mapping_add=[], BB2CA=True, translate=[0.0, 0.0, 0.0], add_notavailableAAAtoms=False):

        self.u        = Universe(structure)
        self.mapping  = ReadMappings(mapping, mapping_add)
        self.resnames = set(self.u.atoms.resname)
        self.BB2CA    = BB2CA

        # construct
        firstdata = self.add_availableCGAtoms(add_notavailableAAAtoms)
        super().__init__(data=firstdata)
        self.cell = self.u.cell
        self.dimensions = self.u.dimensions

        # add missing atoms
        seconddata = self.add_missingCGAtoms()
        df = Universe(data=seconddata)
        self.atoms = pd.concat([self.atoms, df.atoms], ignore_index=True)
        self.sort()

        # translate
        self.atoms[['x','y','z']] += translate

        # save
        if out is not None:
            self.write(out)


    def add_availableCGAtoms(self, add_notavailableAAAtoms):
        dfresids   = []; dfresnames = [];
        dfchains   = []; dfnames    = [];
        dfx = []; dfy = []; dfz = [];       

        for resname in self.resnames:

            # skip cg mapping if the residue not exists in the mapping scheme
            if resname not in self.mapping.RESI.keys():
                if add_notavailableAAAtoms:
                    for index, atom in self.u.atoms[self.u.atoms.resname == resname].iterrows():
                        dfchains.append(atom.chain)
                        dfresids.append(atom.resid)
                        dfresnames.append(resname)
                        dfnames.append(atom['name'])
                        dfx.append(atom.x)
                        dfy.append(atom.y)
                        dfz.append(atom.z)
                continue


            # map
            for cgname, aanames in self.mapping.RESI[resname]['CGAtoms'].items():
                bA1 = self.u.atoms.resname == resname

                if cgname == 'BB' and self.BB2CA:
                    bA2 = self.u.atoms.name == 'CA'
                else:
                    bA2 = self.u.atoms.name.isin(self.mapping.RESI[resname]['CGAtoms'][cgname])
                
                dg  = self.u.atoms[bA1 & bA2].groupby(['chain', 'resid']).mean(numeric_only=True)
                dr  = np.zeros(3)

                for index, atom in dg.iterrows():
                    dfchains.append(index[0])
                    dfresids.append(index[1])
                    dfresnames.append(resname)
                    dfnames.append(cgname)
                    dfx.append(atom.x)
                    dfy.append(atom.y)
                    dfz.append(atom.z)


        data = {'name': dfnames, 'resname': dfresnames,
                'resid': dfresids, 'chain': dfchains,
                'x': dfx, 'y': dfy, 'z': dfz,
                'bfactor': 1.0}

        return data


    def add_missingCGAtoms(self):
        
        # no aa atoms; so no cg atoms
        addresids = []; addresnames = [];
        addchains = []; addnames = []; addpos = [];
        saved = []

        for resname in self.resnames:
            residues = self.atoms[ self.atoms.resname == resname ]
            
            for index, atom in residues.iterrows():
                chain = atom.chain
                resid = atom.resid
                
                save  = f"{chain}_{resname}_{resid}"

                if save in saved: continue

                residue = residues[(residues.chain == chain) & (residues.resid == resid)]
                if len(residue) == len(self.mapping.RESI[resname]['CGAtoms'].keys()):
                    # this residue has all CG atoms
                    saved.append(save)
                    continue

                for cgname in self.mapping.RESI[resname]['CGAtoms'].keys():
                    if cgname not in residue['name'].tolist():
                        addresids.append(resid)
                        addresnames.append(resname)
                        addchains.append(chain)
                        addnames.append(cgname)
                        addpos.append(residue[['x','y','z']].to_numpy()[0] + np.random.rand(3) - 0.5)
                saved.append(save)
                print(f"Some CGAtoms of {save} were randomly built.")

        
        addpos    = np.array(addpos)
        addresids = np.array(addresids).astype(int)
        data   = {'name':  addnames,  'resname': addresnames, 
                  'resid': addresids, 'chain': addchains}

        if len(addpos) == 0:
            data['x'] = []
            data['y'] = []
            data['z'] = []

        else:
            data['x'] = addpos[:,0]
            data['y'] = addpos[:,1]
            data['z'] = addpos[:,2]

        return data

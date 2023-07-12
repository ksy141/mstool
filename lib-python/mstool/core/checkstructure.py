import os
import numpy as np
from  .universe           import Universe
from  .readmappings       import ReadMappings
from  ..utils.protein_sel import three2one

class CheckStructure:
    def __init__(self, structure, mapping=[], mapping_add=[], log=None):
        self.structure = structure
        self.u         = Universe(self.structure)
        self.resnames  = set(self.u.atoms.resname)
        self.mapping   = ReadMappings(mapping, mapping_add)
        self.log       = log
        self.reports   = []
        self.counts    = []

        self.u.sort()
        self.peptide()
        self.chirals()
        self.cistrans()
        self.printout()


    def chirals(self):
        for resname in self.resnames:
            if resname not in self.mapping.RESI.keys():
                print(f'{resname} is unrecognized residue. Skipping a chiral examination for this residue')
                continue

            chirals = self.mapping.RESI[resname]['chiral']
            for chiral in chirals:
                self.counts.append(['chiral', resname, chiral[0], chiral[1], chiral[2], chiral[3], chiral[4]])

                if resname in three2one.keys() and chiral[0] == 'HA':
                    bA1 = ((self.u.atoms.name == 'HA') | (self.u.atoms.name == 'HA1') | (self.u.atoms.name == 'HA2'))
                else:
                    bA1 = self.u.atoms.name == chiral[0]

                bA  = self.u.atoms.resname == resname
                bA2 = self.u.atoms.name == chiral[1]
                bA3 = self.u.atoms.name == chiral[2]
                bA4 = self.u.atoms.name == chiral[3]
                bA5 = self.u.atoms.name == chiral[4]

                # atomC -> atomD -> atomE
                target = self.u.atoms[bA & bA1]
                center = self.u.atoms[bA & bA2]
                atomC  = self.u.atoms[bA & bA3]
                atomD  = self.u.atoms[bA & bA4]
                atomE  = self.u.atoms[bA & bA5]

                assert len(target) == len(center) == \
                    len(atomC) == len(atomD) == len(atomE), \
                    "the length of atoms is different"

                BA = target[['x','y','z']].to_numpy() - center[['x','y','z']].to_numpy()
                CD = atomD[['x','y','z']].to_numpy()  - atomC[['x','y','z']].to_numpy()
                DE = atomE[['x','y','z']].to_numpy()  - atomD[['x','y','z']].to_numpy()

                v   = np.cross(CD, DE)
                cos = np.sum(BA * v, axis=-1) / np.linalg.norm(BA, axis=-1) / np.linalg.norm(v, axis=-1)
                TF  = cos < 0.0

                if np.sum(TF) != 0:
                    self.reports.append(['chiral', resname, target[TF].resid.to_numpy(), target[TF].chain.to_numpy()])


    def cistrans(self):
        for resname in self.resnames:
            if resname not in self.mapping.RESI.keys():
                print(f'{resname} is unrecognized residue. Skipping a cistrans examination for this residue')
                continue

            for isomer in ['cis', 'trans']:
                atomset = self.mapping.RESI[resname][isomer]

                for conf in atomset:
                    self.counts.append([isomer, resname, conf[0], conf[1], conf[2], conf[3]])

                    bA  = self.u.atoms.resname == resname
                    bA1 = self.u.atoms.name    == conf[0]
                    bA2 = self.u.atoms.name    == conf[1]
                    bA3 = self.u.atoms.name    == conf[2]
                    bA4 = self.u.atoms.name    == conf[3]

                    a1  = self.u.atoms[bA & bA1]
                    a2  = self.u.atoms[bA & bA2]                
                    a3  = self.u.atoms[bA & bA3]
                    a4  = self.u.atoms[bA & bA4]

                    dr1 = a2[['x','y','z']].to_numpy() - a1[['x','y','z']].to_numpy()
                    dr2 = a3[['x','y','z']].to_numpy() - a2[['x','y','z']].to_numpy()
                    dr3 = a4[['x','y','z']].to_numpy() - a3[['x','y','z']].to_numpy()

                    plane1 = np.cross(dr1, dr2)
                    plane2 = np.cross(dr2, dr3)

                    plane1 /= np.linalg.norm(plane1, axis=-1)[:, np.newaxis]
                    plane2 /= np.linalg.norm(plane2, axis=-1)[:, np.newaxis]

                    # cis:   cos > 0.0
                    # trans: cos < 0.0
                    cos = np.sum(plane1 * plane2, axis=-1)

                    # trans... but if cis, something went wrong
                    if isomer == 'trans':
                        TF = cos > 0.0
                        if np.sum(TF) != 0:
                            self.reports.append(['trans', resname, a1[TF].resid.to_numpy(), a1[TF].chain.to_numpy()])

                    if isomer == 'cis':
                        TF = cos < 0.0
                        if np.sum(TF) != 0:
                            self.reports.append(['cis', resname, a1[TF].resid.to_numpy(), a1[TF].chain.to_numpy()])


    def peptide(self):
        protein_resnames = three2one.keys()
        protein_chains   = sorted(list(set(self.u.atoms[ self.u.atoms.resname.isin(protein_resnames) ].chain)))
        for protein_chain in protein_chains:
            self.counts.append(['peptide bond cis/trans: chain ', protein_chain])

        bA  = self.u.atoms.resname.isin(protein_resnames)
        bAN = self.u.atoms.name == 'N'
        bAH = (self.u.atoms.name.isin(['HN', 'H'])) | ((self.u.atoms.resname == 'PRO') & (self.u.atoms.name == 'CD'))
        bAC = self.u.atoms.name == 'C'
        bAO = self.u.atoms.name == 'O'

        Natoms = self.u.atoms[bA & bAN]
        Hatoms = self.u.atoms[bA & bAH]
        Catoms = self.u.atoms[bA & bAC]
        Oatoms = self.u.atoms[bA & bAO]

        for i in range(len(Catoms)):
            # Catom is a series object; other atoms are DataFrame objects
            Catom = Catoms.iloc[i]
            resid = Catom.resid
            chain = Catom.chain
            Cresname = Catom.resname

            Oatom = Oatoms[ (Oatoms.resid == resid)   & (Oatoms.chain == chain)]
            Natom = Natoms[ (Natoms.resid == resid+1) & (Natoms.chain == chain)]
            Hatom = Hatoms[ (Hatoms.resid == resid+1) & (Hatoms.chain == chain)]
            if len(Oatom) * len(Natom) * len(Hatom) != 1: continue
            Nresname = Natom['resname'].iloc[0]

            name   = Hatom['name'].iloc[0]
            rCatom = Catom[['x','y','z']].to_numpy()
            rOatom = Oatom[['x','y','z']].to_numpy()
            rNatom = Natom[['x','y','z']].to_numpy()
            rHatom = Hatom[['x','y','z']].to_numpy()

            dr1 = rCatom - rOatom
            dr2 = rNatom - rCatom
            dr3 = rHatom - rNatom

            plane1 = np.cross(dr1, dr2)
            plane2 = np.cross(dr2, dr3)

            if np.sum(plane1 * plane2) > 0:
                self.reports.append(['peptide bond cis/trans:', f'(chain {chain} and name C O and resid {resid} and resname {Cresname})', 
                    f'(chain {chain} and name {name} N and resid {resid+1} and resname {Nresname})'])


    def printout(self):
        spacedivide = "###############################################\n"
        w = ''
        w += spacedivide
        w += f"{self.structure} was reviewed\n"
        w += spacedivide + "\n"


        ### Print what isomers were reviewed
        w += "The following isomers were reviewed:\n"
        for count in self.counts:
            if len(count) == 2:
                # peptide bond
                w += f"{count[0]} {count[1]}\n"
            if len(count) == 7:
                # chiral
                w += f"{count[0]}: resname {count[1]} - {count[2]} {count[3]} {count[4]} {count[5]} {count[6]}\n"
            if len(count) == 6:
                # cis/trans
                w += f"{count[0]}: resname {count[1]} - {count[2]} {count[3]} {count[4]} {count[5]}\n"
        w += spacedivide + '\n'


        ### Print results
        w += spacedivide
        if len(self.reports) == 0:
            w += 'No molecules had flipped isomers\n'

        else:
            w += 'The following molecules had flipped isomers:\n'
            for r in self.reports:
                if r[0] in ['chiral', 'cis', 'trans']:
                    resname = r[1]
                    for resid, chain in zip(r[2], r[3]):
                        w += f'{r[0]}: resname {resname} and resid {resid} and chain {chain}\n'
                else:
                    # peptide bond
                    w += f'{r[0]} {r[1]} or {r[2]}\n'

        w += spacedivide

        print(w)
        if self.log is not None:
            with open(self.log, 'w') as f:
                f.write(w)




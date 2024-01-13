import os
import numpy as np
from  .universe           import Universe
#from  .readxml           import ReadXML
from  .readmappings       import ReadMappings
from  ..utils.protein_sel import three2one


class Backmap(Universe):
    def __init__(self, structure, out=None, refstructure=None, 
        ff=[], ff_add=[],
        mapping=[], mapping_add=[],
        backbone=True,
        guess_atomic_number=False):

        self.data = {'resid': [], 'resname': [], 'chain': [], 'name': [], 'x': [], 'y': [], 'z': []}
        #self.xml           = ReadXML(ff, ff_add)
        self.mapping       = ReadMappings(mapping, mapping_add)
        self.u             = Universe(structure)
        self.prot_resnames = list(three2one.keys())
        self.bond          = []
        self.backbone      = backbone

        ### Change BB -> CA
        bA2 = self.u.atoms.resname.isin(self.prot_resnames)
        bA3 = self.u.atoms.name == 'BB'
        self.chains = set(self.u.atoms[bA2].chain)
        self.resnames = set(self.u.atoms.resname)
        self.u.atoms.loc[bA2 & bA3, 'name'] = 'CA'

        ### Construct
        if self.backbone:
            # Prebuild backbone according to cg2aa paper
            self.construct_backbone()
            self.construct()
        else:
            # do not prebuild backbone
            self.construct()

        self.data['bfactor'] = 0.0
        super().__init__(data=self.data, guess_atomic_number=guess_atomic_number)

        if refstructure != None:
            self.r          = Universe(refstructure)
            self.cell       = self.r.cell
            self.dimensions = self.r.dimensions
        else:
            self.cell       = self.u.cell
            self.dimensions = self.u.dimensions

        self.remove_duplicate_coordinates()
        self.sort()

        if out is not None:
            self.write(out)


    def remove_duplicate_coordinates(self):
        '''
        In case there are some atoms that have the exactly same coordinates
        '''
        
        xyz = self.atoms[['x','y','z']]
        uni, count = np.unique(xyz, axis=0, return_counts=True)
        dup = uni[count > 1]

        while len(dup) > 1:
            for d in dup:
                print('atoms that have the exactly same coordinates found... moving')
                x = d[0]; y = d[1]; z = d[2]
                bAx = self.u.atoms['x'] == x
                bAy = self.u.atoms['y'] == y
                bAz = self.u.atoms['z'] == z

                natoms = self.u.atoms[bAx & bAy & bAz]
                self.u.atoms.loc[bAx & bAy & bAz, ['x','y','z']] = np.array([x,y,z]) + np.random.rand(natoms, 3)

                xyz = self.atoms[['x','y','z']]
                uni, count = np.unique(xyz, axis=0, return_counts=True)
                dup = uni[count > 1]  


    def construct_backbone(self):
        '''
        Construct a backbone according to the following CG2AA paper.
        CG2AA: backmapping protein coarse-grained structures
        Leandro E. Lombardi, Marcelo A. Marti, and Luciana Capece
        https://www.ic.fcen.uba.ar/cg2aa/cg2aa.py

        ALA and GLY are the residues that do not have SC1.
        Therefore, the software builds GLY's HA2 and ALA's CB+HB1+HB2+HB3

        PRO is the residue that does not have HN.
        '''

        dcc   = 1.522
        dco   = 1.229
        dcn   = 1.335
        dhn   = 1.100

        cosg1 = 0.935 #g1 = 20.7 CA CA CO
        sing1 = 0.353
        cosg2 = 0.969 #g2 = 14.2 CA CA N
        sing2 = 0.245
        cos1166 = 0.44775908783 #CCN
        sin1166 = 0.89415423683 


        for chain in self.chains:
            bA1 = self.u.atoms.chain == chain
            bA2 = self.u.atoms.resname.isin(self.prot_resnames)
            bA3 = self.u.atoms.name  == 'CA'

            BB  = self.u.atoms[bA1 & bA2 & bA3]
            #assert (BB.resid.to_numpy()[1:] - BB.resid.to_numpy()[:-1] == 1).all(), 'not continuous?'
            
            resname = BB['resname'].tolist()
            resid   = BB['resid'].tolist()
            r       = BB[['x','y','z']].to_numpy()

            ### Add BB = CA atoms
            self.list2data(name    = ['CA'] * len(BB),
                           chain   = [chain] * len(BB),
                           resname = BB.resname.tolist(),
                           resid   = BB.resid.tolist(),
                           pos     = BB[['x','y','z']])

            ### BACKBONE
            ca1 = r[1:-1] - r[:-2]
            ca2 = r[2:]   - r[:-2]
            p_vector = np.cross(ca1, ca2)
            e1  = ca1 / np.linalg.norm(ca1, axis=1)[:,None]
            e2  = p_vector / np.linalg.norm(p_vector, axis=1)[:,None]


            ### CA's hydrogen: HA (N->CB->C)
            tca1 = r[1:] - r[:-1]
            te1  = tca1 / np.linalg.norm(tca1, axis=1)[:,None]
            w1, w2 = self.construct_tetra(te1[:-1], te1[1:])

            for i, (index, atom) in enumerate(BB[1:-1].iterrows()):
                if atom.resname == 'GLY':
                    self.list2data(name    = ['HA1', 'HA2'],
                                   chain   = [chain]        * 2,
                                   resname = [atom.resname] * 2,
                                   resid   = [atom.resid]   * 2,
                                   pos     = [r[i+1] + w1[i], r[i+1] + w2[i]])

                elif atom.resname == 'ALA':
                    posHBs = np.tile(r[i+1] + dcc * w2[i], (3,1)) + np.random.rand(3,3) - 0.5
                    pos    = [r[i+1] + w1[i], r[i+1] + dcc * w2[i], posHBs[0], posHBs[1], posHBs[2]]
                    self.list2data(name    = ['HA', 'CB', 'HB1', 'HB2', 'HB3'],
                                   chain   = [chain]        * 5, 
                                   resname = [atom.resname] * 5,
                                   resid   = [atom.resid]   * 5, 
                                   pos     = pos)


                else:
                    self.list2data(name    = ['HA'],
                                   chain   = [chain],
                                   resname = [atom.resname],
                                   resid   = [atom.resid], 
                                   pos     = [r[i+1] + w1[i]])


            ### Construct Backbone: C, O, nextN
            posC = r[:-2] + dcc * cosg1 * e1 + dcc * sing1 * e2
            posO = posC + dco * e2
            posnextN = posC + dcn * sin1166 * e1 - dcn * cos1166 * e2
            posnextHN = posnextN - dhn * e2
            pos = {'C': posC, 'O': posO, 'N': posnextN}

            resnames = [resname[:-2], resname[:-2], resname[1:-1]]
            resids   = [resid[:-2],   resid[:-2],   resid[1:-1]]

            for rn, ri, (name, p) in zip(resnames, resids, pos.items()):
                assert len(rn) == len(p)
                self.list2data(name    = [name]  * len(p),
                               chain   = [chain] * len(p),
                               resname = rn,
                               resid   = ri,
                               pos     = p)


            ### Construct HN: Need to separate it from the above as PRO does not have HN.
            for i, (index, atom) in enumerate(BB[1:-1].iterrows()):
                if atom.resname == 'PRO': continue
                self.list2data(name    = ['HN'],
                               chain   = [chain],
                               resname = [atom.resname],
                               resid   = [atom.resid],
                               pos     = [posnextHN[i]])



            ### FIRST residue's N, HN1, HN2, HN3, HA
            if resname[0] == 'GLY':
                name = ['N', 'HT1', 'HT2', 'HT3', 'HA1', 'HA2']
            elif resname[0] == 'ALA':
                name = ['N', 'HT1', 'HT2', 'HT3', 'HA', 'CB', 'HB1', 'HB2', 'HB3']
            elif resname[0] == 'PRO':
                name = ['N', 'HT1', 'HT2', 'HA']
            else:
                name = ['N', 'HT1', 'HT2', 'HT3', 'HA']

            n_atoms = len(name)
            nterpos = np.tile(r[0], (n_atoms, 1)) + np.random.rand(n_atoms, 3) - 0.5

            self.list2data(name    = name,
                           chain   = [chain]      * n_atoms,
                           resname = [resname[0]] * n_atoms, 
                           resid   = [resid[0]]   * n_atoms,
                           pos     = nterpos)
            

            ### LAST - 1 residue's C O
            dr = r[-1] - r[-2]
            tmp = np.random.rand(3)
            p_vector = np.cross(dr, tmp)
            e1 = dr / np.linalg.norm(dr)
            e2 = p_vector / np.linalg.norm(p_vector)
            posC = r[-2] + dcc * cosg1 * e1 + dcc * sing1 * e2
            posO = posC + dco * e2

            self.list2data(name    = ['C', 'O'], 
                           chain   = [chain]       * 2, 
                           resname = [resname[-2]] * 2, 
                           resid   = [resid[-2]]   * 2,
                           pos     = [posC, posO])



            ### LAST residues' N HN C OT1 OT2
            name = ['N', 'C', 'OT1', 'OT2']
            if resname[-1] == 'GLY':
                name += ['HN', 'HA1', 'HA2']
            elif resname[-1] == 'ALA':
                name += ['HN', 'HA', 'CB', 'HB1', 'HB2', 'HB3']
            elif resname[-1] == 'PRO':
                name += ['HA']
            else:
                name += ['HN', 'HA']

            posnextN    = posC + dcn * sin1166 * e1 - dcn * cos1166 * e2
            posnextC    = r[-1] + dcc * cosg1 * e1 + dcc * sing1 * e2
            posnextOT1  = posnextC + dco * e2
            posnextOT2  = r[-1] + (dcc * cosg1 * e1 + dcc * sing1 * e2) * 2

            n_atoms_add = len(name) - 4
            posadd      = np.random.rand(n_atoms_add, 3) - 0.5
            poscter     = np.concatenate([[posnextN, posnextC, posnextOT1, posnextOT2], posadd], axis=0)
            n_atoms     = len(name)

            self.list2data(name    = name,
                           chain   = [chain] * n_atoms,
                           resname = [resname[-1]] * n_atoms,
                           resid   = [resid[-1]]   * n_atoms,
                           pos     = poscter)
          


    def construct(self):
        for resname in self.resnames:
            # if resname not in self.xml.RESI.keys():
            #     if resname != 'HIS':
            #         print(f'Warning: openMM xml does not have {resname}. Skipping backmapping for this molecule.')
            #         continue
            if resname not in self.mapping.RESI.keys():
                print(f'Warning: mapping does not have {resname}. Skipping backmapping for this molecule.')
                continue

            for CGAtom, AAAtoms in self.mapping.RESI[resname]['CGAtoms'].items():
                # if you already construct backbones accroding to cg2aa, skip backbone
                if self.backbone and CGAtom == 'BB': continue
                if CGAtom == 'BB': CGAtom = 'CA'

                bA1 = self.u.atoms.resname == resname
                bA2 = self.u.atoms.name    == CGAtom
                CG  = self.u.atoms[bA1 & bA2]

                n_atoms  = len(AAAtoms)
                name     = AAAtoms * len(CG)

                chain    = np.repeat(CG['chain'].tolist(),         n_atoms)
                resname2 = np.repeat(CG['resname'].tolist(),       n_atoms)
                resid    = np.repeat(CG['resid'].tolist(),         n_atoms)
                pos      = np.repeat(CG[['x','y','z']].to_numpy(), n_atoms, axis=0)
                pos     += np.random.rand(*pos.shape) - 0.5

                self.list2data(name    = name,
                               chain   = chain,
                               resname = resname2,
                               resid   = resid,
                               pos     = pos)
                

    def construct_tetra(self, v1, v2):
        #v1 = v1 / np.linalg.norm(v1)
        #v2 = v2 / np.linalg.norm(v2)
        nv = np.cross(v1,v2)
        nv /= np.linalg.norm(nv, axis=1)[:,None]
        nv *= 0.8164965809277259
        v = (-v1-v2)/2
        w1 = v + nv
        w2 = v - nv
        return  (w1,w2)

    def list2data(self, name, chain, resname, resid, pos):
        pos = np.array(pos)
        self.data['name'].extend(name)
        self.data['chain'].extend(chain)
        self.data['resname'].extend(resname)
        self.data['resid'].extend(resid)
        self.data['x'].extend(list(pos[:,0]))
        self.data['y'].extend(list(pos[:,1]))
        self.data['z'].extend(list(pos[:,2]))


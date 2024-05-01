import os
import numpy as np
import glob
from  .universe           import Universe
from  ..lib.align         import rotation_matrix
#from  .readxml           import ReadXML
from  .readmappings       import ReadMappings
from  ..utils.protein_sel import three2one
from  ..utils.util        import fibo

pwd = os.path.dirname(os.path.realpath(__file__))
AAdefaults = glob.glob(pwd + '/../mapping/structures/*_AA.pdb')

class Ungroup(Universe):
    def __init__(self, structure, out=None, refstructure=None, 
        ff=[], ff_add=[],
        mapping=[], mapping_add=[],
        backbone=True,
        sort=True,
        guess_atomic_number=False, fibor=0.5, version='v1',
        water_resname='W', water_chain=None, water_number=4, water_fibor=2.0, water_chain_dms=False,
        use_AA_structure=False, AA_structure=[], AA_structure_add=[], AA_shrink_factor=0.7):

        self.data = {'resid': [], 'resname': [], 'chain': [], 'name': [], 'x': [], 'y': [], 'z': []}
        #self.xml           = ReadXML(ff, ff_add)
        self.mapping       = ReadMappings(mapping, mapping_add)
        self.u             = Universe(structure)
        self.prot_resnames = list(three2one.keys())
        self.bond          = []
        self.backbone      = backbone
        self.fibor         = fibor
        self.version       = version
        self.water_resname = water_resname
        self.water_chain   = water_chain
        self.water_number  = water_number
        self.water_fibor   = water_fibor
        self.water_chain_dms = water_chain_dms

        ### Change BB -> CA
        bA2 = self.u.atoms.resname.isin(self.prot_resnames)
        bA3 = self.u.atoms.name == 'BB'
        self.chains = set(self.u.atoms[bA2].chain)
        self.resnames = set(self.u.atoms.resname)
        self.u.atoms.loc[bA2 & bA3, 'name'] = 'CA'

        ### version
        if version == 'v1':
            print('using Ungroup version 1: randomize positions (stochastic)')
        elif version == 'v2':
            print('using Ungroup version 2: place atoms on sphere with a radius of {:.2f} A (deterministic)'.format(fibor))
        else:
            assert 0 == 1, 'provide either version = "v1" or version = "v2"'

        ### use_AA
        self.exclude_residues = []
        self.AA_shrink_factor = AA_shrink_factor
        if use_AA_structure:
            if not isinstance(AA_structure, list):
                AA_structure = [AA_structure]

            if not isinstance(AA_structure_add, list):
                AA_structure_add = [AA_structure_add]

            if AA_structure == []:
                AA_structure = AAdefaults

            self.AA_structure = [*AA_structure, *AA_structure_add]
            self.construct_using_AA()


        ### Construct
        if self.backbone:
            # Prebuild backbone according to cg2aa paper
            self.construct_backbone()
            self.construct()
        else:
            # do not prebuild backbone
            self.construct()

        ### Construct water
        self.construct_water()

        self.data['bfactor'] = 0.0
        super().__init__(data=self.data, guess_atomic_number=guess_atomic_number)

        if refstructure != None:
            self.r          = Universe(refstructure)
            self.cell       = self.r.cell
            self.dimensions = self.r.dimensions
        else:
            self.cell       = self.u.cell
            self.dimensions = self.u.dimensions
        
        #if self.version == 'v1': self.remove_duplicate_coordinates()
        if sort: self.sort()
        if out is not None: self.write(out)


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


    def construct_using_AA(self):
        for ifile in self.AA_structure:
            basename = os.path.basename(ifile).split('.')[0]
            ext      = os.path.basename(ifile).split('.')[-1]
            prefix   = '/'.join(os.path.abspath(ifile).split('/')[:-1])
            resname  = basename.split('_')[0]
            path     = prefix + f'/{resname}.{ext}'

            if not os.path.exists(path):
                print(f'{path} does not exist')
                continue

            print('Using AA structure: ' + path)
            self.exclude_residues.append(resname)
            
            bA1 = self.u.atoms.resname == resname
            if bA1.sum() == 0: continue

            resns = self.u.atoms[bA1]['resn'].unique()
            for resn in resns:
                bA2 = self.u.atoms.resn == resn
                refatoms = self.u.atoms[bA1 & bA2]
                resid    = refatoms['resid'].values[0]
                chain    = refatoms['chain'].values[0]

                refpos   = self.u.atoms[bA1 & bA2][['x','y','z']].values
                refcog   = np.average(refpos, axis=0)

                mobatoms = Universe(path).atoms
                mobpos   = mobatoms[['x','y','z']].values
                mobcog   = np.average(mobpos, axis=0)

                R, min_rmsd = rotation_matrix(mobpos - mobcog, refpos - refcog)
                
                aaatoms  = Universe(ifile).atoms
                aapos    = aaatoms[['x','y','z']].values
                aacog    = np.average(aapos, axis=0)
                pos      = refcog + self.AA_shrink_factor * (R @ (aapos - aacog).T).T 
                n_atoms  = len(aaatoms)

                self.list2data(name    = aaatoms['name'].values,
                               chain   = [chain]   * n_atoms,
                               resname = [resname] * n_atoms,
                               resid   = [resid]   * n_atoms,
                               pos     = pos)
     

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
            posadd      = posnextN + np.random.rand(n_atoms_add, 3) - 0.5
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
                if resname != self.water_resname:
                    print(f'Warning: mapping does not have {resname}. Skipping backmapping for this molecule.')
                continue

            if resname in self.exclude_residues:
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

                if self.version == 'v1':
                    kick = np.random.rand(*pos.shape) - 0.5

                if self.version == 'v2':
                    fibopos = fibo(r=self.fibor, N=n_atoms, verbose=False, plot=None)
                    fibopos -= np.average(fibopos, axis=0)
                    fibopos *= np.linspace(0.0, 1.0, n_atoms).reshape(-1, 1)
                    kick = np.tile(fibopos, (len(CG), 1))
                
                pos += kick

                self.list2data(name    = name,
                               chain   = chain,
                               resname = resname2,
                               resid   = resid,
                               pos     = pos)


    def construct_water(self):
        bA = self.u.atoms.resname == self.water_resname
        CG = self.u.atoms[bA]

        fibopos   = fibo(r=self.water_fibor, N=self.water_number, verbose=False, plot=None)
        fibopos  -= np.average(fibopos, axis=0)
        waterbead = np.repeat(fibopos, 3, axis=0) + (np.random.rand(self.water_number * 3, 3) - 0.5) * 0.5
        allwater  = np.tile(waterbead, (len(CG), 1)) + \
                    np.repeat(CG[['x','y','z']].to_numpy(), self.water_number * 3, axis=0)
        
        if self.water_chain:
            water_chain = self.water_chain
        else:
            water_chain = 'ZYXWVUTSRQPONMLKJIHGFEDCBA'
        
        if self.water_chain_dms:
            chains = np.repeat([c for c in water_chain[0:self.water_number]] * len(CG), 3)

        else:    
            water_chain = water_chain[0 : (len(water_chain) // self.water_number) * self.water_number ]
            if len(water_chain) % self.water_number != 0:
                assert 0 == 1, 'len(water_chain) % water_number should be 0'
                    
            chains = []
            for i in range(len(water_chain) // self.water_number):
                chains += [c for c in water_chain[i*self.water_number:(i+1)*self.water_number]] * 9999
            chains  = np.repeat(chains[:len(CG) * self.water_number], 3)

        name    = ['OH2', 'H1', 'H2'] * len(CG) * self.water_number
        resname = ['TIP3'] * 3 * len(CG) * self.water_number
        resid   = np.repeat(CG['resid'].to_list(), self.water_number * 3)

        self.list2data(name    = name,
                       chain   = chains,
                       resname = resname,
                       resid   = resid,
                       pos     = allwater)


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


import sqlite3
import numpy as np
from   .universe           import Universe
from   ..utils.sqlcmd      import *
from   ..utils.protein_sel import *
import time
class MartinizeDMS:
    '''
    Take a structural information (DMS),
    and write a new DMS file.

    1) Update atomic information
       - charges, masses, types, nbtypes
    2) Add "bond" = "exclusion"
    3) Add "bond_term"
       - strech_harm, angle_harm, dihedral_trig, improper_trig
    4) Add "nonbonded_info"
    5) Set "nonbonded" terms to 0.0
    6) Add "nonbonded_combined
    '''

    def __init__(self, dms_in, out, martini, epsilon_r=15.0, 
        fcx=10.0, fcy=10.0, fcz=10.0, bfactor_posre=0.5, helix=False,
        addnbtype='ZCARBON'):

        self.martini       = martini
        self.epsilon_r     = epsilon_r
        self.fcx           = fcx
        self.fcy           = fcy
        self.fcz           = fcz
        self.bfactor_posre = bfactor_posre
        self.helix         = helix
        self.addnbtype     = addnbtype

        if not out:
            sp  = dms_in.split('.')
            out = '.'.join(sp[:-1]) + '.martini.dms'


        ### Read a structrual information
        ### Update atomic information
        ### Write a new DMS file
        t1 = time.time()
        universe = Universe(dms_in)
        self.updateAtoms(universe, out)
        t2 = time.time()
        print(f"MartinizeDMS Updating Atoms: {t2-t1:.3f} s")

        ### Now work with a new DMS file
        self.conn     = sqlite3.connect(out)
        self.cursor   = self.conn.cursor()
        self.u        = Universe(out)
        self.resnames = self.u.atoms.resname.unique()

        ### Bonds + Constraints
        t1 = time.time()
        self.bond_param = -1
        self.exclusions = []
        self.updateBonds()
        t2 = time.time()
        print(f"MartinizeDMS Updating Bonds: {t2-t1:.3f} s")

        ### Exclusions -> add manually defined exclusions
        t1 = time.time()
        self.updateExclusions()
        t2 = time.time()
        print(f"MartinizeDMS Updating Exclusions: {t2-t1:.3f} s")


        ### Angles
        t1 = time.time()
        self.angle_param = -1
        self.updateAngles()
        t2 = time.time()
        print(f"MartinizeDMS Updating Angles: {t2-t1:.3f} s")

        ### Dihedrals
        t1 = time.time()
        self.updateDihedrals()
        t2 = time.time()
        print(f"MartinizeDMS Updating Dihedrals: {t2-t1:.3f} s")

        ### virtual_sites3
        self.updateVirtualSites3()

        ### LJ
        t1 = time.time()
        self.updateLJ()
        t2 = time.time()
        print(f"MartinizeDMS Updating LJ: {t2-t1:.3f} s")

        ### Position restraints
        self.updatePosre()

        ### Protein
        t1 = time.time()
        self.updateProtein2()
        t2 = time.time()
        print(f"MartinizeDMS Updating Protein: {t2-t1:.3f} s")

        ### Save DMS
        self.conn.commit()
        self.conn.close()


    def updateAtoms(self, universe, out):
        nbtypes = []
        martinikeys = self.martini.martini['molecules'].keys()

        for index, atom in universe.atoms.iterrows():
            name    = atom['name']
            resname = atom.resname
            
            if resname in martinikeys:
                i = self.martini.martini['molecules'][resname]['atoms']['name'].index(name)
                t = self.martini.martini['molecules'][resname]['atoms']['type'][i]
                q = self.martini.martini['molecules'][resname]['atoms']['q'][i]
                #m = self.martini.martini['atomtypes'][t]['m']
                m = self.martini.martini['molecules'][resname]['atoms']['m'][i]

            else:
                t = self.addnbtype
                q = 0.0
                m = 72.0

            nbtypes.append(t)
            universe.atoms.loc[index, 'charge'] = q / (self.epsilon_r ** 0.5)
            universe.atoms.loc[index, 'mass']   = m
            universe.atoms.loc[index, 'type']   = t

        # Martini beads need to have anum; Otherwise, openMM will not run a simulation
        universe.atoms['anum'] = 6

        nbtypes = sorted(set(nbtypes))
        universe.t2n = {}
        universe.n2t = {}

        for i, nbtype in enumerate(nbtypes):
            universe.t2n[nbtype] = i
            universe.n2t[i]      = nbtype

        for index, atom in universe.atoms.iterrows():
            t = atom['type']
            universe.atoms.loc[index, 'nbtype'] = universe.t2n[t]
            universe.atoms.loc[index, 'type']   = t

        universe.write(out)


    def updateBonds(self):
        for resname in self.resnames:
            if resname not in self.martini.martini['molecules'].keys(): continue
            bonds = self.martini.martini['molecules'][resname]['bonds']
            for bond in bonds:
                atom1 = bond[0]
                atom2 = bond[1]
                r0    = bond[2] * 10
                k     = bond[3] * 0.5 * 2.39e-3

                bA    = self.u.atoms.resname == resname
                bA1   = self.u.atoms.name    == atom1
                bA2   = self.u.atoms.name    == atom2

                df1   = self.u.atoms[bA & bA1]
                df2   = self.u.atoms[bA & bA2]

                assert (df1.resid.to_numpy() == df2.resid.to_numpy()).all(), 'resid'
                assert (df1.chain.to_numpy() == df2.chain.to_numpy()).all(), 'chain'

                for i1, i2 in zip(df1.index, df2.index):
                    self.bond_param += 1

                    if len(bond) == 4:
                        # normal bond
                        self.cursor.execute(sql_insert_exclusion.format(i1, i2))
                        self.exclusions.append([i1, i2])
                        self.exclusions.append([i2, i1])
                        constrained = 0

                    elif len(bond) == 5 and bond[4] == 1:
                        # type 1 constraint
                        self.cursor.execute(sql_insert_exclusion.format(i1, i2))
                        self.exclusions.append([i1, i2])
                        self.exclusions.append([i2, i1])
                        constrained = 1

                    elif len(bond) == 5 and bond[4] == 2:
                        # type 2 constraint (no exclusion)
                        constrained = 1

                    else:
                        assert 0 == 1, 'Something went wrong with updating bonds'

                    btype = resname + '_' + atom1 + ' ' + resname + '_' + atom2
                    self.cursor.execute(sql_insert_bond.format(i1, i2, 1))
                    self.cursor.execute(sql_insert_stretch_harm_term.format(i1, i2, constrained, self.bond_param))
                    self.cursor.execute(sql_insert_stretch_harm_param.format(btype, r0, k, self.bond_param))


    def updateExclusions(self):
        bA  = self.u.atoms.type == self.addnbtype
        df1 = self.u.atoms[bA]
        df2 = self.u.atoms[bA]
        for i1 in df1.index:
            for i2 in df2.index:
                if i1 >= i2: continue
                self.cursor.execute(sql_insert_exclusion.format(i1, i2))

        for resname in self.resnames:
            if resname not in self.martini.martini['molecules'].keys(): continue
            excls = self.martini.martini['molecules'][resname]['exclusions']
            for excl in excls:
                atom1 = excl[0]
                atom2 = excl[1]

                bA    = self.u.atoms.resname == resname
                bA1   = self.u.atoms.name    == atom1
                bA2   = self.u.atoms.name    == atom2

                df1   = self.u.atoms[bA & bA1]
                df2   = self.u.atoms[bA & bA2]

                assert (df1.resid.to_numpy() == df2.resid.to_numpy()).all(), 'resid'
                assert (df1.chain.to_numpy() == df2.chain.to_numpy()).all(), 'chain'

                for i1, i2 in zip(df1.index, df2.index):
                    if [i1, i2] not in self.exclusions:
                        self.cursor.execute(sql_insert_exclusion.format(i1, i2))



    def updateAngles(self):
        for resname in self.resnames:
            if resname not in self.martini.martini['molecules'].keys(): continue
            angles = self.martini.martini['molecules'][resname]['angles']
            for angle in angles:
                atom1 = angle[0]
                atom2 = angle[1]
                atom3 = angle[2]
                t0    = angle[3]
                k     = angle[4] * 0.5 * 0.239
                func  = angle[5]

                bA    = self.u.atoms.resname == resname
                bA1   = self.u.atoms.name    == atom1
                bA2   = self.u.atoms.name    == atom2
                bA3   = self.u.atoms.name    == atom3

                df1   = self.u.atoms[bA & bA1]
                df2   = self.u.atoms[bA & bA2]
                df3   = self.u.atoms[bA & bA3]

                resid1 = df1.resid.to_numpy()
                resid2 = df2.resid.to_numpy()
                resid3 = df3.resid.to_numpy()

                chain1 = df1.chain.to_numpy()
                chain2 = df2.chain.to_numpy()
                chain3 = df3.chain.to_numpy()

                assert (resid1==resid2).all() and (resid2==resid3).all(), 'resid'
                assert (chain1==chain2).all() and (chain2==chain3).all(), 'chain'

                for i1, i2, i3 in zip(df1.index, df2.index, df3.index):
                    self.angle_param += 1
                    atype = resname + '_' + atom1 + ' ' + resname + '_' + atom2 + ' ' + resname + '_' + atom3

                    if func == 2:
                        # cosine function
                        self.cursor.execute(sql_insert_angle_harmcos_term.format(i1, i2, i3, self.angle_param))
                        self.cursor.execute(sql_insert_angle_harmcos_param.format(atype, np.cos(t0 * np.pi / 180), k, self.angle_param))

                    if func == 1:
                        # harmonic function
                        self.cursor.execute(sql_insert_angle_harm_term.format(i1, i2, i3, 0, self.angle_param))
                        self.cursor.execute(sql_insert_angle_harm_param.format(atype, t0, k, self.angle_param))

    def updateVirtualSites3(self):
        lc3_index  = 0
        out3_index = 0
        for resname in self.resnames:
            if resname not in self.martini.martini['molecules'].keys(): continue
            vsites = self.martini.martini['molecules'][resname]['virtual_sites3']
            for vsite in vsites:
                atom1 = vsite[0]
                atom2 = vsite[1]
                atom3 = vsite[2]
                atom4 = vsite[3]
                func  = vsite[4]
                a     = vsite[5]
                b     = vsite[6]

                bA    = self.u.atoms.resname == resname
                bA1   = self.u.atoms.name    == atom1
                bA2   = self.u.atoms.name    == atom2
                bA3   = self.u.atoms.name    == atom3
                bA4   = self.u.atoms.name    == atom4

                df1   = self.u.atoms[bA & bA1]
                df2   = self.u.atoms[bA & bA2]
                df3   = self.u.atoms[bA & bA3]
                df4   = self.u.atoms[bA & bA4]

                resid1 = df1.resid.to_numpy()
                resid2 = df2.resid.to_numpy()
                resid3 = df3.resid.to_numpy()
                resid4 = df4.resid.to_numpy()

                chain1 = df1.chain.to_numpy()
                chain2 = df2.chain.to_numpy()
                chain3 = df3.chain.to_numpy()
                chain4 = df4.chain.to_numpy()

                assert (resid1==resid2).all() and (resid2==resid3).all() and (resid3==resid4).all(), 'resid'
                assert (chain1==chain2).all() and (chain2==chain3).all() and (chain3==chain4).all(), 'chain'

                for i1, i2, i3, i4 in zip(df1.index, df2.index, df3.index, df4.index):
                    if func == 1:
                        # virtual_lc3
                        self.cursor.execute(sql_insert_vlc3_term.format(i1, i2, i3, i4, lc3_index))
                        self.cursor.execute(sql_insert_vlc3_param.format(a, b, lc3_index))
                        lc3_index += 1

                    elif func == 4:
                        # virtual_out3
                        c = vsite[7] * 0.1 #nm-1 to A-1
                        self.cursor.execute(sql_insert_vout3_term.format(i1, i2, i3, i4, out3_index))
                        self.cursor.execute(sql_insert_vout3_param.format(a, b, c, out3_index))
                        out3_index += 1
                


    def updateDihedrals(self):
        self.proper_param   = -1
        self.improper_param = -1
        for resname in self.resnames:
            if resname not in self.martini.martini['molecules'].keys(): continue
            dihedrals = self.martini.martini['molecules'][resname]['dihedrals']
            for dihedral in dihedrals:
                atom1 = dihedral[0]
                atom2 = dihedral[1]
                atom3 = dihedral[2]
                atom4 = dihedral[3]
                func  = dihedral[4]
                t0    = dihedral[5]
                k     = dihedral[6] * 0.239  # 0.5 for improper 1.0 for proper
                n     = dihedral[7]
                fc    = {'0': k,   '1': 0.0, '2': 0.0,
                         '3': 0.0, '4': 0.0, '5': 0.0, '6': 0.0}
                fc[str(n)] = k


                bA    = self.u.atoms.resname == resname
                bA1   = self.u.atoms.name    == atom1
                bA2   = self.u.atoms.name    == atom2
                bA3   = self.u.atoms.name    == atom3
                bA4   = self.u.atoms.name    == atom4

                df1   = self.u.atoms[bA & bA1]
                df2   = self.u.atoms[bA & bA2]
                df3   = self.u.atoms[bA & bA3]
                df4   = self.u.atoms[bA & bA4]

                resid1 = df1.resid.to_numpy()
                resid2 = df2.resid.to_numpy()
                resid3 = df3.resid.to_numpy()
                resid4 = df4.resid.to_numpy()

                chain1 = df1.chain.to_numpy()
                chain2 = df2.chain.to_numpy()
                chain3 = df3.chain.to_numpy()
                chain4 = df4.chain.to_numpy()

                assert (resid1==resid2).all() and (resid2==resid3).all() and (resid3==resid4).all(), 'resid'
                assert (chain1==chain2).all() and (chain2==chain3).all() and (chain3==chain4).all(), 'chain'

                for i1, i2, i3, i4 in zip(df1.index, df2.index, df3.index, df4.index):
                    atype = resname + '_' + atom1 + ' ' + \
                            resname + '_' + atom2 + ' ' + \
                            resname + '_' + atom3 + ' ' + \
                            resname + '_' + atom4 + ' '

                    if func == 1:
                        # proper dihedral
                        self.proper_param += 1
                        self.cursor.execute(sql_insert_dihedral_trig_term.format(i1, i2, i3, i4, self.proper_param))
                        self.cursor.execute(sql_insert_dihedral_trig_param.format(
                            atype, t0, fc['0'], fc['1'], fc['2'], fc['3'], fc['4'], fc['5'], fc['6'], self.proper_param))

                    if func == 2:
                        # improper function
                        self.improper_param += 1
                        # t0 = t0 - np.floor(t0 / 180) * 180
                        self.cursor.execute(sql_insert_improper_harm_term.format(i1, i2, i3, i4, self.improper_param))
                        self.cursor.execute(sql_insert_improper_harm_param.format(atype, t0, k * 0.5, self.improper_param))


    def updateLJ(self):
        t2n = {}
        n2t = {}

        ### set LJ parmeters to 0.0
        for index, atom in self.u.atoms.iterrows():

            if atom.nbtype not in n2t.keys():
                self.cursor.execute(sql_insert_nonbonded.format(
                    atom.nbtype, 0.0, 0.0, atom.type))

            t2n[atom.type]   = atom.nbtype
            n2t[atom.nbtype] = atom.type


        ### set NBFIX
        for ntype1 in sorted(n2t.keys()):
            for ntype2 in sorted(n2t.keys()):
                if ntype1 > ntype2: continue

                ttype1 = n2t[ntype1]
                ttype2 = n2t[ntype2]

                if ((ttype1, ttype2)) in self.martini.martini['nb'].keys():
                    C6, C12 = self.martini.martini['nb'][(ttype1, ttype2)]

                elif ((ttype2, ttype1)) in self.martini.martini['nb'].keys():
                    C6, C12 = self.martini.martini['nb'][(ttype2, ttype1)]

                else:
                    raise AssertionError('{:s} {:s} not exists'.format(ttype1, ttype2))
                
                if C6 == 0.0 and C12 == 0.0:
                    sigma   = 0.0
                    epsilon = 0.0

                else:
                    if self.martini.martini['energy'].startswith('C12'):
                        sigma6  = C12 / C6
                        sigma   = sigma6 ** (1/6)
                        epsilon = C6 / 4 / sigma6

                    elif self.martini.martini['energy'].startswith('4*eps'):
                        sigma   = C6
                        epsilon = C12

                    else:
                        raise AssertionError('LJ function is not defined')

                    sigma   = sigma   * 10      # nm to A
                    epsilon = epsilon * 0.239   # kJ/mol to kcal/mol

                self.cursor.execute(sql_insert_nonbonded_combined.format(
                    ntype1, ntype2, epsilon, sigma, ttype1 + ' ' + ttype2))


    def updatePosre(self):
        ### Find the preexisting posre_index
        try:
            posre_index = pd.read_sql_query("SELECT param FROM posre_harm_term", 
                self.conn)['param'].to_numpy().max() + 1
        except:
            posre_index = 0

        ### Add restraints on the atoms whose bfactors > bfactor_posre
        bA = self.u.atoms.bfactor > self.bfactor_posre
        df = self.u.atoms[bA]

        num_posre = 0
        for index, row in df.iterrows():
            self.cursor.execute(sql_insert_posre_harm_term.format(
                index, row.x, row.y, row.z, posre_index))
            num_posre += 1

        if num_posre != 0:
            self.cursor.execute(sql_insert_posre_harm_param.format(
                self.fcx, self.fcy, self.fcz, posre_index))



    def updateProtein(self):
        prot_resnames = list(three2one.keys())

        ### Add bonds/angles between the resideus
        chains = set(self.u.atoms.chain)

        for chain in chains:
            bA1 = self.u.atoms.chain == chain
            bA2 = self.u.atoms.resname.isin(prot_resnames)
            bA3 = self.u.atoms.name  == 'BB'
            BB  = self.u.atoms[bA1 & bA2 & bA3]

            ### BB bond
            r0 = 0.35 * 10
            k  = 1250 * 0.5 * 2.39e-3

            for i in range(len(BB) - 1):
                CG1 = BB.iloc[i+0]
                CG2 = BB.iloc[i+1]
                
                if CG2.resid - CG1.resid != 1:
                    continue

                i1 = CG1.id
                i2 = CG2.id
                self.bond_param += 1

                self.cursor.execute(sql_insert_exclusion.format(i1, i2))
                btype = f'BB_{i1}_{i2}'

                self.cursor.execute(sql_insert_bond.format(i1, i2, 1))
                self.cursor.execute(sql_insert_stretch_harm_term.format(i1, i2, 0, self.bond_param))
                self.cursor.execute(sql_insert_stretch_harm_param.format(btype, r0, k, self.bond_param))

            
            if self.helix:
                #### BBBB dihedral
                t0 = -120
                k  = 400 * 1.0 * 0.239

                if len(BB.index) > 3.5:
                    for i1, i2, i3, i4 in zip(BB[:-3].index, BB[1:-2].index, BB[2:-1].index, BB[3:].index):
                        self.proper_param += 1
                        dtype = f'BBBB_{i1}_{i2}_{i3}_{i4}'
                        self.cursor.execute(sql_insert_dihedral_trig_term.format(i1, i2, i3, i4, self.proper_param))
                        self.cursor.execute(sql_insert_dihedral_trig_param.format(
                            dtype, t0, k, k, 0.0, 0.0, 0.0, 0.0, 0.0, self.proper_param))


                ### BBB angle params
                t0 = 96
                k  = 700 * 0.5 * 0.239
            
            else:
                ### BBB angle params
                t0 = 127
                k  = 20 * 0.5 * 0.239
            
            
            ### BBB angle
            for i in range(len(BB) - 2):
                CG1 = BB.iloc[i+0]
                CG2 = BB.iloc[i+1]
                CG3 = BB.iloc[i+2]

                if (CG3.resid - CG2.resid != 1) or (CG2.resid - CG1.resid != 1):
                    continue
                
                i1 = CG1.id
                i2 = CG2.id
                i3 = CG3.id
                self.angle_param += 1

                atype = f'BBB_{i1}_{i2}_{i3}'
                self.cursor.execute(sql_insert_angle_harmcos_term.format(i1, i2, i3, self.angle_param))
                self.cursor.execute(sql_insert_angle_harmcos_param.format(atype, np.cos(t0 * np.pi / 180), k, self.angle_param))
            

            ### BBS angle
            t0 = 100
            k  = 25 * 0.5 * 0.239

            bA4 = self.u.atoms.name == 'SC1'
            SC1 = self.u.atoms[bA4]


            for row_index, (index, SC) in enumerate(SC1.iterrows()):

                # # add SBB for the first one
                # if row_index == 0 and SC.resid == BB.iloc[0].resid:
                #     resid = SC.resid
                #     bA5   = BB.resid == resid
                #     bA6   = BB.resid == resid + 1
                #     i1    = index
                #     i2    = BB[bA5].index[0]
                #     i3    = BB[bA6].index[0]
                #     self.angle_param += 1
                #     atype = f'SBB_{i1}_{i2}_{i3}'
                #     self.cursor.execute(sql_insert_angle_harmcos_term.format(i1, i2, i3, self.angle_param))
                #     self.cursor.execute(sql_insert_angle_harmcos_param.format(atype, np.cos(t0 * np.pi / 180), k, self.angle_param))

                # BBS
                resid = SC.resid
                bA5   = BB.resid == resid - 1
                bA6   = BB.resid == resid

                atom1 = BB[bA5]
                atom2 = BB[bA6]

                if len(atom1) * len(atom2) == 1:
                    i1 = atom1.index[0]
                    i2 = atom2.index[0]
                    i3 = index

                    self.angle_param += 1
                    atype = f'BBS_{i1}_{i2}_{i3}'
                    self.cursor.execute(sql_insert_angle_harmcos_term.format(i1, i2, i3, self.angle_param))
                    self.cursor.execute(sql_insert_angle_harmcos_param.format(atype, np.cos(t0 * np.pi / 180), k, self.angle_param))


    def updateProtein2(self):
        prot_resnames = list(three2one.keys())

        ### Add bonds/angles between the residues
        chains = set(self.u.atoms.chain)

        bA2 = self.u.atoms.resname.isin(prot_resnames)
        bA3 = self.u.atoms.name == 'BB'
        BB  = self.u.atoms[bA2 & bA3].copy()
        bA4 = self.u.atoms.name == 'SC1'
        SC1 = self.u.atoms[bA2 & bA4].copy()

        ### Precompute necessary values
        BB_ids    = BB.index.to_numpy()
        BB_resids = BB.resid.to_numpy()
        BB_chains = BB.chain.to_numpy()

        SC_ids    = SC1.index.to_numpy()
        SC_resids = SC1.resid.to_numpy()
        SC_chains = SC1.chain.to_numpy()


        ### BB bond parameters
        r0 = 0.35 * 10
        k  = 1250 * 0.5 * 2.39e-3

        t0 = 100
        kt = 25 * 0.5 * 0.239

        for i in range(len(BB) - 1):
            if BB_chains[i+1] != BB_chains[i]: continue
            if BB_resids[i+1] - BB_resids[i] != 1: continue

            i1 = BB_ids[i]
            i2 = BB_ids[i+1]
            self.bond_param += 1

            self.cursor.execute(sql_insert_exclusion.format(i1, i2))
            btype = f'BB_{i1}_{i2}'

            self.cursor.execute(sql_insert_bond.format(i1, i2, 1))
            self.cursor.execute(sql_insert_stretch_harm_term.format(i1, i2, 0, self.bond_param))
            self.cursor.execute(sql_insert_stretch_harm_param.format(btype, r0, k, self.bond_param))

            SC_id = SC_ids[(BB_chains[i+1] == SC_chains) & (BB_resids[i+1] == SC_resids)]
            if len(SC_id) != 1: continue
            i3 = SC_id[0]

            self.angle_param += 1
            atype = f'BBS_{i1}_{i2}_{i3}'
            self.cursor.execute(sql_insert_angle_harmcos_term.format(i1, i2, i3, self.angle_param))
            self.cursor.execute(sql_insert_angle_harmcos_param.format(atype, np.cos(t0 * np.pi / 180), kt, self.angle_param))


        ### BBB angle parameters
        t0, k = (96, 700 * 0.5 * 0.239) if self.helix else (127, 20 * 0.5 * 0.239)
        for i in range(len(BB) - 2):
            if BB_chains[i+2] != BB_chains[i+1]: continue
            if BB_chains[i+1] != BB_chains[i+0]: continue
            if (BB_resids[i+2] - BB_resids[i+1] != 1): continue
            if (BB_resids[i+1] - BB_resids[i+0] != 1): continue

            i1 = BB_ids[i+0]
            i2 = BB_ids[i+1]
            i3 = BB_ids[i+2]
            self.angle_param += 1

            atype = f'BBB_{i1}_{i2}_{i3}'
            self.cursor.execute(sql_insert_angle_harmcos_term.format(i1, i2, i3, self.angle_param))
            self.cursor.execute(sql_insert_angle_harmcos_param.format(atype, np.cos(t0 * np.pi / 180), k, self.angle_param))
        

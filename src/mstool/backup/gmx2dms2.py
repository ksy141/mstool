import sqlite3
import numpy as np
import os
from   .universe           import Universe
from   ..utils.sqlcmd      import *

class GMX2DMS2:
    '''
    I thought this would be faster than GMX2DMS...
    But it is surprisingly slightly slower.

    Take structural information + gromacs 
    and write a fully parameterized DMS file.

    1) Update atomic information
       - charges, masses, types, nbtypes
    2) Add "bond" = "exclusion"
    3) Add "bond_term"
       - strech_harm, angle_harm, dihedral_trig, improper_trig
    4) Add "nonbonded_info"
    5) Set "nonbonded" terms to 0.0
    6) Add "nonbonded_combined
    '''

    def __init__(self, structure, martini, out, epsilon_r=15.0, use_structure_resid=True):
        self.universe  = Universe(structure)
        self.martini   = martini
        self.epsilon_r = epsilon_r
        self.use_structure_resid = use_structure_resid

        ### Make a DMS file
        if os.path.exists(out): os.remove(out)
        self.conn     = sqlite3.connect(out)
        self.cursor   = self.conn.cursor()

        ### CREATE TABLES
        self.cursor.executescript(sql_create)

        ### Add Cell
        self.updateCell()

        ### Collect nbtypes
        self.collectNB()

        ### LJ
        self.updateLJ()

        ### Add Atoms
        self.update()

        ### Save DMS
        self.conn.commit()
        self.conn.close()


    def updateCell(self):
        self.cursor.execute(sql_insert_cell.format(1,*self.universe.cell[0]))
        self.cursor.execute(sql_insert_cell.format(2,*self.universe.cell[1]))
        self.cursor.execute(sql_insert_cell.format(3,*self.universe.cell[2]))

    def collectNB(self):
        '''self.t2n = {'Qa': 0, 'Na': 1}
           self.n2t = {0: 'Qa', 1: 'Na'}'''

        nbtypes = []
        for mol in self.martini.martini['mols']:
            nbtypes.extend(self.martini.martini['molecules'][mol[0]]['atoms']['type'])

        nbtypes = sorted(set(nbtypes))
        self.t2n = {}
        self.n2t = {}

        for i, nbtype in enumerate(nbtypes):
            self.t2n[nbtype] = i
            self.n2t[i] = nbtype


    def updateLJ(self):

        ### check nbfunc
        nbfunc = self.martini.martini['nbfunc']
        assert nbfunc == 'LJ', '{nbfunc:s} is not supported'

        if 'C6' in self.martini.martini['energy']:
            comb = 1
        elif 'sigma' in self.martini.martini['energy']:
            comb = 2
        else:
            assert 0 == 1, 'check the nonbonded function'


        ### set all sigma and epsilon to zero
        for nbtype, attype in self.n2t.items():
            self.cursor.execute(sql_insert_nonbonded.format(
                nbtype, 0.0, 0.0, attype))

        ### set NBFIX
        for ntype1 in sorted(self.n2t.keys()):
            for ntype2 in sorted(self.n2t.keys()):
                if ntype1 > ntype2: continue

                ttype1 = self.n2t[ntype1]
                ttype2 = self.n2t[ntype2]

                if ((ttype1, ttype2)) in self.martini.martini['nb'].keys():
                    Vii, Wii = self.martini.martini['nb'][(ttype1, ttype2)]

                elif ((ttype2, ttype1)) in self.martini.martini['nb'].keys():
                    Vii, Wii = self.martini.martini['nb'][(ttype2, ttype1)]

                else:
                    raise AssertionError('{:s} {:s} not exists'.format(ttype1, ttype2))
                
                if comb == 1:
                    C6      = Vii
                    C12     = Wii
                    sigma6  = C12 / C6
                    sigma   = sigma6 ** (1/6)
                    epsilon = C6 / 4 / sigma6

                elif comb == 2:
                    sigma   = Vii
                    epsilon = Wii

                sigma   = sigma   * 10     # nm to A
                epsilon = epsilon * 0.239  # kJ/mol to kcal/mol

                self.cursor.execute(sql_insert_nonbonded_combined.format(
                    ntype1, ntype2, epsilon, sigma, ttype1 + ' ' + ttype2))



    def update(self):
        insertion = ''; msys_ct = 0.0
        self.bond_param = 0
        self.angle_param = 0
        proper_param = 0
        improper_param = 0
        posre_index = 0
        index = -1

        ### loop over molecule types
        for mol in self.martini.martini['mols']:
            molecule = self.martini.martini['molecules'][mol[0]]
            nmols    = mol[1]

            ### loop over molecule copies
            for i in range(nmols):
                natoms = len(molecule['atoms']['name'])


                ### atoms
                '''Combining information from gromacs and structure.
                Gromacs: name, resname, mass, charge, nbtype, type, resid
                Structure: x, y, z, anum, chain, bfactor

                resid is tricky. 
                If there is only one molecule, follow gromacs, as it will be likely a protein chain.
                If there are more than one molecule, recalculate resid.
                '''
                for j in range(natoms):
                    i1     = index + j + 1
                    aname  = molecule['atoms']['name'][j]
                    atype  = molecule['atoms']['type'][j]
                    charge = molecule['atoms']['q'][j] / self.epsilon_r ** 0.5
                    nbtype = self.t2n[atype]
                    mass   = self.martini.martini['atomtypes'][atype]['m']
                    resn   = molecule['atoms']['resname'][j]

                    anum, chain, resid_s, bfactor, x, y, z = \
                    self.universe.atoms[['anum', 'chain', 'resid', 'bfactor', 
                                         'x',  'y',  'z']].loc[i1]

                    # Martini beads need to have anum; Otherwise, openMM will not run a simulation
                    anum = 6

                    # assigning resid
                    if nmols == 1:
                        resid = molecule['atoms']['resid'][j]
                    elif nmols != 1 and self.use_structure_resid:
                        resid = resid_s
                    elif nmols != 1 and not self.use_structure_resid:
                        resid = i + 1

                    # save
                    self.cursor.execute(sql_insert_particle, (
                        i1, anum, aname, resn, chain, resid, mass, charge,
                        x, y, z, 0.0, 0.0, 0.0, chain, 
                        insertion, msys_ct, nbtype, atype, bfactor))


                ### Bonds
                bonds = molecule['idx_bonds']
                for bond in bonds:

                    i1 = index + bond[0]
                    i2 = index + bond[1]
                    r0 = bond[2] * 10.0
                    k  = bond[3] * 0.5 * 2.39e-3

                    btype = mol[0] + '_' + \
                            molecule['atoms']['name'][bond[0] - 1] + '_' + \
                            molecule['atoms']['name'][bond[1] - 1] 

                    if len(bond) == 4:
                        # normal bond
                        self.cursor.execute(sql_insert_exclusion.format(i1, i2))

                    elif len(bond) == 5 and bond[4] == 1:
                        # type 1 constraint
                        self.cursor.execute(sql_insert_exclusion.format(i1, i2))

                    else:
                        # type 2 constraint (no exclusion)
                        pass

                    self.cursor.execute(sql_insert_bond.format(i1, i2, 1))
                    self.cursor.execute(sql_insert_stretch_harm_term.format(i1, i2, 0, self.bond_param))
                    self.cursor.execute(sql_insert_stretch_harm_param.format(btype, r0, k, self.bond_param))
                    self.bond_param += 1



                ### Exclusions
                exclusions = molecule['idx_exclusions']

                for excl in exclusions:
                    i1 = index + excl[0]
                    i2 = index + excl[1]
                    self.cursor.execute(sql_insert_exclusion.format(i1, i2))



                ### Angles
                angles = molecule['idx_angles']

                for angle in angles:
                    i1   = index + angle[0]
                    i2   = index + angle[1]
                    i3   = index + angle[2]
                    t0   = angle[3]
                    k    = angle[4] * 0.5 * 0.239
                    func = angle[5]

                    atype = mol[0] + '_' + \
                            molecule['atoms']['name'][angle[0] - 1] + '_' + \
                            molecule['atoms']['name'][angle[1] - 1] + '_' + \
                            molecule['atoms']['name'][angle[2] - 1]

                    if func == 2:
                        # cosine function
                        self.cursor.execute(sql_insert_angle_harmcos_term.format(i1, i2, i3, self.angle_param))
                        self.cursor.execute(sql_insert_angle_harmcos_param.format(atype, np.cos(t0 * np.pi / 180), k, self.angle_param))

                    if func == 1:
                        # harmonic function
                        self.cursor.execute(sql_insert_angle_harm_term.format(i1, i2, i3, 0, self.angle_param))
                        self.cursor.execute(sql_insert_angle_harm_param.format(atype, t0, k, self.angle_param))

                    self.angle_param += 1



                ### Dihedrals
                dihes  = molecule['idx_dihedrals']

                for dihe in dihes:
                    i1 = index + dihe[0]
                    i2 = index + dihe[1]
                    i3 = index + dihe[2]
                    i4 = index + dihe[3]

                    func  = dihe[4]
                    t0    = dihe[5]
                    k     = dihe[6] * 0.239  # 0.5 for improper 1.0 for proper
                    n     = dihe[7]
                    fc    = {'0': k,   '1': 0.0, '2': 0.0,
                             '3': 0.0, '4': 0.0, '5': 0.0, '6': 0.0}
                    fc[str(n)] = k

                    dtype = mol[0] + '_' + \
                            molecule['atoms']['name'][dihe[0] - 1] + '_' + \
                            molecule['atoms']['name'][dihe[1] - 1] + '_' + \
                            molecule['atoms']['name'][dihe[2] - 1] + '_' + \
                            molecule['atoms']['name'][dihe[3] - 1]

                    if func == 1:
                        # proper dihedral
                        self.cursor.execute(sql_insert_dihedral_trig_term.format(i1, i2, i3, i4, proper_param))
                        self.cursor.execute(sql_insert_dihedral_trig_param.format(
                            dtype, t0, fc['0'], fc['1'], fc['2'], fc['3'], fc['4'], fc['5'], fc['6'], proper_param))
                        proper_param += 1

                    if func == 2:
                        # improper dihedral
                        self.cursor.execute(sql_insert_improper_harm_term.format(i1, i2, i3, i4, improper_param))
                        self.cursor.execute(sql_insert_improper_harm_param.format(dtype, t0, k * 0.5, improper_param))
                        improper_param += 1



                ### Positional Restraints
                posres = molecule['idx_position_restraints']

                for posre in posres:
                    assert int(posre[1]) == 1, 'Only harmonic positional restraints are supported'

                    i1 = index + posre[0]

                    ### CHECK THE FUNCTION (UNCHECKED)
                    fcx = posre[2] * 0.5 * 2.39e-3
                    fcy = posre[3] * 0.5 * 2.39e-3
                    fcz = posre[4] * 0.5 * 2.39e-3

                    x, y, z = self.universe.atoms[['x',  'y',  'z']].loc[index]

                    self.cursor.execute(sql_insert_posre_harm_term.format(
                        i1, x, y, z, posre_index))

                    self.cursor.execute(sql_insert_posre_harm_param.format(
                        fcx, fcy, fcz, posre_index))

                    posre_index += 1


                ### Go to a next molecule copy
                index += natoms



    def updateVirtualSites3(self):
        virtual_index = 0
        index = -1

        # loop over molecule types
        for mol in self.martini.martini['mols']:
            molecule = self.martini.martini['molecules'][mol[0]]
            nmols = mol[1]

            # loop over molecule copies
            for i in range(nmols):
                natoms = len(molecule['atoms']['name'])
                vsites = molecule['idx_virtual_sites3']

                for vsite in vsites:
                    i1 = index + vsite[0]
                    i2 = index + vsite[1]
                    i3 = index + vsite[2]
                    i4 = index + vsite[3]
                    func = vsite[4]
                    a  = vsite[5]
                    b  = vsite[6]

                index += natoms

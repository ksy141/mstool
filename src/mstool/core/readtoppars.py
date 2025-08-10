import numpy  as np
import pandas as pd
import glob
import os
import networkx as nx
pwd  = os.path.dirname(os.path.realpath(__file__))

elements = ['LP', 'H',  'He', 'Li', 'Be', 'B',  'C',  'N',  'O',  'F',  'Ne',
            'Na', 'Mg', 'Al', 'Si', 'P',  'S',  'Cl', 'Ar', 'K',  'Ca', 'Sc',
            'Ti', 'V',  'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge',
            'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y',  'Zr', 'Nb', 'Mo', 'Tc',
            'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I',  'Xe',
            'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb',
            'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W',  'Re', 'Os',
            'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr',
            'Ra', 'Ac', 'Th', 'Pa', 'U',  'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf',
            'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt',
            'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']

masses = np.array([  0.        ,   1.008     ,   4.002602  ,   6.94      ,
                     9.0121831 ,  10.81      ,  12.011     ,  14.007     ,
                    15.999     ,  18.99840316,  20.1797    ,  22.98976928,
                    24.305     ,  26.9815385 ,  28.085     ,  30.973762  ,
                    32.06      ,  35.45      ,  39.948     ,  39.0983    ,
                    40.078     ,  44.955908  ,  47.867     ,  50.9415    ,
                    51.9961    ,  54.938044  ,  55.845     ,  58.933194  ,
                    58.6934    ,  63.546     ,  65.38      ,  69.723     ,
                    72.63      ,  74.921595  ,  78.971     ,  79.904     ,
                    83.798     ,  85.4678    ,  87.62      ,  88.90584   ,
                    91.224     ,  92.90637   ,  95.95      ,  97.90721   ,
                   101.07      , 102.9055    , 106.42      , 107.8682    ,
                   112.414     , 114.818     , 118.71      , 121.76      ,
                   127.6       , 126.90447   , 131.293     , 132.90545196,
                   137.327     , 138.90547   , 140.116     , 140.90766   ,
                   144.242     , 144.91276   , 150.36      , 151.964     ,
                   157.25      , 158.92535   , 162.5       , 164.93033   ,
                   167.259     , 168.93422   , 173.045     , 174.9668    ,
                   178.49      , 180.94788   , 183.84      , 186.207     ,
                   190.23      , 192.217     , 195.084     , 196.966569  ,
                   200.592     , 204.38      , 207.2       , 208.9804    ,
                   209.        , 210.        , 222.        , 223.        ,
                   226.        , 227.        , 232.0377    , 231.03588   ,
                   238.02891   , 237.        , 244.        , 243.        ,
                   247.        , 247.        , 251.        , 252.        ,
                   257.        , 258.        , 259.        , 262.        ,
                   267.        , 268.        , 271.        , 274.        ,
                   269.        , 276.        , 281.        , 281.        ,
                   285.        , 286.        , 289.        , 288.        ,
                   293.        , 294.        , 294.        ])

comment = '!'
kcal2kJ = 4.184
items   = ('RESI', 'ATOMS', 'BONDS', 'ANGLES', 'DIHEDRALS', 'IMPROPER', 'NONBONDED', 'NBFIX', 'PRES', 'CMAP')
items1  = (        'ATOMS', 'BONDS', 'ANGLES', 'DIHEDRALS', 'IMPROPER', 'NONBONDED', 'NBFIX',               )


def DFS(G):
    nodes   = [node for node in G.nodes()]
    visited = set()
    queue   = [nodes[0]]

    while queue:
        current = queue.pop()
        if current not in visited:
            visited.add(current)
            for neighbor in G.neighbors(current):
                if neighbor not in visited:
                    queue.append(neighbor)

    return nodes, visited


class ReadToppars:

    def __init__(self, toppars=pwd + '/../FF/charmm36.toppar/toppar.str', create_group=True, verbose=False):

        # a list of files to read
        if isinstance(toppars, list):
            pass

        # indicate a path
        elif os.path.isdir(toppars):
            toppars = glob.glob(toppars + '/*rtf') + \
                      glob.glob(toppars + '/*prm') + \
                      glob.glob(toppars + '/*str')

        elif os.path.isfile(toppars):
            p = os.path.dirname(toppars)
            if p == '': p = '.'

            with open(toppars) as f:
                tops = []

                for line in f.readlines():
                    if line.startswith(('open', 'stream')):
                        tops.append(p + '/' + line.split()[-1])

            toppars = tops


        if not toppars:
            raise ValueError('please provide a list of top/prm/str or a path to them')

        self.create_group = create_group
        self.GROUP = []
        self.verbose = verbose
        print(pd.DataFrame([toppar.split('/')[-1] for toppar in toppars]))
        self.toppars = toppars
        self.topparmols = {}
        self.parse()
        self.collect()
        if self.create_group: self.remove_duplicate_graphs()


    def remove_duplicate_graphs(self):
        hashes = set()
        final = []
        for subG in self.GROUP:
            if len(subG.nodes) < 3.5:
                continue
            h = nx.weisfeiler_lehman_graph_hash(subG)
            if h not in hashes:
                hashes.add(h)
                final.append(subG)
        self.GRAPH = sorted(final, key=lambda g: len(g.nodes))


    def parse(self):
        '''Remove comments and parse useful lines'''

        save  = ''
        saves = []
        atoms = []

        read_type = False
        for toppar in self.toppars:
            toppar_name = toppar.split('/')[-1]
            self.topparmols[toppar_name] = []

            for line in open(toppar, 'r'):

                rline = line.rstrip()
                if rline: # Pass empty line

                    sline = rline.split(comment)[0].rstrip()
                    if sline: # Pass a line starting with !

                        if sline.startswith('HBOND'):
                            continue

                        if sline.startswith('MASS'):
                            atoms.append(sline.split())
                            continue

                        if sline.startswith(items):
                            if sline.startswith(('RESI', 'PRES')):
                                self.topparmols[toppar_name].append(sline.split()[1])
                                bA_CMAP = False

                            # Hoping CMAP is not the first entity in any files.
                            if sline.startswith(items1):
                                bA_CMAP = True

                            if sline.startswith('CMAP') and not bA_CMAP:
                                save += sline + '\n'
                                continue

                            else:
                                read_type = sline.split()[0]
                                saves.append(save)
                                save = ''


                        if sline.startswith(('END', 'end')):
                            read_type = False
                            saves.append(save)
                            save = ''

                        # READ EACH COMPONENT
                        if read_type:
                            if sline.endswith('-'):
                                save += sline
                            else:
                                save += sline + '\n'

        self.atoms = atoms
        self.saves = saves



    def collect(self):
        self.RESI      = {}
        self.ATOMS     = {}
        self.BONDS     = {}
        self.ANGLES    = {}
        self.DIHEDRALS = {}
        self.IMPROPER  = {}
        self.NONBONDED = {}
        self.NB14      = {}
        self.NBFIX     = {}
        self.CMAP      = []

        for atom in self.atoms:
            #if len(atom) == 5: continue # don't read an atom line with element

            atomtype = atom[2]
            atommass = float(atom[3])

            diffmass = np.abs(masses - atommass)
            numb = int(np.argmin(diffmass))
            elem = elements[numb]
            if diffmass[numb] > 1.5:
                print('confirm: atomtype (%s) and atomelem (%s)' %(atomtype, elem))

            self.ATOMS[atomtype] = {}
            self.ATOMS[atomtype]['mass'] = atommass
            self.ATOMS[atomtype]['elem'] = elem
            self.ATOMS[atomtype]['numb'] = numb

        for save in self.saves:
            if not save: continue
            ssave = save.rstrip().split('\n')

            if save.startswith(('RESI', 'PRES')):
                self._read_RESI(ssave)

            #elif save.startswith('ATOMS'):
            #    self._read_ATOMS(ssave)

            elif save.startswith('BONDS'):
                self._read_BONDS(ssave)

            elif save.startswith('ANGLES'):
                self._read_ANGLES(ssave)

            elif save.startswith('DIHEDRALS'):
                self._read_DIHEDRALS(ssave)

            elif save.startswith('IMPROPER'):
                self._read_IMPROPER(ssave)

            elif save.startswith('NONBONDED'):
                self._read_NONBONDED(ssave)

            elif save.startswith('NBFIX'):
                self._read_NBFIX(ssave)

            elif save.startswith('CMAP'):
                self._read_CMAP(ssave)



    def _read_RESI(self, ssave):
        resname = ssave[0].split()[1].upper()
        rtype   = ssave[0].split()[0].upper()

        names     = []
        types     = []
        charges   = []
        masses    = []
        bonds     = []
        imprs     = []
        angles    = []
        dihedrals = []
        cmaps     = []
        groups    = []
        g         = []

        for line in ssave:
            segments = line.split()

            if line.startswith('GROUP'):
                if g:
                    groups.append(g)
                    g = []

            elif line.startswith('ATOM'):
                type = segments[2].upper()
                names.append(segments[1])
                types.append(type)
                charges.append(float(segments[3]))
                g.append([segments[1], type, float(segments[3])])

                if type not in self.ATOMS.keys():
                    masses.append(0)
                else:
                    masses.append(self.ATOMS[type]['mass'])

            elif line.startswith(('BOND', 'DOUB')):
                for i in range(0, len(segments[1:]), 2):
                    bonds.append([segments[i+1], segments[i+2]])

            elif line.startswith('ANGL'):
                for i in range(0, len(segments[1:]), 3):
                    angles.append([segments[i+1], segments[i+2], segments[i+3]])

            elif line.startswith('DIHE'):
                for i in range(0, len(segments[1:]), 4):
                    dihedrals.append([segments[i+1], segments[i+2], segments[i+3], segments[i+4]])

            elif line.startswith('IMPR'):
                for i in range(0, len(segments[1:]), 4):
                    imprs.append([segments[i+1], segments[i+2], segments[i+3], segments[i+4]])

            elif line.startswith('CMAP'):
                cmaps.append([segments[1], segments[2], segments[3], segments[4], segments[8]])

        if resname in self.RESI.keys() and self.verbose:
            out = 'Duplicated residue name: ' + resname + '\n'
            for key, value in self.topparmols.items():
                if resname in value:
                    out += '{:<50}: {:<10}'.format(key, resname)
            print(out + '\n')

        if g:
            groups.append(g)

        if self.create_group and rtype == 'RESI':
            for group in groups:
                if len(group) < 3.5:
                    continue

                try:
                    G = nx.Graph()
                    G.resname = resname
                    for atom in group:
                        G.add_node(atom[0], type=atom[1], charge=atom[2],
                                   atomic_number=self.ATOMS[atom[1]]['numb'])
                    for atomi, atomj in bonds:
                        if atomi in G and atomj in G:
                            G.add_edge(atomi, atomj)

                    if len(G.edges) == 0:
                        continue

                    DFS_nodes, DFS_visited = DFS(G)
                    if len(DFS_nodes) == len(DFS_visited):
                        self.GROUP.append(G)

                except Exception as e:
                    print(e)
                    print(resname, group)
                    pass

        self.RESI[resname] = {'names':     np.array(names),
                              'types':     np.array(types),
                              'masses':    np.array(masses),
                              'bonds':     np.array(bonds),
                              'imprs':     np.array(imprs),
                              'cmaps':     np.array(cmaps),
                              'charges':   np.array(charges),
                              'angles':    np.array(angles),
                              'dihedrals': np.array(dihedrals)}



    def _read_ATOMS(self, ssave):
        for line in ssave:
            if not line.startswith('MASS'): continue

            segments = line.split()
            if segments[2] not in self.ATOMS.keys():
                self.ATOMS[segments[2]] = float(segments[3])


    def _read_BONDS(self, ssave):
        for line in ssave[1:]:
            segments = line.split()
            type1 = segments[0]
            type2 = segments[1]

            if type1 > type2:
                type1, type2 = type2, type1

            # kcal/mol/A**2 -> kJ/mol/nm**2 incl factor 2
            Kb    = float(segments[2]) * 2 * kcal2kJ * 1000 / 10
            b0    = float(segments[3]) / 10

            self.BONDS[(type1, type2)] = [b0, Kb]



    def _read_ANGLES(self, ssave):
        for line in ssave[1:]:
            segments = line.split()

            type1 = segments[0]
            type2 = segments[1]
            type3 = segments[2]

            if type1 > type3:
                type1, type3 = type3, type1

            cth   = float(segments[3]) * 2 * kcal2kJ
            th0   = float(segments[4])

            Kub   = 0.0
            S0    = 0.0

            if len(segments) > 6: # check for Urey-Bradley parameters
                Kub = float(segments[5]) * 2 * kcal2kJ * 1000 / 10
                S0  = float(segments[6]) / 10


            self.ANGLES[(type1, type2, type3)] = [th0, cth, S0, Kub]




    def _read_DIHEDRALS(self, ssave):
        for line in ssave[1:]:
            segments = line.split()

            type1 = segments[0]
            type2 = segments[1]
            type3 = segments[2]
            type4 = segments[3]
            cp   = float(segments[4]) * kcal2kJ
            mult = int(segments[5])
            phi0 = float(segments[6])


            #if type1 > type4:
            #    type1, type2, type3, type4 = type4, type3, type2, type1

            #if type1 == type4 and type2 > type3:
            #    type2, type3 = type3, type2

            i4   = (type1, type2, type3, type4)
            j4   = (type4, type3, type2, type1)

            if i4 in self.DIHEDRALS:
                for imult, dihedral in enumerate(self.DIHEDRALS[i4]):
                    if dihedral[2] == mult: del self.DIHEDRALS[i4][imult]
                self.DIHEDRALS[i4].append([phi0, cp, mult])
            elif j4 in self.DIHEDRALS:
                for imult, dihedral in enumerate(self.DIHEDRALS[j4]):
                    if dihedral[2] == mult: del self.DIHEDRALS[j4][imult]
                self.DIHEDRALS[j4].append([phi0, cp, mult])
            else:
                self.DIHEDRALS[i4] = [[phi0, cp, mult]]

            #if i4 in self.DIHEDRALS:
            #    for imult, dihedral in enumerate(self.DIHEDRALS[i4]):
            #        if dihedral[2] == mult: del self.DIHEDRALS[i4][imult]
            #    self.DIHEDRALS[i4].append([phi0, cp, mult])

            #else:
            #    self.DIHEDRALS[i4] = [[phi0, cp, mult]]



    def _read_IMPROPER(self, ssave):
        for line in ssave[1:]:
            segments = line.split()

            type1 = segments[0]
            type2 = segments[1]
            type3 = segments[2]
            type4 = segments[3]

            if type1 > type4:
                type1, type2, type3, type4 = type4, type3, type2, type1

            cq = float(segments[4]) * 2 * kcal2kJ
            q0 = float(segments[6])

            self.IMPROPER[(type1, type2, type3, type4)] = [q0, cq]



    def _read_NONBONDED(self, ssave):
        ### LENNARD JONES
        ### [CHARMM]   V_LJ = 1 * eps * ( (r_0/ r)^12  - 2 * (r_0 / r)^6)
        ### [GROMACS]  V_LJ = 4 * eps * ( (sigma/r)^12 - (sigma/r)^6 )
        ### r_0 = RminHalf * 2 = 2^(1/6) * sigma

        for line in ssave[1:]:
            segments = line.split()

            attype   = segments[0]
            epsilon  = float(segments[2])
            RminHalf = float(segments[3])
            eps      = abs(epsilon*kcal2kJ)                # conversion to kJ and positive
            sigma    = 2*RminHalf/(10.0*2.0**(1.0/6.0))    # -> nm, double distance and rmin2sigma factor

            self.NONBONDED[attype] = [sigma, eps]

            if len(segments)> 6:        # test length to avoid IndexError
                try:                    # if possible, convert element 5 to float
                    segments[5] = float(segments[5])
                except:
                    None

                # is segment 5 and 6 floats => there's 1-4 defined
                if not isinstance(segments[5], str):                  # not string?
                    epsilon14  = float(segments[5])                   # read charmm epsilon
                    eps14      = abs(epsilon14*kcal2kJ)               # conversion to gromacs units
                    Rmin14Half = float(segments[6])                   # read charmm Rmin*1/2
                    sigma14    = 2*Rmin14Half/(10.0*2.0**(1.0/6.0))   # conversion to gromacs units
                    self.NB14[attype] = [sigma14, eps14]              # add to list



    def _read_NBFIX(self, ssave):
        for line in ssave[1:]:
            if not line: continue
            if line.startswith('HBOND'): continue

            segments = line.split()
            type1 = segments[0]
            type2 = segments[1]

            if type1 > type2:
                type1, type2 = type2, type1

            epsilon = float(segments[2])
            Rmin    = float(segments[3])
            eps     = abs(epsilon*kcal2kJ)          # conversion to kJ and positive
            sigma   = Rmin/(10.0*2.0**(1.0/6.0))    # -> nm, double distance and rmin2sigma factor

            self.NBFIX[(type1, type2)] = [sigma, eps]


    def _read_CMAP(self, ssave):
        for line in ssave[1:]:
            segments = line.split()
            try:
                segments[0] = float(segments[0])
                for icmap in segments:
                    icmap = float(icmap)*kcal2kJ
                    self.CMAP[-1]['data'].append(icmap)

            except:
                type1 = segments[0]
                type2 = segments[1]
                type3 = segments[2]
                type4 = segments[3]
                type5 = segments[7]
                ncmap = int(segments[8])
                self.CMAP.append(dict(atoms=[type1, type2, type3, type4, type5], ncmap=ncmap, data=[]))

        for i in range(len(self.CMAP)):
            self.CMAP[i]['data'] = np.array(self.CMAP[i]['data'])


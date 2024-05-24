import pandas as pd
import numpy  as np
import sqlite3
from   ..utils.mdautils import *
from   ..utils.atomic_data import anum_mass_vdw_from_name
from   ..utils.amberselection import *
from   ..utils.sqlcmd  import *
from   .readxml        import ReadXML
from   ..lib.distance  import distance_matrix
import os
pd.set_option('display.max_rows', 500)

# caution: name is a predefined variable for pandas
# use atom['name'] instead of atom.name

# to iterate the data frame
# for index, atom in universe.atoms.iterrows():
#     universe.atoms.loc[index, 'charge'] = 2
# to change the values
# u.atoms.loc[u.atoms.resid == 1, 'resid'] = 2

# bfactor
# if dms file doesn't have bfactor, make a column with a default value of 0.0
# if dms file has bfactor, read values from there
# if pdb file is provided, reset bfactor value to 0.0


def custom_sort(value):
    # 0 ... 9 A ... Z '1' ... '9'
    if isinstance(value, int):
        return (0, value)
    elif value.isdigit():
        return (2, value)
    else:
        return (1, value)


class Universe:
    def __init__(self, data=None, sort=False, fixResid=False):
        '''Create an universe object. 
        Provide a structure file (.pdb or .dms),
        dictionary as data argument, or
        pandas DataFrame as data argument.
        If nothing is provided, an empty universe will be created.

        Parameters
        ----------
        data : str or dict or pd.DataFrame
            Filename of a structure file (*.dms or *.pdb). 
            If a structure file is not provided, use this to create a universe object.
            The default is None.
        sort : bool
            Sort before adding resn. The default is False.
        
        Attributes
        ----------
        atoms : pd.DataFrame
        cell : np.array
            a numpy array of ``(3, 3)``.
        dimensions : np.array
            a numpy array of ``(6,)``.
        bonds : list
            a list of bonds in the universe. 
            Initialized to [] unless DMS file contains bonds or 
            you create bonds
        '''
        
        ### INITIAL VALUES
        self.cell       = [[1,0,0], [0,1,0], [0,0,1]] #triclinic_vectors
        self.dimensions = [0, 0, 0, 90, 90, 90]       #triclinic_box
        self.bonds = []
        self.cols  = {'anum': -1, 'name': 'tbd', 'charge': 0.0, 'mass': 0.0,
                      'type': 'tbd', 'nbtype': 0, 'resname': 'tbd', 'resid': 0, 
                      'segname': '', 'chain': 'S', 
                      'x': 0.0, 'y': 0.0, 'z': 0.0, 'bfactor': 0.0, 'vdw': 0.0}


        ### READ DATA
        self.ext = None
        if isinstance(data, str):
            if not os.path.exists(data):
                raise ValueError(data, " does not exist")
            ext = data.split('.')[-1]
            if ext.lower() not in ['pdb', 'dms', 'gro']:
                raise ValueError('only pdb, dms, and gro are supported')
            if ext == 'pdb' or ext == 'PDB':
                self.readPDB(data)
            if ext == 'dms' or ext == 'DMS':
                self.readDMS(data)
            if ext == 'gro' or ext == 'GRO':
                self.readGRO(data)
            self.ext = ext.lower()

        elif isinstance(data, dict):
            self.construct_from_dict(data)

        elif isinstance(data, pd.DataFrame):
            self.construct_from_df(data)

        else:
            self.atoms = pd.DataFrame(columns = self.cols)
        self.atoms.universe = self
        

        ### IF atomic number / mass / vdw is not in input
        anum_sum = np.sum(self.atoms.anum == -1)
        mass_sum = np.sum(self.atoms.mass == 0.0)
        vdw_sum  = np.sum(self.atoms.vdw  == 0.0)
        anum_mass_vdw_sum  = anum_sum + mass_sum + vdw_sum

        if anum_mass_vdw_sum > 0.5:
            anum, mass, vdw = self.update_anum_mass_vdw_from_name()
            if anum_sum > 0.5: self.atoms.anum = anum
            if mass_sum > 0.5: self.atoms.mass = mass
            if vdw_sum  > 0.5: self.atoms.vdw  = vdw


        # if sort is needed before adding resn
        if sort: self.sort()

        # add "resn"
        self.addResidues()
        
        # fix resid
        if fixResid: self.fixResid()


    @property
    def n_atoms(self):
        return len(self.atoms)

    @property
    def N(self):
        return len(self.atoms)

    @property
    def positions(self):
        return self.atoms[['x','y','z']].to_numpy()

    def setpositions(self, pos):
        self.atoms.loc[:,['x','y','z']] = pos

    def construct_from_dict(self, data):
        for col in self.cols:
            if col not in data.keys():
                data[col] = self.cols[col]
        self.atoms = pd.DataFrame(data)

    def construct_from_df(self, data):
        self.atoms = data.copy().reset_index()
        for col in self.cols:
            if col not in self.atoms.columns:
                self.atoms[col] = self.cols[col]

    def select(self, string, returnbA=False):
        if '~' in string or '>' in string or '<' in string:
            print('~ > <  currently not supported')
            return np.array(len(self.atoms) * [False])

        # first check out how many amber selection blocks are
        amberblocks = strsplitbysepexc(string)
        #print(amberblocks)

        # amberselect for each amberblock
        bAdict = {}
        for i in range(len(amberblocks)):
            if amberblocks[i] == 'all':
                bAdict[f'bA{i}'] = np.array(len(self.atoms) * [True])

            else:
                dictselect = amberSelection(amberblocks[i])
                #print(dictselect)
                bAchain   = np.array(len(self.atoms) * [True])
                bAresid   = np.array(len(self.atoms) * [True])
                bAresname = np.array(len(self.atoms) * [True])
                bAname    = np.array(len(self.atoms) * [True])

                for key, value in dictselect.items():
                    if key == 'resid':
                        bAresid = self.atoms.resid.isin(value)
                    if key == 'chain':
                        value = [str(val) for val in value]
                        valstar, valnostar = separateStarNoStar(value)
                        bAstar   = self.atoms.chain.str.startswith(tuple(valstar))
                        bAnostar = self.atoms.chain.isin(valnostar)
                        bAchain  = bAstar | bAnostar
                    if key == 'resname':
                        valstar, valnostar = separateStarNoStar(value)
                        bAstar    = self.atoms.resname.str.startswith(tuple(valstar))
                        bAnostar  = self.atoms.resname.isin(valnostar)
                        bAresname = bAstar | bAnostar
                    if key == 'name':
                        valstar, valnostar = separateStarNoStar(value)
                        bAstar   = self.atoms['name'].str.startswith(tuple(valstar))
                        bAnostar = self.atoms['name'].isin(valnostar)
                        bAname   = bAstar | bAnostar
                
                bAdict[f'bA{i}'] = bAchain & bAresid & bAresname & bAname
        
        # now replace amberblock with bA0, bA1, ...
        bAstring = ''
        for s in strsplitbysepinc(string):
            if s in amberblocks:
                bAstring += 'bA' + str(amberblocks.index(s))
            else:
                bAstring += s
        
        if returnbA:
            return eval(bAstring, bAdict)
        else:
            return self.atoms[eval(bAstring, bAdict)]

    def changeName(self, change):
        for key, value in change.items():
            bA = self.select(key, returnbA=True)

            for key2, value2 in amberSelection(value).items():
                self.atoms.loc[bA, key2] = value2[0]

    def clone(self, string):
        u = Universe(self.select(string))
        u.dimensions = self.dimensions
        u.cell = self.cell
        return u

    def write(self, ofile, wrap=False):
        if wrap: self.wrapMolecules()

        ext = ofile.split('.')[-1]
        if ext == 'pdb' or ext == 'PDB':
            self.writePDB(ofile)
        elif ext == 'dms' or ext == 'DMS':
            self.writeDMS(ofile)
        elif ext == 'gro' or ext == 'GRO':
            self.writeGRO(ofile)
        else:
            print('the file format should be either pdb or dms')

    def wrapMolecules(self):
        xavg, yavg, zavg = self.atoms[['x','y','z']].mean()
        xeval = np.round(xavg / (self.dimensions[0]/2))
        yeval = np.round(yavg / (self.dimensions[1]/2))
        zeval = np.round(zavg / (self.dimensions[2]/2))

        if xeval == 0 and yeval == 0 and zeval == 0:
            shiftfunc = lambda vec: np.round(vec)
        elif xeval == 1 and yeval == 1 and zeval == 1:
            shiftfunc = lambda vec: np.floor(vec)
        else:
            print("wrapping molecules could be strange")
            shiftfunc = lambda vec: np.round(vec)
        

        resnxyz = self.atoms[['x','y','z','resn']].groupby(['resn']).mean()
        dim     = self.dimensions[0:3]
        for index, value in resnxyz.iterrows():
            pos     = np.array([value['x'], value['y'], value['z']])
            unit    = shiftfunc(pos / dim)
            
            # translate
            if not np.all(unit == 0):
                self.atoms.loc[self.atoms.resn == index, ['x','y','z']] -= dim * unit

    def readPDB(self, ifile):
        '''will read only the lines that start with ATOM or CRYST1'''
        data = {'name': [], 'resname': [], 'resid': [], 'chain': [], 
                'x': [], 'y': [], 'z': [], 'bfactor': [], 'segname': []}
        
        with open(ifile) as W:
            for line in W.readlines():
                if line.startswith('CRYST1'):
                    self.dimensions = np.array([line[6:15], line[15:24],
                                                line[24:33], line[33:40],
                                                line[40:47], line[47:54]],
                                                dtype=np.float32)
                    self.cell = triclinic_vectors(self.dimensions)

                if line.startswith(('ATOM', 'HETATM')):
                    name = line[12:16].strip()
                    data['name'].append(name)
                    data['resname'].append(line[17:21].strip())
                    data['resid'].append(int(line[22:26].strip()))
                    data['chain'].append(line[21].strip())
                    data['x'].append(float(line[30:38]))
                    data['y'].append(float(line[38:46]))
                    data['z'].append(float(line[46:54]))
                    data['bfactor'].append(float(line[60:66]))
                    try:
                        data['segname'].append(line[66:76].strip())
                    except:
                        data['segname'].append('')
        
        self.construct_from_dict(data)

    def writePDB(self, ofile, spacegroup='P 1', zvalue=1):
        #fmt is copied from MDAnalysis

        fmt = {
            'ATOM': (
                "ATOM  {:5d} {:<4s}{:<1s}{:<4s}"
                "{:1s}{:4d}{:1s}"
                "   {:8.3f}{:8.3f}{:8.3f}{:6.2f}"
                "{:6.2f}      {:<4s}{:>2s}\n"),
            'CRYST1': ("CRYST1{:9.3f}{:9.3f}{:9.3f}"
                "{:7.2f}{:7.2f}{:7.2f} "
                "{:<11s}{:4d}\n"),
        }
 
        with open(ofile, 'w') as W:
            W.write(fmt['CRYST1'].format(
                *self.dimensions,
                spacegroup, zvalue))
            
            i = 1
            for index, atom in self.atoms.iterrows():
                name = atom['name']

                # if len(name) = 1,2,3, add space
                if len(name) in [1,2,3]: name = f' {name:s}'
                W.write(fmt['ATOM'].format(
                    ltruncate_int(i, 5), name[:4], '', atom.resname[:4],
                    atom.chain[:1] if atom.chain else '', ltruncate_int(atom.resid, 4), '',
                    atom.x, atom.y, atom.z, 0.0, atom.bfactor,
                    atom.segname, ''))
                i += 1
    
    def readGRO(self, ifile):
        data = {'name': [], 'resname': [], 'resid': [], 
                'x': [], 'y': [], 'z': []}

        with open(ifile) as fp:
            lines = fp.readlines()
            for i in range(len(lines)):
                if i == 0:
                    continue
                elif i == 1:
                    n_atoms = int(lines[i].strip())
                else:
                    sl = lines[i].strip().split()
                    if len(sl) < 3:
                        continue

                    if isfloat(sl[0]) and isfloat(sl[1]) and isfloat(sl[2]):
                        if len(sl) == 3:
                            dimx = float(sl[0].strip()) * 10.0
                            dimy = float(sl[1].strip()) * 10.0
                            dimz = float(sl[2].strip()) * 10.0
                            self.dimensions = [dimx, dimy, dimz, 90, 90, 90]
                            self.cell       = triclinic_vectors(self.dimensions)

                        if len(sl) == 9:
                            v1x = float(sl[0].strip()) * 10.0
                            v2y = float(sl[1].strip()) * 10.0
                            v3z = float(sl[2].strip()) * 10.0
                            v1y = float(sl[3].strip()) * 10.0
                            v1z = float(sl[4].strip()) * 10.0
                            v2x = float(sl[5].strip()) * 10.0
                            v2z = float(sl[6].strip()) * 10.0
                            v3x = float(sl[7].strip()) * 10.0
                            v3y = float(sl[8].strip()) * 10.0
                            x = [v1x, v1y, v1z]
                            y = [v2x, v2y, v2z]
                            z = [v3x, v3y, v3z]
                            self.dimensions = triclinic_box(x, y, z)
                            self.cell       = triclinic_vectors(self.dimensions)

                    else:
                        data['resid'].append(int(lines[i][0:5].strip()))
                        data['resname'].append(lines[i][5:10].strip())
                        data['name'].append(lines[i][10:15].strip())
                        data['x'].append(float(lines[i][20:28].strip()) * 10.0)
                        data['y'].append(float(lines[i][28:36].strip()) * 10.0)
                        data['z'].append(float(lines[i][36:44].strip()) * 10.0)
        
        if n_atoms != len(data['x']):
            print(f"""WARNING: n_atoms is not consistent: {n_atoms} and {len(data['x'])}""")
        self.construct_from_dict(data)

    def writeGRO(self, ofile):
        with open(ofile, 'w') as W:
            # n_atoms (the note field should not be empty; otherwise VMD does not like it)
            W.write("   \n{0:5d}\n".format(len(self.atoms)))
            
            # atoms
            i = 1
            for index, atom in self.atoms.iterrows():
                idx     = ltruncate_int(i, 5)
                name    = atom['name']
                resid   = ltruncate_int(atom['resid'], 5)
                resname = atom['resname']
                x       = atom['x'] * 0.1
                y       = atom['y'] * 0.1
                z       = atom['z'] * 0.1
                W.write(f"{resid:>5d}{resname:<5.5s}{name:>5.5s}{idx:>5d}{x:8.3f}{y:8.3f}{z:8.3f}\n")
                i += 1

            # cell
            b = np.array(self.cell).flatten() * 0.1
            if b[1] == 0 and b[2] == 0 and b[5] == 0:
                W.write(f"{b[0]:10.5f} {b[4]:9.5f} {b[8]:9.5f}")
            else:
                W.write(f"{b[0]:10.5f} {b[4]:9.5f} {b[8]:9.5f} {b[1]:9.5f} {b[2]:9.5f} {b[3]:9.5f} {b[5]:9.5f} {b[6]:9.5f} {b[7]:9.5f}")
 

    def readDMS(self, ifile):
        conn = sqlite3.connect(ifile, isolation_level=None, detect_types=sqlite3.PARSE_COLNAMES)
        self.atoms = pd.read_sql_query("SELECT * FROM PARTICLE", conn)
        for col in self.cols:
            if col not in self.atoms.columns:
                self.atoms[col] = self.cols[col]

        try:
            self.bonds = conn.execute('SELECT p0,p1 FROM bond').fetchall()
        except:
            pass

        try:        
            for i, (x, y, z) in enumerate(conn.execute('SELECT x, y, z FROM global_cell')):
                self.cell[i][0] = x
                self.cell[i][1] = y
                self.cell[i][2] = z

            self.dimensions = triclinic_box(self.cell[0], self.cell[1], self.cell[2])
        except:
            pass

        conn.close()


    def writeDMS(self, ofile):
        if os.path.exists(ofile): os.remove(ofile)
        conn   = sqlite3.connect(ofile)
        cursor = conn.cursor()

        ### CREATE TABLES
        cursor.executescript(sql_create)
        # for line in sql_create.split(';'):
        #     cursor.execute(line.strip())

        ### CELLS
        cursor.execute(sql_insert_cell.format(1,*self.cell[0]))
        cursor.execute(sql_insert_cell.format(2,*self.cell[1]))
        cursor.execute(sql_insert_cell.format(3,*self.cell[2]))
        
        ### PARTICLES        
        msys_ct = 0
        vx = 0.0; vy=0.0; vz=0.0;
        insertion = ''
        
        i = 0
        for index, atom in self.atoms.iterrows():
            cursor.execute(sql_insert_particle, (
                i, atom.anum, atom['name'], atom.resname,
                atom.chain, atom.resid, atom.mass, atom.charge,
                atom.x, atom.y, atom.z, vx, vy, vz, atom.segname,
                insertion, msys_ct, atom.nbtype, atom.type, atom.bfactor,
                ))
            i += 1

        ### BONDS
        for bond in self.bonds:
            cursor.execute(sql_insert_bond.format(bond[0], bond[1], 1))
        
        ### Save
        conn.commit()
        conn.close()


    def addBondFromXML(self, ff=[], ff_add=[]):
        self.bonds = []
        xml = ReadXML(ff=ff, ff_add=ff_add)
        resnames = set(self.atoms.resname)

        for resname in set(self.atoms.resname):
            if resname not in xml.RESI.keys():
                print(f'Warning: openMM xml does not have {resname}')
                continue

            bonds = xml.RESI[resname]['bonds']

            for bond in bonds:
                bA    = self.atoms.resname == resname
                bA0   = self.atoms.name == bond[0]
                bA1   = self.atoms.name == bond[1]

                atomA = self.atoms[bA & bA0].index
                atomB = self.atoms[bA & bA1].index

                assert len(atomA) == len(atomB), f":{resname}@{bond[0]},{bond[1]}={len(atomA)},{len(atomB)}"

                for a, b in zip(atomA, atomB):
                    self.bonds.append([a, b])
    
    def addBondFromDistance(self):
        # two atoms are bonded if 0.1 < d < 0.55 * (R1 + R2)
        if self.dimensions[0] * self.dimensions[1] * self.dimensions[2] < 1.5:
            dist = distance_matrix(self.positions, self.positions, dimensions=None)
        else:
            dist = distance_matrix(self.positions, self.positions, dimensions=self.dimensions)

        vdw = self.atoms.vdw.values
        bA1 = dist < 0.55 * (np.array([vdw]).T + np.array([vdw]))
        bA2 = 0.2  < dist
        bA3 = np.tri(*bA1.shape).astype(bool)
        self.bonds = np.transpose(np.where(bA1 & bA2 & bA3))
        
    def addResidues(self):
        '''Add residue number. The atoms must be sorted before this.
        Universe.atoms do not have a residue level unlike MDAnalysis.
        While this has been useful for adding atoms to the Universe,
        selecting atoms based on residues is difficult.
        To make residue-based selection easier,
        this function adds residue numbers.'''

        if len(self.atoms) == 0: return
        
        # remove a residue number
        if 'resn' in self.atoms.keys():
            self.atoms.drop('resn', axis=1, inplace=True)
        
        dg = self.atoms
        
        # fast
        chains = dg.chain.values
        resids = dg.resid.values
        resnas = dg.resname.values

        bAchain = chains[1:] == chains[:-1]
        bAresid = resids[1:] == resids[:-1]
        bAresna = resnas[1:] == resnas[:-1]
        bAfinal = bAchain & bAresna & bAresid
        dg['resn'] = np.insert(np.cumsum((~bAfinal)), 0, 0)
        self.atoms['resn'] = dg['resn']


    def fixResid(self):
        lastres = max(self.atoms.resn)
        grouped = self.atoms.groupby('resn')
        modify  = False
        diff    = 0
        previous_group = None
        updates = {}

        if self.ext == 'pdb':
            limit = 10000
        elif self.ext == 'gro':
            limit = 100000
        else:
            return

        for resn, current_group in grouped:
            if previous_group is not None:
                # new chain
                if previous_group['chain'].iloc[0] != current_group['chain'].iloc[0]:
                    modify = False
                    diff = 0
                else:
                    resid0 = previous_group['resid'].iloc[0]
                    resid1 = current_group['resid'].iloc[0]
    
                    # Check for transition in resid numbers
                    if resid0 % limit == (limit-1) and (resid1 == 0 or resid1 == 1):
                        modify = True
                        diff = limit * (resid0 // limit + 1)
    
                # Apply modification if flagged
                if modify:
                    updates[tuple(current_group.index)] = resid1 + diff
    
            # Move to next group
            previous_group = current_group
    
        # Apply all accumulated changes in one operation
        for idx_tuple, new_resid in updates.items():
            self.atoms.loc[list(idx_tuple), 'resid'] = new_resid


    def sort(self, by=['chain', 'resid'], ignore_index=True):
        #LoopModeler is not happy with it
        #self.atoms = self.atoms.sort_values(by=by)
        #self.atoms['newindex'] = np.arange(0, len(self.atoms))
        #
        #bonds_tmp = []
        #for bond in self.bonds:
        #    i1 = bond[0]
        #    i2 = bond[1]
        #    bonds_tmp.append([self.atoms.loc[i1].newindex, 
        #                      self.atoms.loc[i2].newindex])
        #self.bonds = bonds_tmp

        #self.atoms = self.atoms.sort_values(by=by, ignore_index=ignore_index)
        self.atoms = self.atoms.sort_values(by=by, ignore_index=ignore_index, key=lambda x: x.apply(custom_sort))
        self.bonds = []

    
    def update_anum_mass_vdw_from_name(self):
        names  = self.atoms['name'].str.upper().values
        masses = np.zeros(len(names))
        anums  = np.zeros(len(names), dtype=np.int64) 
        vdws   = np.zeros(len(names))

        for i, name in enumerate(names):
            upper_name = name.upper()
            if upper_name in anum_mass_vdw_from_name.keys():
                anum, mass, vdw = anum_mass_vdw_from_name[upper_name]

            elif upper_name[0] in anum_mass_vdw_from_name.keys():
                anum, mass, vdw = anum_mass_vdw_from_name[upper_name[0]]

            else:
                anum, mass, vdw = 6, 12.0, 1.0

            anums[i]  = anum
            masses[i] = mass
            vdws[i]   = vdw
                
        return anums, masses, vdws

    def box_volume(self, dimensions):
        """Return the volume of the unitcell described by `dimensions`.

        The volume is computed as the product of the box matrix trace, with the
        matrix obtained from :func:`triclinic_vectors`.
        If the box is invalid, i.e., any box length is zero or negative, or any
        angle is outside the open interval ``(0, 180)``, the resulting volume will
        be zero.

        Parameters
        ----------
        dimensions : array_like
            Unitcell dimensions provided in the same format as returned by
            :attr:`MDAnalysis.coordinates.timestep.Timestep.dimensions`:\n
            ``[lx, ly, lz, alpha, beta, gamma]``.

        Returns
        -------
        volume : float
            The volume of the unitcell. Will be zero for invalid boxes.


        .. versionchanged:: 0.20.0
            Calculations are performed in double precision and zero is returned
            for invalid dimensions.
        """
        dim = np.asarray(dimensions, dtype=np.float64)
        lx, ly, lz, alpha, beta, gamma = dim
        if alpha == beta == gamma == 90.0 and lx > 0 and ly > 0 and lz > 0:
            # valid orthogonal box, volume is the product of edge lengths:
            volume = lx * ly * lz
        else:
            # triclinic or invalid box, volume is trace product of box matrix
            # (invalid boxes are set to all zeros by triclinic_vectors):
            tri_vecs = triclinic_vectors(dim, dtype=np.float64)
            volume = tri_vecs[0, 0] * tri_vecs[1, 1] * tri_vecs[2, 2]
        return volume


    def makeTopology(self, martini=None, top=None):
        counts = []; c = 0
        residue_previous = -1
        resname_previous = 'PREVIOUS'
        for index, atom in self.atoms.iterrows():
            resname = atom['resname']
            residue = atom['resn']

            if resname_previous != resname:
                # adding previous residue information
                counts.append([resname_previous, c])

                # new resname is introduced!
                c = 1
                resname_previous = resname
                residue_previous = residue

            elif resname_previous == resname and residue_previous != residue:
                # same resname but different residue
                c += 1
                resname_previous = resname
                residue_previous = residue

            elif resname_previous == resname and residue_previous == residue:
                # same resname and residue
                pass

        # the final residue
        counts.append([resname, c])

        # remove ['PREVIOUS', 0]
        counts.pop(0)

        stdout  = '; \n'
        stdout += '; topol.top made by mstool\n'
        stdout += '; You should modify this file based on the below '
        stdout += 'because it counts only correctly for '
        stdout += 'lipids / waters / ions.\n'
        stdout += '; It does not work for protein chain. \n'
        stdout += '; You need to run martinize.py anyway to parameterize your protein. \n'
        stdout += '; \n\n'
        
        if martini:
            for ifile in martini.ifiles:
                stdout += f'#include "{ifile}"\n'

        stdout += '\n[ system ]\nMartini\n\n[ molecules ]\n'
        for count in counts:
            stdout += '{:10s} {:10d}\n'.format(count[0], count[1])

        print(stdout)
        if top:
            with open(top, 'w') as W:
                W.write(stdout)

        return counts



def Merge(*args, sort=False, ignore_index=True):
    '''Merge atomic groups.

    Parameters
    ----------
    args : pd.DataFrame atom objects
    sort : bool
        The default is False.
    ignore_index : bool
        The default is True.

    Returns
    -------
    u : Universe
        a merged Universe.

    Examples
    --------
    >>> u = mstool.Merge(u1.atoms, u2.atoms, u3.atoms)
    >>> u.dimensions = u1.dimensions
    >>> u.cell       = u1.cell
    >>> u.write('final.pdb')
    '''

    #index = 0
    #for arg in args:
    #    bonds.extend(list(np.array(arg.universe.bonds) + index))
    #    index += len(arg)
    
    u = Universe()

    atoms = [arg for arg in args if len(arg) > 0]
    if len(atoms) == 0:
        return u
    
    u.atoms = pd.concat(atoms, ignore_index=ignore_index)
    if sort: u.sort()
    return u


def RemoveOverlappedResidues(atoms1, atoms2, rcut, dimensions=None, returnoverlapped=False):
    '''Remove molecules that have steric clashes.
    As Universe object consists of atom-based objects rather than residue-based,
    residue information is constructed based on resid, resname, and chain.
    Equivalent line in MDAnalysis will be 
    ``u.select_atoms('(sel2) or (not (same residue as (sel1 and within rcut of sel2)))')``.''

    Parameters
    ----------
    atoms1 : pd.DataFrame or Universe.atoms
        Molecules that are close to atoms2 will be removed. (e.g., solvent)
    atoms2 : pd.DataFrame or Universe.atoms
        Molecules in this group will be kept if they are close to atoms1 (e.g., protein)
    rcut : float
        Cutoff distance in A.
    dimensions : np.array
        If provided, PBC will be taken into account.
        Plus, the returned Universe object will have this dimensions.
        The default is None. (PBC will not be taken into account)

    Returns
    -------
    u : Universe
        Universe that has atoms2 and some of the molecules of atoms1 
        that do not have close contacts with atoms2.
        Please note that in the returned Universe,
        atoms2 come first and then a part of atoms1.
    '''

    from ..lib.distance import distance_overlap
    u1 = Universe(data=atoms1)
    u2 = Universe(data=atoms2)

    if len(u2.atoms) == 0:
        if returnoverlapped:
            return [False] * len(u1.atoms)

        else:
            return u1

    pos1 = u1.atoms[['x','y','z']].to_numpy()
    pos2 = u2.atoms[['x','y','z']].to_numpy()

    if dimensions is not None:
        far = distance_overlap(pos1, pos2, rcut, dimensions)
    else:
        far = distance_overlap(pos1, pos2, rcut)

    # find residue numbers that are close to atoms2
    resn = u1.atoms.loc[~far, 'resn']
    bA   = ~u1.atoms['resn'].isin(resn)

    if returnoverlapped:
        return ~bA

    u = Merge(atoms2, atoms1[bA], ignore_index=True)

    if dimensions is not None:
        u.dimensions = dimensions
        u.cell       = triclinic_vectors(dimensions)

    return u



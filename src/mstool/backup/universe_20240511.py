import pandas as pd
import numpy  as np
import sqlite3
from   ..utils.amberselection import *
from   ..utils.sqlcmd  import *
from   .readxml        import ReadXML
import os
import re
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

anum_mass_from_name = {
   'H':   [1,   1.00],
   'C':   [6,  12.00],
   'N':   [7,  14.00],
   'O':   [8,  16.00],
   'F':   [9,  19.00],
   'NA':  [11, 23.00],
   'SOD': [11, 23.00],  
   'MG':  [12, 24.30],
   'P':   [15, 31.00],
   'S':   [16, 32.00],
   'CL':  [17, 35.00],
   'CLA': [17, 35.00],
   'ZN':  [30, 65.38],
   'W':   [8,  16.00],
}


update_anum_mass_cmdline = '''
UPDATE particle SET 

anum = CASE 
WHEN name = "CL"     THEN 17
WHEN name = "Cl"     THEN 17
WHEN name = "CLA"    THEN 17
WHEN name = "NA"     THEN 11
WHEN name = "Na"     THEN 11
WHEN name = "SOD"    THEN 11
WHEN name = "Zn"     THEN 30
WHEN name = "ZN"     THEN 30
WHEN name = "MG"     THEN 12
WHEN name = "Mg"     THEN 12
WHEN name = "BB"     THEN 6
WHEN name = "W"      THEN 8
WHEN name LIKE "SC%" THEN 6 
WHEN name LIKE "H%"  THEN 1
WHEN name LIKE "C%"  THEN 6
WHEN name LIKE "N%"  THEN 7
WHEN name LIKE "O%"  THEN 8
WHEN name LIKE "F%"  THEN 9
WHEN name LIKE "P%"  THEN 15
WHEN name LIKE "S%"  THEN 16
ELSE 1
END, 

mass = CASE 
WHEN name = "CL"     THEN 35.00
WHEN name = "Cl"     THEN 35.00
WHEN name = "CLA"    THEN 35.00
WHEN name = "NA"     THEN 23.00
WHEN name = "Na"     THEN 23.00
WHEN name = "SOD"    THEN 23.00
WHEN name = "Zn"     THEN 65.38
WHEN name = "ZN"     THEN 65.38
WHEN name = "MG"     THEN 24.30
WHEN name = "Mg"     THEN 24.30
WHEN name = "BB"     THEN 12.00
WHEN name = "W"      THEN 16.00
WHEN name LIKE "SC%" THEN 12.00
WHEN name LIKE "H%"  THEN 1.00
WHEN name LIKE "C%"  THEN 12.00
WHEN name LIKE "N%"  THEN 14.00
WHEN name LIKE "O%"  THEN 16.00
WHEN name LIKE "F%"  THEN 19.00
WHEN name LIKE "P%"  THEN 31.00
WHEN name LIKE "S%"  THEN 32.00
ELSE 0.0
END;
'''


class Universe:
    def __init__(self, ifile=None, data=None, create_bonds=False, ff=[], ff_add=[], 
                 guess_atomic_number=False):
        '''Create an universe object. 
        Provide a structure file (.pdb or .dms),
        dictionary as data argument, or
        pandas DataFrame as data argument.
        If nothing is provided, an empty universe will be created.

        Parameters
        ----------
        ifile : str
            Filename of a structure file (*.dms or *.pdb). 
            The default is None.
        data : dict or pd.DataFrame
            If a structure file is not provided, use this to create a universe object.
            The default is None.
        create_bonds : bool
            create bonds for non-protein molecules using forcefield information.
            The default is None.
        ff : list
        ff_add : list
        guess_atomic_number : bool
            Guess atomic number and mass from atom names. 
            This part is significantly slow. Do not use it unless you have to.
            The default is ``False``.
        wrap : bool
            Wrap moleculees when writing a structure. Default is ``False``.

        
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
            you create bonds via ``create_bonds=True`` 
        '''
        
        #if not isinstance(add_columns, list): add_columns = [add_columns]
        self.cell       = [[1,0,0], [0,1,0], [0,0,1]] #triclinic_vectors
        self.dimensions = [0, 0, 0, 90, 90, 90]       #triclinic_box
        self.bonds = []
        self.cols  = {'anum': -1, 'name': 'tbd', 'charge': 0.0, 'mass': 0.0,
                      'type': 'tbd', 'nbtype': 0, 'resname': 'tbd', 'resid': 0, 
                      'segname': '', 'chain': 'A', 
                      'x': 0.0, 'y': 0.0, 'z': 0.0, 'bfactor': 0.0}
        self.create_bonds = create_bonds

        if ifile:
            if not os.path.exists(ifile):
                raise ValueError(ifile, " does not exist")
            ext = ifile.split('.')[-1]
            if ext not in ['pdb', 'dms']:
                raise ValueError('only pdb and dms are supported')
            if ext == 'pdb' or ext == 'PDB':
                self.readPDB(ifile)
            if ext == 'dms' or ext == 'DMS':
                self.readDMS(ifile)

        elif isinstance(data, dict):
            self.construct_from_dict(data)

        elif isinstance(data, pd.DataFrame):
            self.construct_from_df(data)

        else:
            self.atoms = pd.DataFrame(columns = self.cols)
        
        self.atoms.universe = self

        # Update atomic numbers and masses based on names
        if guess_atomic_number: self.update_anum_mass()

        self.xml      = ReadXML(ff=ff, ff_add=ff_add)
        self.resnames = set(self.atoms.resname)
        self.addResidues()


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

    def select(self, string):
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

        return self.atoms[eval(bAstring, bAdict)]


    def write(self, ofile, guess_atomic_number=True, wrap=False):
        if wrap: self.wrapMolecules()

        ext = ofile.split('.')[-1]
        if ext == 'pdb' or ext == 'PDB':
            self.writePDB(ofile)
        elif ext == 'dms' or ext == 'DMS':
            if self.create_bonds: self.addBonds()
            self.writeDMS(ofile, guess_atomic_number)
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
                    self.cell = self.triclinic_vectors(self.dimensions)

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
                    self.ltruncate_int(i, 5), name[:4], '', atom.resname[:4],
                    atom.chain[:1], self.ltruncate_int(atom.resid, 4), '',
                    atom.x, atom.y, atom.z, 0.0, atom.bfactor,
                    atom.segname, ''))
                i += 1


    def ltruncate_int(self, value, ndigits):
        """Truncate an integer, retaining least significant digits
    
        Parameters
        ----------
        value : int
          value to truncate
        ndigits : int
          number of digits to keep
    
        Returns
        -------
        truncated : int
          only the `ndigits` least significant digits from `value`
    
        Examples
        --------
        >>> ltruncate_int(123, 2)
        23
        >>> ltruncate_int(1234, 5)
        1234
        """
        return int(str(value)[-ndigits:])

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

            self.dimensions = self.triclinic_box(self.cell[0], self.cell[1], self.cell[2])
        except:
            pass

        conn.close()


    def writeDMS(self, ofile, guess_atomic_number):
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
        
        ### Update anum and mass
        if guess_atomic_number: cursor.execute(update_anum_mass_cmdline)
        
        ### Save
        conn.commit()
        conn.close()


    def addBonds(self):
        for resname in self.resnames:
            if resname not in self.xml.RESI.keys():
                print(f'Warning: openMM xml does not have {resname}')
                continue

            bonds = self.xml.RESI[resname]['bonds']

            for bond in bonds:
                bA    = self.atoms.resname == resname
                bA0   = self.atoms.name == bond[0]
                bA1   = self.atoms.name == bond[1]

                atomA = self.atoms[bA & bA0].index
                atomB = self.atoms[bA & bA1].index

                assert len(atomA) == len(atomB), "the length of atoms for bonds is different"

                for a, b in zip(atomA, atomB):
                    self.bonds.append([a, b])



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
        
        # make a tmp instance that is sorted
        # dg = self.atoms.sort_values(by=['chain', 'resid', 'resname'])
        dg = self.atoms
        
        # slow
        #self.atoms['resn'] = -1
        #resna_prev = '-1000'
        #resid_prev = -1000
        #chain_prev = '-1000'
        #value = -1

        #for i in range(len(dg)):
        #    resna = dg.iloc[i].resname
        #    resid = dg.iloc[i].resid
        #    chain = dg.iloc[i].chain
        #    index = dg.iloc[i].name

        #    if not (resid == resid_prev and chain == chain_prev and resna == resna_prev):
        #        value += 1
        #        resid_prev = resid
        #        resna_prev = resna
        #        chain_prev = chain

        #    self.atoms.loc[index, 'resn'] = value

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


    def guess_anum_mass(self, name):
        anum = 0
        mass = 0.0

        try:
            return anum_mass_from_name[name.upper()]
        except:
            pass

        try:
            return anum_mass_from_name[name.upper()[0]]
        except:
            pass

        return anum, mass


    def update_anum_mass(self):
        for index, atom in self.atoms.iterrows():
            name = atom['name']
            anum, mass = self.guess_anum_mass(name)
            self.atoms.loc[index, 'mass'] = mass
            self.atoms.loc[index, 'anum'] = anum


    ### the below is copied from MDAnalysis 2.4.1
    def triclinic_box(self, x, y, z):
        """Convert the three triclinic box vectors to
        ``[lx, ly, lz, alpha, beta, gamma]``.

        If the resulting box is invalid, i.e., any box length is zero or negative,
        or any angle is outside the open interval ``(0, 180)``, a zero vector will
        be returned.

        All angles are in degrees and defined as follows:

        * ``alpha = angle(y,z)``
        * ``beta  = angle(x,z)``
        * ``gamma = angle(x,y)``

        Parameters
        ----------
        x : array_like
            Array of shape ``(3,)`` representing the first box vector
        y : array_like
            Array of shape ``(3,)`` representing the second box vector
        z : array_like
            Array of shape ``(3,)`` representing the third box vector

        Returns
        -------
        numpy.ndarray
            A numpy array of shape ``(6,)`` and dtype ``np.float32`` providing the
            unitcell dimensions in the same format as returned by
            :attr:`MDAnalysis.coordinates.timestep.Timestep.dimensions`:\n
            ``[lx, ly, lz, alpha, beta, gamma]``.\n
            Invalid boxes are returned as a zero vector.

        Note
        ----
        Definition of angles: http://en.wikipedia.org/wiki/Lattice_constant

        See Also
        --------
        :func:`~MDAnalysis.lib.mdamath.triclinic_vectors`


        .. versionchanged:: 0.20.0
           Calculations are performed in double precision and invalid box vectors
           result in an all-zero box.
        """
        x = np.asarray(x, dtype=np.float64)
        y = np.asarray(y, dtype=np.float64)
        z = np.asarray(z, dtype=np.float64)
        lx = np.linalg.norm(x)
        ly = np.linalg.norm(y)
        lz = np.linalg.norm(z)
        if lx == 0 or ly ==0 or lz == 0:
            return np.zeros(6, dtype=np.float32)
        alpha = np.rad2deg(np.arccos(np.dot(y, z) / (ly * lz)))
        beta = np.rad2deg(np.arccos(np.dot(x, z) / (lx * lz)))
        gamma = np.rad2deg(np.arccos(np.dot(x, y) / (lx * ly)))
        box = np.array([lx, ly, lz, alpha, beta, gamma], dtype=np.float32)
        # Only positive edge lengths and angles in (0, 180) are allowed:
        if np.all(box > 0.0) and alpha < 180.0 and beta < 180.0 and gamma < 180.0:
            return box
        # invalid box, return zero vector:
        return np.zeros(6, dtype=np.float32)



    def triclinic_vectors(self, dimensions, dtype = np.float32):
        """Convert ``[lx, ly, lz, alpha, beta, gamma]`` to a triclinic matrix
        representation.

        Original `code by Tsjerk Wassenaar`_ posted on the Gromacs mailinglist.

        If `dimensions` indicates a non-periodic system (i.e., all lengths are
        zero), zero vectors are returned. The same applies for invalid `dimensions`,
        i.e., any box length is zero or negative, or any angle is outside the open
        interval ``(0, 180)``.

        .. _code by Tsjerk Wassenaar:
           http://www.mail-archive.com/gmx-users@gromacs.org/msg28032.html

        Parameters
        ----------
        dimensions : array_like
            Unitcell dimensions provided in the same format as returned by
            :attr:`MDAnalysis.coordinates.timestep.Timestep.dimensions`:\n
            ``[lx, ly, lz, alpha, beta, gamma]``.
        dtype: numpy.dtype
            The data type of the returned box matrix.

        Returns
        -------
        box_matrix : numpy.ndarray
            A numpy array of shape ``(3, 3)`` and dtype `dtype`,
            with ``box_matrix[0]`` containing the first, ``box_matrix[1]`` the
            second, and ``box_matrix[2]`` the third box vector.

        Notes
        -----
        * The first vector is guaranteed to point along the x-axis, i.e., it has the
          form ``(lx, 0, 0)``.
        * The second vector is guaranteed to lie in the x/y-plane, i.e., its
          z-component is guaranteed to be zero.
        * If any box length is negative or zero, or if any box angle is zero, the
          box is treated as invalid and an all-zero-matrix is returned.


        .. versionchanged:: 0.7.6
           Null-vectors are returned for non-periodic (or missing) unit cell.
        .. versionchanged:: 0.20.0
           * Calculations are performed in double precision and zero vectors are
             also returned for invalid boxes.
           * Added optional output dtype parameter.
        """
        dim = np.asarray(dimensions, dtype=np.float64)
        lx, ly, lz, alpha, beta, gamma = dim
        # Only positive edge lengths and angles in (0, 180) are allowed:
        if not (np.all(dim > 0.0) and
                alpha < 180.0 and beta < 180.0 and gamma < 180.0):
            # invalid box, return zero vectors:
            box_matrix = np.zeros((3, 3), dtype=dtype)
        # detect orthogonal boxes:
        elif alpha == beta == gamma == 90.0:
            # box is orthogonal, return a diagonal matrix:
            box_matrix = np.diag(dim[:3].astype(dtype, copy=False))
        # we have a triclinic box:
        else:
            box_matrix = np.zeros((3, 3), dtype=np.float64)
            box_matrix[0, 0] = lx
            # Use exact trigonometric values for right angles:
            if alpha == 90.0:
                cos_alpha = 0.0
            else:
                cos_alpha = np.cos(np.deg2rad(alpha))
            if beta == 90.0:
                cos_beta = 0.0
            else:
                cos_beta = np.cos(np.deg2rad(beta))
            if gamma == 90.0:
                cos_gamma = 0.0
                sin_gamma = 1.0
            else:
                gamma = np.deg2rad(gamma)
                cos_gamma = np.cos(gamma)
                sin_gamma = np.sin(gamma)
            box_matrix[1, 0] = ly * cos_gamma
            box_matrix[1, 1] = ly * sin_gamma
            box_matrix[2, 0] = lz * cos_beta
            box_matrix[2, 1] = lz * (cos_alpha - cos_beta * cos_gamma) / sin_gamma
            box_matrix[2, 2] = np.sqrt(lz * lz - box_matrix[2, 0] ** 2 -
                                       box_matrix[2, 1] ** 2)
            # The discriminant of the above square root is only negative or zero for
            # triplets of box angles that lead to an invalid box (i.e., the sum of
            # any two angles is less than or equal to the third).
            # We don't need to explicitly test for np.nan here since checking for a
            # positive value already covers that.
            if box_matrix[2, 2] > 0.0:
                # all good, convert to correct dtype:
                box_matrix = box_matrix.astype(dtype, copy=False)
            else:
                # invalid box, return zero vectors:
                box_matrix = np.zeros((3, 3), dtype=dtype)
        return box_matrix



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
        u.cell       = u.triclinic_vectors(dimensions)

    return u



def UpdateAtomicNumberAndMassDMS(dms):
    conn = sqlite3.connect(dms)
    conn.execute(update_anum_mass_cmdline)
    conn.commit()
    conn.close()


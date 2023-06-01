import pandas as pd
import numpy  as np
import sqlite3
from   .sqlcmd  import *
from   .readxml import ReadXML
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

class Universe:
    def __init__(self, ifile=None, data=None, add_columns=[], anum=None, create_bonds=False, ff=[], ff_add=[]):
        '''
        create an universe object.
        it will be made through ifile (either dms or pdb) or
        with data provided as a dictionary.
        if none of the above is provided, create an empty universe.
        self.atoms is a pandas dataframe.
        self.cell  will be initialized to a 3x3 matrix with 0.
        self.bonds will be initialized to []. 
        '''
        
        #if not isinstance(add_columns, list): add_columns = [add_columns]
        self.cell       = [[0,0,0], [0,0,0], [0,0,0]] #triclinic_vectors
        self.dimensions = [0, 0, 0, 90, 90, 90]       #triclinic_box
        self.bonds = []
        self.cols  = {'anum': -1, 'name': 'tbd', 'charge': 0.0, 'mass': 0.0,
                      'type': 'tbd', 'nbtype': 0, 'resname': 'tbd', 'resid': 0,
                      'chain': 'A', 'x': 0.0, 'y': 0.0, 'z': 0.0, 'bfactor': 0.0}


        # Martini beads need to have anum; Otherwise, openMM will not run a simulation
        self.anum         = anum
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

        # Update atomic numbers and masses based on names
        self.update_anum_mass()

        self.xml      = ReadXML(ff=ff, ff_add=ff_add)
        self.resnames = set(self.atoms.resname)


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

            # if col not in data.keys():
            #     data[col] = 0
            # elif isinstance(data[col], list) and len(data[col]) == 0:
            #     data[col] = 0
        self.atoms = pd.DataFrame(data)

    def construct_from_df(self, data):
        self.atoms = data.copy().reset_index()
        for col in self.cols:
            if col not in self.atoms.columns:
                self.atoms[col] = self.cols[col]

    def write(self, ofile):
        ext = ofile.split('.')[-1]
        if ext == 'pdb' or ext == 'PDB':
            self.writePDB(ofile)
        elif ext == 'dms' or ext == 'DMS':
            if self.create_bonds: self.addBonds()
            self.writeDMS(ofile)
        else:
            print('the file format should be either pdb or dms')

    def readPDB(self, ifile):
        '''will read only the lines that start with ATOM or CRYST1'''
        data = {'name': [], 'resname': [], 'resid': [], 'chain': [], 
                'x': [], 'y': [], 'z': [], 'bfactor': []}

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

            for index, atom in self.atoms.iterrows():
                name = atom['name']

                # if len(name) = 1,2,3, add space
                if len(name) in [1,2,3]: name = f' {name:s}'
                W.write(fmt['ATOM'].format(
                    self.ltruncate_int(index, 5), name[:4], '', atom.resname[:4],
                    atom.chain[:1], self.ltruncate_int(atom.resid, 4), '',
                    atom.x, atom.y, atom.z, 0.0, atom.bfactor,
                    atom.chain, ''))

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
            self.bonds = pd.read_sql_query("SELECT * FROM BOND", conn).to_numpy()
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

        if self.anum:
            self.atoms.loc[:, 'anum'] = 6

        for index, atom in self.atoms.iterrows():
            # cursor.execute(sql_insert_particle.format(
            #     index, atom.anum, atom['name'], atom.resname, 
            #     atom.chain, atom.resid, atom.mass, atom.charge, 
            #     atom.x, atom.y, atom.z, vx, vy, vz, 
            #     atom.chain, insertion, msys_ct, atom.nbtype, atom.type, atom.bfactor
            #     ))

            cursor.execute(sql_insert_particle, (
                index, atom.anum, atom['name'], atom.resname, 
                atom.chain, atom.resid, atom.mass, atom.charge, 
                atom.x, atom.y, atom.z, vx, vy, vz, 
                atom.chain, insertion, msys_ct, atom.nbtype, atom.type, atom.bfactor
                ))

        ### BONDS
        for bond in self.bonds:
            cursor.execute(sql_insert_bond.format(bond[0], bond[1], 1))

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

    def sort(self, by=['chain', 'resid'], ignore_index=True):
        self.atoms = self.atoms.sort_values(by=by, ignore_index=True)
        self.bonds = []


    def guess_anum_mass(self, name):
        anum = 0
        mass = 0.0
        
        if name.startswith('C'):
            anum = 6
            mass = 12.000

        if name.startswith('N'):
            anum = 7
            mass = 14.000

        if name.startswith('H'):
            anum = 1
            mass = 1.000

        if name.startswith('O'):
            anum = 8
            mass = 16.000

        if name.startswith('S'):
            anum = 16
            mass = 32.000

        if name.startswith('P'):
            anum = 15
            mass = 31.000

        if name.startswith('F'):
            anum = 9
            mass = 19.000

        if name == 'Cl' or name == 'CL':
            anum = 17
            mass = 35.000

        # CG
        if name == 'BB' or name.startswith('SC'):
            anum = 6
            mass = 72.00

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

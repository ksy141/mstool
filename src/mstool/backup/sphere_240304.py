import numpy as np
from   .lipid          import Lipid
from   ..utils.util    import fibo
from   ..core.universe import Universe
import os
import glob
pwd = os.path.dirname(os.path.realpath(__file__))

class Sphere(Lipid):
    def __init__(self, r, dr=0.0, chain='UPPER', monolayer={}, 
                 martini=None, inverse=1, plot=None, out=None,
                 alpha=0.0, beta=0.0, gamma=0.0,
                 lipidpath=pwd + '/../../../FF/martini2.2/structures/'):
        '''Place protein or lipids on a sphere surface with a radius of ``r``.
        The center of the sphere will be approximately (0,0,0).

        Parameters
        ----------
        r : float
            The radius of the sphere.
        dr : float
            The translation amount toward the outer direction.
            +15 is recommended for outer leaflet lipids.
            -15 is recommended for inner leaflet lipids.
            The default is 0.0.
        chain : str or list
            The name of the chain of the molecules.
        monolayer : dict
            Specify what types of molecules and how many of the molecules you want.
            If you are constructing lipids, you should use the lipid resnames provided in
            ``martini.molecules.keys()``.
            e.g., {'POPC': 100, 'DOPC': 100}
        martini : ReadMartini object
            Resnames and atomic names are parsed from here. The default is None.
            e.g., martini = mstool.ReadMartini()
        inverse : float
            1 for the outer leaflet, and -1 for the inner leaflet. The default is 1.
        plot : str
            Provide the name of the plot that shows the points on the Fibonacci sphere.
            The default is None.
        out : str
            Provide the name of the structure file to be saved. The default is None.
        alpha : float
            Angle in degree for yaw to change the initial point as 
            the Fibonacci sphere is deterministic. The default is 0.0.
        beta : float
            Angle in degree for pitch to change the initial point as 
            the Fibonacci sphere is deterministic. The default is 0.0.
        gamma : float
            Angle in degree for pitch to change the initial point as 
            the Fibonacci sphere is deterministic. The default is 0.0.
        lipidpath : str
            Path to a folder that contains the structures of lipids.
            Phospholipids that have names of GL*/C*/D* can be internally constructed.
            However, cholesterols and other molecules are not internally constructed.
            Therefore, users should put lipid structures that cannot be internally constructed to this path.
            The default is $mstool/FF/martini2.2/structures/
            The filename should be ``name_of_moleculetype.pdb`` or ``name_of_moleculetype.dms``.
            The lipid should have an orientation of +z, and approximately centered at the origin.
            Imagine that this lipid is located at the upper leaflet of a bilayer whose center is 0.

        
        Attributes
        ----------
        data : dict
            Contains atomic information.
        universe : Universe
            Contains Universe information.
        atoms : pd.DataFrame
            self.universe.atoms

        Examples
        --------
        >>> import mstool
        >>> from mstool.membrane.sphere import Sphere
        >>> martini = mstool.ReadMartini()
        >>> upper = Sphere(r=100, dr=15, inverse=1, 
        ...                chain='U', monolayer={'POPC':1300, 'DOPC': 1300}, 
        ...                martini=martini, out='upper.pdb')
        >>> lower = Sphere(r=100, dr=-15, inverse=-1,
        ...                chain='L', monolayer={'POPC':700, 'DOPC': 700},
        ...                martini=martini, out='lower.pdb')

        >>> u = mstool.Merge(upper.universe.atoms, lower.universe.atoms)
        >>> u.write('final.pdb')
        '''

        Lipid.__init__(self, martini=martini, lipidpath=lipidpath, hydrophobic_thickness=dr*2)
        monolayerN = np.sum(list(monolayer.values()))
        assert monolayerN > 0, 'Please provide outer and/or inner, e.g., outer = {"POPC": 100}'
        assert r > 0, 'Please provide the radius (r) > 0 A'


        ### monolayer = {'POPC': 3, 'DOPE': 2, 'SAPI': 1, 'protein.pdb': 5}
        ### monolayer_keys = ['POPC', 'DOPE', 'SAPI', Universe('protein.pdb')]
        ### monolayer_list = [0, 0, 0, 1, 1, 2, 3, 3, 3, 3, 3]

        monolayer_keys = []
        for key in monolayer.keys():
            try:
                # protein or a structure file
                monolayer_keys.append(Universe(key))
            except:
                # lipid that needs to be constructed
                monolayer_keys.append(key)

        monolayer_list = []
        for idx, number in enumerate(monolayer.values()):
            monolayer_list += [idx] * number
        

        ### Shuffle so that lipids are not exactly evenly distributed
        np.random.shuffle(monolayer_list)
        

        ### Construct Fibonacchi sphere.
        finalr = r + dr
        points = fibo(r=finalr, N=monolayerN, alpha=alpha, beta=beta, gamma=gamma, plot=plot)
        

        data = {'name': [], 'resname': [], 'resid': [], 'chain': [], 'x': [], 'y': [], 'z': []}
        for i in range(len(points)):
            
            key = monolayer_keys[monolayer_list[i]]
            if isinstance(key, str):
                if key in self.structures.keys():
                    # structure information exists.
                    lipid_u      = self.structures[key]
                    positions    = lipid_u.atoms[['x','y','z']].to_numpy()
                    save_name    = lipid_u.atoms['name'].tolist()
                    save_resname = lipid_u.atoms['resname'].tolist()
                    save_resid   = [i+1]   * len(positions)
                    save_chain   = [chain] * len(positions)


                else:
                    # structure information does not exist. Construct.
                    positions, names = self.construct_molecule(key)
                    save_resname     = [key]   * len(positions)
                    save_resid       = [i+1]   * len(positions)
                    save_chain       = [chain] * len(positions)
                    save_name        = names


            elif isinstance(key, Universe):
                # a structure e.g., protein or lipid with a structure
                positions     = key.atoms[['x','y','z']].to_numpy()
                save_resname  = key.atoms['resname'].tolist()
                save_name     = key.atoms['name']
                save_resid    = key.atoms['resid'].tolist()

                if len(set(save_resid)) == 1:
                    save_resid = [i+1]   * len(positions)
                    save_chain = [chain] * len(positions)

                else:
                    save_chain = chain[i] * len(positions)

            else:
                assert 0 == 1, 'resname does not exist or input structure cannot be found'

            finalpositions = self.place(positions, r=finalr, r_vector=points[i], inverse=inverse)
            data['chain'].extend(save_chain)
            data['resname'].extend(save_resname)
            data['resid'].extend(save_resid)
            data['name'].extend(save_name)
            data['x'].extend(finalpositions[:,0])
            data['y'].extend(finalpositions[:,1])
            data['z'].extend(finalpositions[:,2])
        
        self.universe = Universe(data=data)
        self.universe.sort(by=['chain', 'resname','resid'])
        self.atoms = self.universe.atoms
        if out: self.universe.write(out)


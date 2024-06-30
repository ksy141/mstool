import numpy as np
from   .lipid          import Lipid
from   ..utils.util    import fibo
from   ..core.universe import Universe, Merge, RemoveOverlappedResidues
import os
import glob
pwd = os.path.dirname(os.path.realpath(__file__))

class Sphere(Lipid):
    def __init__(self, r, protein=None, upper={}, lower={}, trans={}, between={}, sep=0.0,
                 water=40, rcut=3.0,
                 martini=None, out=None,
                 hydrophobic_thickness=30.0,
                 alpha=0.0, beta=0.0, gamma=0.0,
                 lipidpath=pwd + '/../../../FF/martini2.2/structures/'):

        '''Make a sphere bilayer (or monolayer) with the provided number of lipids.
        The center of the sphere will be approximately (0,0,0).

        Parameters
        ----------
        r : float
            The radius of the sphere.
        protein : str
            Filename for protein structure. The default is None (lipid-only systems).
            The center of the bilayer will be approxiamtely (0,0,0).
            The protein structure should be prealigned / oriented with respect to this.
        martini : ReadMartini object
            Resnames and atomic names are parsed from here. The default is None.
            e.g., martini = mstool.ReadMartini()
        hydrophobic thickness : float
            Hydrophobic thickness of a bilayer in A.
        water : float
            Thickness of water in A
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
        '''

        Lipid.__init__(self, martini=martini, lipidpath=lipidpath, hydrophobic_thickness=hydrophobic_thickness)

        upperN   = int(np.sum(list(upper.values())))
        lowerN   = int(np.sum(list(lower.values())))
        transN   = int(np.sum(list(trans.values())))
        betweenN = int(np.sum(list(between.values())))
        Nmol     = len(list(upper.values()) + list(lower.values())) + len(trans.values())
        assert Nmol > 0, 'Please provide upper and/or lower, e.g., upper = {"POPC": 100}'

        ### monolayer = {'POPC': 3, 'DOPE': 2, 'SAPI': 1}
        ### monolayer_keys = ['POPC', 'DOPE', 'SAPI']
        ### monolayer_list = [0, 0, 0, 1, 1, 2]

        monolayers     = {'upper': upper, 'lower': lower, 'trans': trans, 'between': between}
        monolayer_keys = {'upper': [], 'lower': [], 'trans': [], 'between': []}
        monolayer_list = {'upper': [], 'lower': [], 'trans': [], 'between': []}

        for layerkey, monolayer in monolayers.items():
            for key in monolayer.keys():
                monolayer_keys[layerkey].append(key)

            for idx, number in enumerate(monolayer.values()):
                monolayer_list[layerkey] += [idx] * number


        ### Shuffle so that lipids are not exactly evenly distributed
        for key in monolayer_list.keys():
            np.random.shuffle(monolayer_list[key])


        ### Construct Fibonacchi sphere.
        hht = hydrophobic_thickness * 0.5

        # upper
        upperU = self.make_sphere(r+hht+sep/2, upperN,
                    monolayer_keys['upper'], monolayer_list['upper'], chain='0',
                    inverse=1, alpha=alpha, beta=beta, gamma=gamma)
                        
        # lower
        lowerU = self.make_sphere(r-hht-sep/2, lowerN,
                    monolayer_keys['lower'], monolayer_list['lower'], chain='1',
                    inverse=-1, alpha=alpha, beta=beta, gamma=gamma)

        # trans
        transU = self.make_sphere(r, transN,
                    monolayer_keys['trans'], monolayer_list['trans'], chain='2',
                    inverse=1, alpha=alpha, beta=beta, gamma=gamma)

        # between
        betweenP = []
        while len(betweenP) < betweenN:
            radius1 = r - sep/2 + 5
            radius2 = r + sep/2 - 5
            # Generate random points within the annular region
            betweenr = np.cbrt(np.random.uniform(radius1**3, radius2**3))
            theta = np.random.uniform(0, 2*np.pi)
            phi = np.arccos(np.random.uniform(-1, 1))
            
            betweenx = betweenr * np.sin(phi) * np.cos(theta)
            betweeny = betweenr * np.sin(phi) * np.sin(theta)
            betweenz = betweenr * np.cos(phi)
            betweenP.append([betweenx, betweeny, betweenz])
            
        betweenU = self.makeUfromPoints(betweenP, 
            monolayer_keys['between'], monolayer_list['between'], chain='3')

        # pbc
        pbc = (r + hht + water + sep/2) * 2

        # protein
        if protein:
            try:
                protein = Universe(protein)
            except:
                protein = protein

        else:
            # construct an empty protein universe
            protein = Universe()
 

        ### New Universe
        #proteinU = Merge(protein.atoms, transU.atoms) if len(transU.atoms) > 0 else protein
        proteinU = Merge(protein.atoms, transU.atoms)
        lipidU   = Merge(upperU.atoms,  lowerU.atoms, betweenU.atoms)
        u        = RemoveOverlappedResidues(lipidU.atoms, proteinU.atoms, rcut=rcut)

        ### Setting pbc
        u.dimensions = np.array([pbc, pbc, pbc, 90, 90, 90])
        u.cell       = np.array([[pbc, 0, 0], [0, pbc, 0], [0, 0, pbc]])

        if out: u.write(out)
        self.universe = u


    def makeUfromPoints(self, points, monolayer_key, monolayer_list, chain):
        data = {'name': [], 'resname': [], 'resid': [], 'chain': [], 'x': [], 'y': [], 'z': []}
        for i in range(len(points)):
            
            key = monolayer_key[monolayer_list[i]]
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
                    save_chain = [chain] * len(positions)

            else:
                assert 0 == 1, 'resname does not exist or input structure cannot be found'
            
            finalr = np.linalg.norm(points[i])
            finalpositions = self.place(positions, r=finalr, r_vector=points[i], inverse=np.random.choice([1, -1]))
            data['chain'].extend(save_chain)
            data['resname'].extend(save_resname)
            data['resid'].extend(save_resid)
            data['name'].extend(save_name)
            data['x'].extend(finalpositions[:,0])
            data['y'].extend(finalpositions[:,1])
            data['z'].extend(finalpositions[:,2])
        
        universe = Universe(data=data)
        # when merging, if there is empty universe, resid changes to float.
        universe.atoms = universe.atoms.astype({'resid': int})
        universe.sort(by=['chain', 'resname','resid'])
        return universe

    def make_sphere(self, finalr, finalN, monolayer_key, monolayer_list, chain, inverse, alpha, beta, gamma):
        if finalN == 0:
            return Universe()
        
        print('leaflet ' + chain)
        points = fibo(r=finalr, N=finalN, alpha=alpha, beta=beta, gamma=gamma, plot=None)
        data = {'name': [], 'resname': [], 'resid': [], 'chain': [], 'x': [], 'y': [], 'z': []}

        for i in range(len(points)):
            
            key = monolayer_key[monolayer_list[i]]
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

        universe = Universe(data=data)
        # when merging, if there is empty universe, resid changes to float.
        universe.atoms = universe.atoms.astype({'resid': int})
        universe.sort(by=['chain', 'resname','resid'])
        return universe


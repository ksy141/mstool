import numpy as np
from   ..utils.util    import fibo, align_a_b
from   ..core.universe import Universe
import os
import glob
pwd = os.path.dirname(os.path.realpath(__file__))

class SphereProtein:
    def __init__(self, radius, protein={}, dr=0.0,
                 alpha=0.0, beta=0.0, gamma=0.0,
                 out=None):

        '''Place proteins at the sphere surface.
        The center of the sphere will be approximately (0,0,0).

        Parameters
        ----------
        radius : float
            The radius of the sphere.
        protein : dict
            Filename in key and number in value.
            The protein structure should be prealigned / oriented 
            with respect to the point [0,0,0] and norm [0,0,1].
        dr : float=0.0
            Systematically shift protein.
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
        '''

        proteinN = int(np.sum(list(protein.values())))
        assert proteinN > 0, 'provide filename and number of each protein'

        ### monolayer = {'POPC': 3, 'DOPE': 2, 'SAPI': 1}
        ### monolayer_keys = ['POPC', 'DOPE', 'SAPI']
        ### monolayer_list = [0, 0, 0, 1, 1, 2]

        monolayer_key  = [key for key in protein.keys()]
        monolayer_list = []
        for idx, number in enumerate(protein.values()):
            monolayer_list += [idx] * number
        np.random.shuffle(monolayer_list)

        # sphere
        proteinU = self.make_sphere(radius+dr, proteinN,
                    monolayer_key, monolayer_list, chain='PROTEIN',
                    alpha=alpha, beta=beta, gamma=gamma)
                        
        if out: proteinU.write(out)
        self.universe = proteinU


    def make_sphere(self, finalr, finalN, monolayer_key, monolayer_list, chain, alpha, beta, gamma):
        points = fibo(r=finalr, N=finalN, alpha=alpha, beta=beta, gamma=gamma, plot=None)
        data = {'name': [], 'resname': [], 'resid': [], 'chain': [], 'x': [], 'y': [], 'z': []}

        for i in range(len(points)):
            
            key = monolayer_key[monolayer_list[i]]
            if isinstance(key, str):
                # structure information exists.
                univ         = Universe(key)
                positions    = univ.atoms[['x','y','z']].to_numpy()
                save_resname = univ.atoms['resname'].tolist()
                save_name    = univ.atoms['name'].tolist()
                save_resid   = univ.atoms['resid'].tolist()
                save_chain   = np.char.add(univ.atoms['chain'].tolist(), [str(i)] * len(positions))

            elif isinstance(key, Universe):
                # a structure e.g., protein or lipid with a structure
                positions     = key.atoms[['x','y','z']].to_numpy()
                save_resname  = key.atoms['resname'].tolist()
                save_name     = key.atoms['name'].tolist()
                save_resid    = key.atoms['resid'].tolist()
                save_chain    = np.char.add(key.atoms['chain'].tolist(), [str(i)] * len(positions))

            else:
                assert 0 == 1, 'input structure cannot be found'

            finalpositions = self.place(positions, r=finalr, r_vector=points[i])
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
        #universe.sort(by=['chain', 'resid'])
        return universe


    def place(self, positions, r=0, r_vector=[0,0,1]):
        positions  = np.array(positions)
        positions += np.array([0, 0, r])
        R          = align_a_b(np.array([0, 0, 1]), r_vector)
        finalpos   = np.matmul(R, positions.T).T
        return finalpos


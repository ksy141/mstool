import numpy as np
from   .lipid import Lipid
from   ..core.universe import Universe, Merge, RemoveOverlappedResidues
import os
import glob
pwd = os.path.dirname(os.path.realpath(__file__))

class Bilayer(Lipid):
    def __init__(self, protein=None, upper={}, lower={}, trans={},
                 dx=8.0, waterz=25.0, rcut=3, out=None, mode='shift',
                 dN=5, martini=None, hydrophobic_thickness=30.0, 
                 lipidpath=pwd + '/../../../FF/martini2.2/structures/'):

        '''Make a plane bilayer (or monolayer) with the provided numbers of lipids.
        The center of the bilyaer will be approximately (0,0,0).

        Parameters
        ----------
        protein : str
            Filename for protein structure. The default is None (lipid-only systems).
            The center of the bilayer will be approxiamtely (0,0,0).
            The protein structure should be prealigned / oriented with respect to this.
        chain : str
            Chain name of the lipid molecules. The defualt is PLANE.
        upper : dict
            Specify what types of molecules and how many of the molecules you want in the upper leaflet.
            If you are constructing lipids, you should use the lipid resnames provided in
            ``martini.molecules.keys()``.
            e.g., {'POPC': 100, 'DOPC': 100}
        lower : dict
            The types (key) and the numbers (value) of molecules in the lower leaflet.
            e.g., {'POPC': 100, 'DOPC': 100}
        trans : dict
            Some molecules span from the upper to lower leaflets. e.g., water pores.
        dx : float
            Initial distance between lipids in A. The default is 8.0.
        waterz : float
            Thickness of water in A. 
            If waterz = 0, the Z dimension of the final structure will be 
            [0, pbcz] where pbcz is the absolute highest Z position of (protein + lipid atoms).
            The default is 25 A, which will increase the Z dimension by 50 compared to waterz = 0
            (+25 A toward +Z and -25 A toward -Z)
        rcut : float
            Cutoff distance in A of removal of lipids that have a close contact with protein (provided via --protein).
        out : str
            Filename for the final structure.
        mode : str
            Options to handle overlap between (protein + trans atoms) and (upper + lower atoms).
            "remove" will remove overlapped lipids.
            "shift" will shift overlapped lipids (no removal of lipids).
            The former is useful when a bilayer is big so that you do not care some removal of lipids.
            The latter is useufl when you want exactly the same number of lipids as you provide.
        dN : int
            Number of additional layers in XY dimensions if you shift overlapped lipids.
            dN = 5 (default) is usually fine, but if you have a protein-crowded membrane structure,
            you should increase this value.
        hydrophobic thickness : float
            Hydrophobic thickness of a bilayer in A.
        lipidpath : str
            Path to a folder that contains the structures of lipids.
            Phospholipids that have names of GL*/C*/D* can be internally constructed.
            However, cholesterols and other molecules are not internally constructed.
            Therefore, users should put lipid structures that cannot be internally constructed to this path.
            The default is $mstool/FF/martini2.2/structures/
            The filename should be ``name_of_moleculetype.pdb`` or ``name_of_moleculetype.dms``.
            The lipid should have an orientation of +z, and approximately centered at the origin.
            Imagine that this lipid is located at the upper leaflet of a bilayer whose center is 0.

        Examples
        --------
        >>> import mstool
        >>> from mstool.membrane.bilayer import Bilayer
        >>> martini = mstool.ReadMartini()
        >>> u = Bilayer(protein='protein.pdb', 
        ...             upper={'POPC':100, 'DOPC':100},
        ...             lower={'POPC':100, 'DOPC':100},
        ...             trans={'P006':1}, martini=martini)
        '''

        Lipid.__init__(self, martini=martini, lipidpath=lipidpath, hydrophobic_thickness=hydrophobic_thickness)

        upperN = int(np.sum(list(upper.values())))
        lowerN = int(np.sum(list(lower.values())))
        transN = int(np.sum(list(trans.values())))
        Nmol   = len(list(upper.values()) + list(lower.values())) + len(trans.values())
        assert Nmol > 0, 'Please provide upper and/or lower, e.g., upper = {"POPC": 100}'



        ### monolayer = {'POPC': 3, 'DOPE': 2, 'SAPI': 1}
        ### monolayer_keys = ['POPC', 'DOPE', 'SAPI']
        ### monolayer_list = [0, 0, 0, 1, 1, 2]

        monolayers     = {'upper': upper, 'lower': lower, 'trans': trans}
        monolayer_keys = {'upper': [], 'lower': [], 'trans': []}
        monolayer_list = {'upper': [], 'lower': [], 'trans': []}

        for layerkey, monolayer in monolayers.items():
            for key in monolayer.keys():
                monolayer_keys[layerkey].append(key)

            for idx, number in enumerate(monolayer.values()):
                monolayer_list[layerkey] += [idx] * number


        ### Shuffle so that lipids are not exactly evenly distributed
        for key in monolayer_list.keys():
            np.random.shuffle(monolayer_list[key])


        ### Construct plane monolayer
        # upper
        upperP, unused_upperP = self.make_rect2(upperN, dx, dN)
        upperU = self.make_monolayer(upperP, 
            monolayer_keys['upper'], monolayer_list['upper'], chain='UPPER', dz=+hydrophobic_thickness/2, inverse=+1.0)

        # lower
        lowerP, unused_lowerP = self.make_rect2(lowerN, dx, dN)
        lowerU = self.make_monolayer(lowerP, 
            monolayer_keys['lower'], monolayer_list['lower'], chain='LOWER', dz=-hydrophobic_thickness/2, inverse=-1.0)

        # pbc
        half_pbcx = max(upperP.max(), abs(upperP.min()), lowerP.max(), abs(lowerP.min()))
        pbcx = 2 * half_pbcx

        # trans
        transP = np.random.rand(transN, 3) - 0.5
        transP[:,2] = 0.0
        transP *= half_pbcx
        transU = self.make_monolayer(transP, 
            monolayer_keys['trans'], monolayer_list['trans'], chain='TRANS', dz=0.0, inverse=+1.0)

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
        proteinU = Merge(protein.atoms, transU.atoms)
        lipidU   = Merge(upperU.atoms,  lowerU.atoms)

       
        if mode == 'remove':
            u = RemoveOverlappedResidues(
                    lipidU.atoms, proteinU.atoms, 
                    rcut=rcut)

        elif mode == 'shift':
            unused_upperP_index = 0
            unused_lowerP_index = 0
            lipidU.addResidues()

            bA = RemoveOverlappedResidues(
                    lipidU.atoms, proteinU.atoms, 
                    rcut=rcut,
                    returnoverlapped=True)
            
            move_atoms = lipidU.atoms[bA]
            move_resns = set(move_atoms['resn'])

            for move_resn in move_resns:
                sel = lipidU.atoms.resn == move_resn
                move_xyz = lipidU.atoms.loc[sel, ['x','y','z']].to_numpy()

                if np.average(move_xyz[:,2]) > 0:
                    #print('UPPER lipid needed to move')
                    target_xyz = unused_upperP[unused_upperP_index] 
                    unused_upperP_index += 1
                else:
                    #print('LOWER lipid needed to move')
                    target_xyz = unused_lowerP[unused_lowerP_index]
                    unused_lowerP_index += 1
                
                move_dx = target_xyz[0] - np.average(move_xyz[:,0])
                move_dy = target_xyz[1] - np.average(move_xyz[:,1])

                lipidU.atoms.loc[sel, ['x','y','z']] = move_xyz + np.array([move_dx, move_dy, 0.0])

            u = Merge(proteinU.atoms, lipidU.atoms)



        ### Setting pbc
        pbcx         = max(u.atoms['x'].max(), u.atoms['y'].max(), abs(u.atoms['x'].min()), abs(u.atoms['y'].min())) * 2 + 3.0
        pbcz         = max(u.atoms['z'].max(), abs(u.atoms['z'].min())) * 2 + waterz * 2
        u.dimensions = np.array([pbcx, pbcx, pbcz, 90, 90, 90])
        u.cell       = np.array([[pbcx, 0, 0], [0, pbcx, 0], [0, 0, pbcz]])

        if out: u.write(out)
        self.universe = u
        self.upperN = upperN
        self.lowerN = lowerN



    def make_rect2(self, N, dx, dN=5):
        '''
        Try to make sqrt(N) x sqrt(N) rectangular points.

          x

         ccc
         cxc
         ccc

        ooooo
        occco
        ocxco
        occco
        ooooo
        '''

        Nx = int(np.sqrt(N)) + dN

        collect = []
        for i in range(Nx):
            pts = np.arange(-i, i+1)

            # edges
            for j in range(1, len(pts)-1):
                collect.append([pts[0], pts[j], 0.0])

            for j in range(1, len(pts)-1):
                collect.append([pts[j], pts[0], 0.0])

            for j in range(1, len(pts)-1):
                collect.append([pts[-1], pts[j], 0.0])

            for j in range(1, len(pts)-1):
                collect.append([pts[j], pts[-1], 0.0])

            # corners
            if len(pts) == 1:
                collect.append([0.0,0.0,0.0])
            
            else:
                collect.append([pts[0], pts[0],   0.0])
                collect.append([pts[0], pts[-1],  0.0])
                collect.append([pts[-1], pts[0],  0.0])
                collect.append([pts[-1], pts[-1], 0.0])

        collect = np.array(collect) * dx
        return collect[0:N], collect[N:]



    def make_monolayer(self, points, monolayer_key, monolayer_list, chain, dz, inverse=1.0):
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

            finalpositions = inverse * np.array(positions) + points[i] + np.array([0,0,dz])
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


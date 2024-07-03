from   ..utils.util import align_a_b
from   ..core.universe import Universe
import numpy as np
import os
import glob
pwd = os.path.dirname(os.path.realpath(__file__))

class Lipid:
    def __init__(self, martini, lipidpath, hydrophobic_thickness):
        '''
        Parameters
        ----------
        martini : mstool.core.ReadMartini
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
        
        #             PC,    PE,    PS,    PG,    PI,               PA    cholesterol
        self.head = ['NC3', 'NH3', 'CNO', 'GL0', 'C1', 'C2', 'C3', 'PO4', 'R1', 'ROH']

        self.chain1 = []
        for i in range(1, 7):
            self.chain1.append('C%dA' %i)
            self.chain1.append('D%dA' %i)

        self.chain2 = []
        for i in range(1, 7):
            self.chain2.append('C%dB' %i)
            self.chain2.append('D%dB' %i)

        self.chain3 = []
        for i in range(1, 7):
            self.chain3.append('C%dC' %i)
            self.chain3.append('D%dC' %i)
        
        if martini: self.MOLS = martini.martini['molecules']


        ### Add lipids that have structure information
        ### structures = {'CHL1': Universe(structure + '/CHL1.pdb'), 'P008': Universe(structure + '/P008.pdb')}  
        structures    = {}
        lipid_ifiles  = []
        for ext in ['*.dms', '*.pdb']:
            lipid_ifiles += glob.glob(lipidpath + '/' + ext)

        for ifile in lipid_ifiles:
            moleculetype = ifile.split('/')[-1].split('.')[0]
            structures[moleculetype] = Universe(ifile)

        self.structures = structures
        self.hydrophobic_thickness = hydrophobic_thickness
        

    def construct_molecule(self, resname):
        names = list(self.MOLS[resname]['atoms']['name'])
        positions = []

        chain1 = list(set(names) & set(self.chain1))
        chain2 = list(set(names) & set(self.chain2))
        chain3 = list(set(names) & set(self.chain3))
        Nmax   = max(len(chain1), len(chain2), len(chain3))
        zz     = np.linspace(0, -(self.hydrophobic_thickness/2 - 2), Nmax + 1)
        
        for name in names:
            added = False

            if name in self.head:
                if name == 'PO4':
                    positions.append([0, 0, 5])
                
                # cholesterol
                elif 'R1' in names:
                    # R1 or ROH
                    if name == 'R1' or name == 'ROH':
                        positions.append(np.random.rand(3) - 0.5 + np.array([0, 0, +5.0]))
                    # C1 or C2
                    else:
                        positions.append(np.random.rand(3) - 0.5 + np.array([0, 0, -5.0]))
                
                # NC3
                else:
                    positions.append([0, 0, 10])
                added = True

            if name == 'GL1':
                positions.append([0, 0, 0])
                added = True

            if name == 'GL2':
                positions.append([2.5, 0, 0])
                added = True

            if name == 'GL3':
                positions.append([2.6, 0, 0])
                added = True

            if name in self.chain1:
                posz = zz[int(name[1])]
                positions.append([0.0, 0, posz])
                added = True

            if name in self.chain2:
                posz = zz[int(name[1])]
                positions.append([2.5, 0, posz])
                added = True

            if name in self.chain3:
                posz = zz[int(name[1])]
                positions.append([2.6, 0, posz])
                added = True

            if not added:
                positions.append(np.random.rand(3) - 0.5)

        return positions, names
    

    def place(self, positions, r=0, r_vector=[0,0,1], inverse=1):
        positions  = inverse * np.array(positions)
        positions += np.array([0, 0, r])
        R          = align_a_b(np.array([0, 0, 1]), r_vector)
        finalpos   = np.matmul(R, positions.T).T
        return finalpos


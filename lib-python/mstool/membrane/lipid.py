from   ..utils.util import align_a_b
import numpy as np

class Lipid:
    def __init__(self, martini=None):

        #             PC,    PE,    PS,    PG,    PI,               PA    cholesterol
        self.head = ['NC3', 'NH3', 'CNO', 'GL0', 'C1', 'C2', 'C3', 'PO4']

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
        

    def construct_molecule(self, resname):

        names = list(self.MOLS[resname]['atoms']['name'])
        positions = []

        chain1 = list(set(names) & set(self.chain1))
        chain2 = list(set(names) & set(self.chain2))
        chain3 = list(set(names) & set(self.chain3))
        Nmax   = max(len(chain1), len(chain2), len(chain3))
        zz     = np.linspace(0, -13, Nmax + 1)
        dz     = zz[0] - zz[1] # e.g., +2.6

        ### Position
        # NC3 = [0,   0, 45]
        # PO4 = [0,   0, 40]
        # GL1 = [0,   0, 35]
        # GL2 = [2.5, 0, 35]
        # GL3 = [5.0, 0, 35]
        
        for name in names:
            added = False

            if name in self.head:
                if name == 'PO4':
                    positions.append([0, 0, 5])
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
                positions.append([5.0, 0, 0])
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
                positions.append([5.0, 0, posz])
                added = True

            if not added:
                positions.append(np.random.rand(3) - 0.5)

        return positions, names
    

#    def place(self, positions, resname, chain, resid, 
#              r=0, r_vector=[0,0,1], inverse=1):
#
#        R          = align_a_b(np.array([0, 0, 1]), r_vector)
#        positions  = inverse * np.array(positions)
#        positions += np.array([0, 0, r])
#        finalpos   = np.matmul(R, positions.T).T
#        n_atoms    = len(names)
#
#        data = {'name':    name,
#                'resname': resname,
#                'resid':   resid,
#                'chain':   chain,
#                'x':       finalpos[:,0],
#                'y':       finalpos[:,1],
#                'z':       finalpos[:,2]}
#        
#        return data


    def place(self, positions, r=0, r_vector=[0,0,1], inverse=1):
        R          = align_a_b(np.array([0, 0, 1]), r_vector)
        positions  = inverse * np.array(positions)
        positions += np.array([0, 0, r])
        finalpos   = np.matmul(R, positions.T).T
        return finalpos


import MDAnalysis as mda
import numpy as np
from .lipid import Lipid
from .solvate_martini import *

class Bilayer(Lipid):

    def __init__(self, protein = None, upper = {}, lower = {},
            waterz = 25, dx = 8.0, fc = 1000, rcut = 3, ff_add=[], structure_add=[]):

        Lipid.__init__(self, ff_add, structure_add)

        upperN = np.sum(list(upper.values()))
        lowerN = np.sum(list(lower.values()))
        Nmol = len(list(upper.values()) + list(lower.values()))
        assert Nmol > 0, 'Please provide upper and/or lower, e.g., upper = {"POPC": 100}'

        try:
            self.protein = mda.Universe(protein)
        except:
            self.protein = protein

        if isinstance(self.protein, mda.Universe):
            self.protein.add_TopologyAttr('segids', ['PROT'] * self.protein.segments.n_segments)

        self.upper = upper
        self.lower = lower
        self.waterz = waterz
        self.dx    = dx
        self.fc    = fc
        self.rcut  = rcut

        backmapped = list(self.upper.keys()) + list(self.lower.keys())

        self.log = open('log.bilayer.txt', 'w')
        self.log.write('waterz: %10.6f\n' %waterz)
        self.log.write('dx:     %10.6f\n' %dx)
        self.log.write('BACKMAPPED: ' + ' '.join(set(backmapped)))
        self.log.write('\n')

        self.universe = self.plane_bilayer()
        self.log.close()


    def make_rect(self, monolayerN):
        '''
        try to make Nx x Nx rectangular points,
        however, stop when the index == monolayerN

        monolayerN <- 9
        Nx <- 4

        procedure:
        1. x

        2. x x
           x x

        3. x x x
           x x x
           x x x


        You shouldn't do:
        Nx <- 4
        make a 4 x 4 grid
        and then fill it up

        x x x x
        x x x x
        x e e e
        e e e e
        '''


        Nx = int(np.sqrt(monolayerN)) + 1
        collect = []
        k = 0

        for i in range(Nx):
            for j in range(i + 1):

                if i != j:
                    for ii, jj in [[i, j], [j, i]]:
                        collect.append([ii, jj])
                        k += 1

                        if k == monolayerN:
                            return collect

                elif i == j:
                    collect.append([i, i])
                    k += 1

                    if k == monolayerN:
                        return collect



    def make_rect_cont(self, prev, addN):
        k  = len(prev)
        curr = self.make_rect(len(prev) + addN)
        new = []
        for i in curr:
            if i not in prev:
                new.append(i)

        return new


    def plane_monolayer(self, monolayer, direction):
        monolayerN = np.sum(list(monolayer.values()))
        Nx = int(np.sqrt(monolayerN)) + 1
        self.xx = self.dx * Nx / 2

        monolayer_keys = list(monolayer.keys())
        monolayer_list = []
        for idx, number in enumerate(monolayer.values()):
            monolayer_list += [idx] * number

        # monolayer = {'POPC': 3, 'DOPE': 2, 'SAPI': 1}
        # monolayer_keys = ['POPC', 'DOPE', 'SAPI']
        # monolayer_list = [0, 0, 0, 1, 1, 2]
        np.random.shuffle(monolayer_list)

        #monolayerU = []; kk = 0; break_flag = False
        #for i in range(Nx):
        #    for j in range(Nx):
        #        if kk == monolayerN:
        #            break_flag = True
        #            break

        #        # molecule information
        #        shift = np.array([-xx + i * self.dx, -xx + j * self.dx, 0])
        #        resname = monolayer_keys[monolayer_list[kk]]
        #        single  = self.construct_molecule(resname, dr=15, r=direction)
        #        single.residues.resids += kk + 1
        #        single.atoms.positions += shift
        #        monolayerU.append(single.atoms)
        #        kk += 1
        #
        #    if break_flag:
        #        break

        monolayerU = []
        rect = self.make_rect(monolayerN)

        for k in range(len(rect)):
            i, j = rect[k]
            shift = np.array([-self.xx + i * self.dx, -self.xx + j * self.dx, 0])
            resname = monolayer_keys[monolayer_list[k]]

            single  = self.construct_molecule(resname, dr=15, r=direction)
            single.atoms.positions += shift
            single.residues.resids = k + 1
            monolayerU.append(single.atoms)

        newU = mda.Merge(*monolayerU)

        if np.sum(direction) == 1:
            newU.add_TopologyAttr('segids', ['UPPER'] * newU.segments.n_segments)

        elif np.sum(direction) == -1:
            newU.add_TopologyAttr('segids', ['LOWER'] * newU.segments.n_segments)

        return newU, rect



    def shallow(self, ag, monolayerR):

        rect_new = self.make_rect_cont(monolayerR, ag.residues.n_residues)

        for k in range(len(rect_new)):
            single = ag.residues[k]
            if single.atoms.select_atoms('name GL1').n_atoms != 0:
                posx, posy, posz = single.atoms.select_atoms('name GL1')[0].position
            elif single.atoms.select_atoms('name R1').n_atoms != 0:
                posx, posy, posz = single.atoms.select_atoms('name R1')[0].position

            single.atoms.positions -= np.array([posx, posy, 0])

            i, j = rect_new[k]
            shift = np.array([-self.xx + i * self.dx, -self.xx + j * self.dx, 0])
            single.atoms.positions += shift

        return ag



    def plane_bilayer(self):

        collect = []
        if len(self.upper) != 0:
            upperU, upperR = self.plane_monolayer(self.upper, direction = [0, 0, +1])
            collect.append(upperU.atoms)

        if len(self.lower) != 0:
            lowerU, lowerR = self.plane_monolayer(self.lower, direction = [0, 0, -1])
            collect.append(lowerU.atoms)

        memb = mda.Merge(*collect)


        # REMOVE lipids close to the protein atoms
        if isinstance(self.protein, mda.Universe):
            u    = mda.Merge(self.protein.atoms, memb.atoms)
            newu = mda.Merge(u.select_atoms('protein or not byres (around %f protein)' %self.rcut))

            if u.select_atoms('byres (around %f protein)' %self.rcut).n_atoms != 0:
                left = mda.Merge(u.select_atoms('byres (around %f protein)' %self.rcut))
                upleft = self.shallow(left.select_atoms('segid UPPER'), upperR)
                loleft = self.shallow(left.select_atoms('segid LOWER'), lowerR)
                self.log.write('MOVED UPPER LEAFLET LIPIDS: %d\n' %(upleft.n_residues))
                self.log.write('MOVED LOWER LEAFLET LIPIDS: %d\n' %(loleft.n_residues))

                if loleft.n_residues == 0:
                    newu = mda.Merge(newu.atoms, upleft)

                elif upleft.n_residues == 0:
                    newu = mda.Merge(newu.atoms, loleft)

                else:
                    newu = mda.Merge(newu.atoms, upleft, loleft)

        else:
            newu = memb


        # RE-ORDER LIPIDS
        resnames = set(self.upper.keys()).union(set(self.lower.keys()))
        #self.make_topol(resnames = resnames, fc = self.fc, resz=True)

        collect = []

        if newu.select_atoms('protein').n_atoms != 0:
            collect.append(newu.select_atoms('protein'))

        for resname in resnames:
            ag = newu.select_atoms('resname ' + resname)
            ag.residues.resids = np.arange(1, ag.residues.n_residues + 1)
            collect.append(newu.select_atoms('resname ' + resname))

        newu = mda.Merge(*collect)

        maxx = max(newu.atoms.positions[:,0])
        minx = min(newu.atoms.positions[:,0])
        maxy = max(newu.atoms.positions[:,1])
        miny = min(newu.atoms.positions[:,1])

        pbcx = (maxx - minx)
        pbcy = (maxy - miny)
        pbcxy = max(pbcx, pbcy)
        pbcz = 70 + self.waterz * 2

        nowater = 5
        newu.dimensions = [pbcxy + nowater, pbcxy + nowater, pbcz, 90, 90, 90]
        # self.shift = np.array([pbcxy, pbcxy, pbcz]) / 2
        # self.log.write('SHIFT: %10.6f %10.6f %10.6f\n'
        #         %(self.shift[0], self.shift[1], self.shift[2]))
        #newu.atoms.positions += self.shift
        newu.atoms.positions += newu.dimensions[0:3] / 2

        s = Solvate_Martini().sol(newu, cutoff=5, zUP = pbcz/2 + 25, zDW = pbcz/2 - 25)
        s.dimensions = newu.dimensions
        return s
        # s2 = Solvate_Martini().ion(s.select_atoms('resname W'), qtot=-16)
        # s2.dimensions = [pbcx, pbcy, pbcz, 90, 90, 90]
        # return s2

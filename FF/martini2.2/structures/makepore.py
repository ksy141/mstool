import numpy as np
import argparse
import mstool
from   mstool.lib.distance import distance_matrix

### an average distance bewteen water atoms in martini is 5.0 A.
### Here, I increased it to 6.0 A, so that waters are not stuck inside the channel.

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('radius', type=float, help='provide a radius of a water cylinder in A')
    parser.add_argument('--height', type=float, default=45.0, help='provide a height of the cylinder in A')
    parser.add_argument('--rcut', type=float, default=10.0, help='provide a cutoff distance of elastic network model in A')
    parser.add_argument('--K', type=float, default=1000.0, help='provide a bond strength for elastic network model in kJ/mol/nm^2')
    parser.add_argument('--fc', type=float, default=1000.0, help='provide a positional restraint strength')
    parser.add_argument('--chain', default='Z', help='Provide a chain name')
    parser.add_argument('--dr', type=float, default=5.0, help='distance between water atoms')
    return parser.parse_args()


def make_structure(radius, height, resname, chain, dr):
    resid   = 1
    data    = {'resname': resname, 'resid': 1, 'chain': chain,
               'name': [], 'x': [], 'y': [], 'z': []}

    Nr = int(2 * np.pi * radius // dr)
    theta = 2 * np.pi / Nr

    Nz = int(height // dr)
    zs = np.linspace(-height/2, height/2, Nz)

    for i in range(Nz):
        for j in range(Nr):
            x = radius * np.cos(theta * j)
            y = radius * np.sin(theta * j)
            z = zs[i]
            name = 'W' + str(Nr * i + j + 1)

            data['name'].append(name)
            data['x'].append(x)
            data['y'].append(y)
            data['z'].append(z)
    
    u = mstool.Universe(data=data)
    u.dimensions = [90] * 6
    u.write(resname + '.pdb')


def enm(u, rcut, K):
    dm = distance_matrix(u.atoms[['x','y','z']].to_numpy(),
                         u.atoms[['x','y','z']].to_numpy())
    
    result = []
    for i in range(len(dm[0])):
        for j in range(i, len(dm[1])):
            if i == j: continue
            if dm[i, j] < rcut:
                result.append([i+1, j+1, dm[i,j]*0.1, K])
    
    return result



def make_itp(bonds, resname, fc):
    itpfile = open(resname + '.itp', 'w')
    itpfile.write('[ moleculetype ]\n')
    itpfile.write('; name  nrexcl\n')
    itpfile.write(f'{resname:s} 1\n\n')

    itpfile.write('[ atoms ]\n')
    n_atoms = len(mstool.Universe(resname + '.pdb').atoms)
    for i in range(1, n_atoms + 1):
        name = f'W{i}'
        itpfile.write(f'{i:10d}  P4  1  {resname:10s} {name:7s} {i:10d} {0.0:10f}\n')
    
    itpfile.write('\n[ bonds ]\n')
    for bond in bonds:
        i1 = bond[0]
        i2 = bond[1]
        l  = bond[2]
        K  = bond[3]
        itpfile.write(f'{i1:10d} {i2:10d}  1  {l:10.3f}  {K:10.3f}\n')
    
    itpfile.write('\n#ifdef POSRES_Z\n')
    itpfile.write('[ position_restraints ]\n')
    for i in range(1, n_atoms + 1):
        itpfile.write(f'{i:10d}  1   {0.0:10.3f}  {0.0:10.3f}  {fc:10.3f}\n')
    itpfile.write('#endif\n')

    itpfile.write('\n#ifdef POSRES_XYZ\n')
    itpfile.write('[ position_restraints ]\n')
    for i in range(1, n_atoms + 1):
        itpfile.write(f'{i:10d}  1   {fc:10.3f}  {fc:10.3f}  {fc:10.3f}\n')
    itpfile.write('#endif\n')

    itpfile.close()


def main():
    args = parse_args()
    resname = 'P{:03d}'.format(int(args.radius))
    
    make_structure(args.radius, args.height, resname, args.chain, args.dr)
    u = mstool.Universe(resname + '.pdb')
    bonds = enm(u, args.rcut, args.K)
    make_itp(bonds, resname, args.fc)

if __name__ == '__main__':
    main()



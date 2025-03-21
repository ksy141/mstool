#!/usr/bin/env python3
import argparse
import mstool
import numpy as np
description = "This script combines several Martini PRO*.itp into one itp"

def parse_args():
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("itp", nargs='*', help="provide itp files")
    parser.add_argument("out", help="the name of the output itp")
    parser.add_argument("-s", "-c", help="provide the structure (required if using elastic network")
    parser.add_argument("-k", default=500.0, type=float, help="a rubber force constant (kJ/nm^2)")
    parser.add_argument("-eu", default=0.9,  type=float, help="a rubber band upper cutoff (nm)")
    parser.add_argument("-el", default=0.4,  type=float, help="a rubber band lower cutoff (nm)")
    return parser.parse_args()

def parse_itp(itp):
    data = {'atoms': [], 'bonds': [], 'constraints': [], 'angles': [], 'dihedrals': [], 'position_restraints': []}
    read = None
    with open(itp) as fin:
        for line in fin:
            line = line.strip()
            if not line: continue
            line = line.split(';')[0]
            if not line: continue
            if line.startswith('#'): continue
            if line.startswith('['):
                read = line.split('[')[1].split(']')[0].strip()
                continue
            if read in data.keys():
                data[read].append(line.split())
    return data

def update_data(data, n_atoms):
    for key in data.keys():
        if key == 'atoms' or key == 'position_restraints':
            for i in range(len(data[key])):
                data[key][i][0] = str(int(data[key][i][0]) + n_atoms)
        if key == 'bonds' or key == 'constraints':
            for i in range(len(data[key])):
                data[key][i][0] = str(int(data[key][i][0]) + n_atoms)
                data[key][i][1] = str(int(data[key][i][1]) + n_atoms)
        if key == 'angles':
            for i in range(len(data[key])):
                data[key][i][0] = str(int(data[key][i][0]) + n_atoms)
                data[key][i][1] = str(int(data[key][i][1]) + n_atoms)
                data[key][i][2] = str(int(data[key][i][2]) + n_atoms)
        if key == 'dihedrals':
            for i in range(len(data[key])):
                data[key][i][0] = str(int(data[key][i][0]) + n_atoms)
                data[key][i][1] = str(int(data[key][i][1]) + n_atoms)
                data[key][i][2] = str(int(data[key][i][2]) + n_atoms)
                data[key][i][3] = str(int(data[key][i][3]) + n_atoms)
    return data

def main():
    args = parse_args()
    n_atoms = 0
    data = {'atoms': [], 'bonds': [], 'constraints': [], 'angles': [], 'dihedrals': [], 'position_restraints': []}
    for itp in args.itp:
        d = parse_itp(itp)
        d = update_data(d, n_atoms)
        for key in data.keys():
            data[key].extend(d[key])
        n_atoms += len(d['atoms'])

    if args.s:
        u = mstool.Universe(args.s)
        atoms = u.select('@BB')
        pos = atoms[['x','y','z']].values
        dist = mstool.distance_matrix(pos, pos, dimensions=None) * 0.1
        rows, cols = np.where((args.el < dist) & (dist < args.eu))
        for row, col in zip(rows, cols):
            if col >= row: continue
            idx1 = atoms.iloc[row].name + 1
            idx2 = atoms.iloc[col].name + 1
            data['bonds'].append([str(idx1), str(idx2), str(1), str(np.round(dist[row, col], 5)), 'RUBBER_FC*1.0'])

    for key, value in data.items():
        save = ''
        for v in value:
            save += ' '.join(v) + '\n'
        if key == 'atoms':
            atoms = save
        if key == 'bonds':
            bonds = save
        if key == 'constraints':
            constraints = save
        if key == 'angles':
            angles = save
        if key == 'dihedrals':
            dihedrals = save
        if key == 'position_restraints':
            position_restraints = save

    with open(args.out, 'w') as W:
        W.write(f'''
[ moleculetype ]
PROTEIN 1

[ atoms ]
{atoms}

[ bonds ]
#ifndef RUBBER_FC
#define RUBBER_FC {args.k}
#endif
{bonds}

[ constraints ]
{constraints}

[ angles ]
{angles}

[ dihedrals ]
{dihedrals}

#ifdef POSRES
#ifndef POSRES_FC
#define POSRES_FC 1000.00
#endif
[ position_restraints ]
{position_restraints}
#endif
''')

if __name__ == '__main__':
    main()



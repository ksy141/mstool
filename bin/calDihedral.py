#!/usr/bin/env python

import argparse
import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

description = """This script collects dihedral angles"""

def cal_Dihedral(pos):
    """pos is Nframe x 4 x 3"""
    dr1 = pos[:,1] - pos[:,0]
    dr2 = pos[:,2] - pos[:,1]
    dr3 = pos[:,3] - pos[:,2]

    vol = np.sum(np.cross(dr1, dr2) * dr3, axis=1)

    v1 = np.cross(dr1, dr2)
    v2 = np.cross(dr2, dr3)

    v1 /= np.linalg.norm(v1, axis=1)[:,None]
    v2 /= np.linalg.norm(v2, axis=1)[:,None]

    x = np.sum(v1 * v2, axis=1)
    # catch roundoffs that lead to nan otherwise
    x = np.arccos(np.clip(x, -1.0, 1.0)) * 180 / np.pi

    bA = vol < 0
    x[bA] *= -1
    return x

def parse_args():
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    
    # Positional arguments
    parser.add_argument("input", help="provide a coarse-grained structure", nargs='*')
    parser.add_argument("sel1",  help="the first selection")
    parser.add_argument("sel2",  help="the first selection")
    parser.add_argument("sel3",  help="the first selection")
    parser.add_argument("sel4",  help="the first selection")

    return parser.parse_args()


def main():
    args = parse_args()

    collect = []

    try:
        u   = mda.Universe(*args.input)
        ag1 = u.select_atoms(args.sel1)
        ag2 = u.select_atoms(args.sel2)
        ag3 = u.select_atoms(args.sel3)
        ag4 = u.select_atoms(args.sel4)

        for ts in u.trajectory:
            pos1 = ag1.positions
            pos2 = ag2.positions
            pos3 = ag3.positions
            pos4 = ag4.positions
            pos  = np.stack([pos1, pos2, pos3, pos4], axis=1)
            collect.extend(list(cal_Dihedral(pos)))
    
    except:
        for ifile in args.input:
            u = mda.Universe(ifile)
            pos1 = u.select_atoms(args.sel1).positions
            pos2 = u.select_atoms(args.sel2).positions
            pos3 = u.select_atoms(args.sel3).positions
            pos4 = u.select_atoms(args.sel4).positions
            pos  = np.stack([pos1, pos2, pos3, pos4], axis=1)
            collect.extend(list(cal_Dihedral(pos)))
    print(collect)

    fig, ax = plt.subplots()
    ax.hist(collect, bins=np.linspace(-180.0, 180.0, 20))
    ax.axvline(np.average(collect))
    ax.set_xlabel('dihedral (degree)')
    ax.set_xlim(-180.0, 180.0)
    fig.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()


import os
import mstool
import argparse
import shutil
import numpy as np
from   openmm.app import *
from   openmm import *
from   openmm.unit import *


### ARGPARSE
# * for several arguments and ? for only one argument
# nargs=*:                   --ff FIRST       SECOND ---> [FIRST, SECOND]
# nargs=*:                   --ff FIRST  --ff SECOND ---> [SECOND]
# nargs=* and action=append: --ff FIRST  --ff SECOND ---> [[FIRST], [SECOND]]

# nargs=?:                   --ff FIRST              ---> FIRST
# nargs=?:                   --ff FIRST  --ff SECOND ---> SECOND
# nargs=? and action=append: --ff FIRST  --ff SECOND ---> [FIRST, SECOND]


description = """This script predicts structures of missing loops"""

def parse_args():
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--protein', nargs='?', required=True,     help='provide a protein structure (pdb/dms)')
    parser.add_argument('--fasta',   nargs='?', required=True,     help='provide a fasta file. A single fasta file can contain multiple chains')
    parser.add_argument('--workdir', nargs='?', default='workdir', help='provide a workdir path')
    parser.add_argument('--A',       nargs='?', type=float,        default=100.0,     help='(soft) repulsion parameter')
    parser.add_argument('--C',       nargs='?', type=float,        default=50.0,      help='(soft) Coulomb parameter')
    parser.add_argument('--fc1',     nargs='?', type=float,        default=50.0,      help='posre strength (kcal/mol/A^2)')
    parser.add_argument('--fc2',     nargs='?', type=float,        default=1000.0,    help='steerMD strength (kJ/mol/nm^2)')
    parser.add_argument('--Kchiral', nargs='?', type=float,        default=300.0,     help='Kchiral')
    parser.add_argument('--Kpeptide',nargs='?', type=float,        default=300.0,     help='Kpeptide')

    # turn on soft interactions by default 
    parser.add_argument('--soft',    action='store_true',  dest='soft')
    parser.add_argument('--no-soft', action='store_false', dest='soft', 
                        help='by default, the tool uses soft EM. Use this option to turn it off')
    parser.set_defaults(soft=True)
    
    # turn on fix mutation by default
    parser.add_argument('--mutate',    action='store_true',  dest='mutate')
    parser.add_argument('--no-mutate', action='store_false', dest='mutate', 
                        help='by default, the tool fixes the mutation. Use this option to turn it off')
    parser.set_defaults(mutate=True)

    # Optional arguments
    parser.add_argument("--mapping",
        nargs   = "*",
        default = [],
        help    = "mapping information of molecules. defaults: $mstool/mapping/martini.protein.c36m.dat and $mstool/mapping/martini.lipid.c36.dat")

    parser.add_argument("--mapping_add",
        nargs   = "*",
        default = [],
        help    = "additional mapping information of molecules")
    
    parser.add_argument("--ff",
        nargs   = "*",
        default = [],
        help    = "all-atom forcefield. defaults: $mstool/FF/charmm36/charmm36.xml and $mstool/FF/charmm36/water.xml")

    parser.add_argument("--ff_add",
        nargs   = "*",
        default = [],
        help    = "additional all-atom forcefield")

    return parser.parse_args()


def main():
    args = parse_args()
    mstool.LoopModeler(protein=args.protein,
                       fasta=args.fasta,
                       workdir=args.workdir,
                       A=args.A,
                       C=args.C,
                       soft=args.soft,
                       mutate=args.mutate,
                       mapping=args.mapping,
                       mapping_add=args.mapping_add,
                       ff=args.ff,
                       ff_add=args.ff_add,
                       fc1=args.fc1,
                       fc2=args.fc2,
                       Kchiral=args.Kchiral,
                       Kpeptide=args.Kpeptide)


if __name__ == '__main__':
    main()


import argparse
import mstool

# Execute from the Python interpreter.
# run the code inside the if statement when the program is run directly by the Python interpreter. 
# The code inside the if statement is not executed when the file's code is imported as a module.

# * for several arguments and ? for only one argument
# nargs=*:                   --ff FIRST       SECOND ---> [FIRST, SECOND]
# nargs=*:                   --ff FIRST  --ff SECOND ---> [SECOND]
# nargs=* and action=append: --ff FIRST  --ff SECOND ---> [[FIRST], [SECOND]]

# nargs=?:                   --ff FIRST              ---> FIRST
# nargs=?:                   --ff FIRST  --ff SECOND ---> SECOND
# nargs=? and action=append: --ff FIRST  --ff SECOND ---> [FIRST, SECOND]

description = """This script backmaps a coarse-grained structure into an atomistic structure"""

def parse_args():
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    
    # Positional arguments
    parser.add_argument("input",   help="provide a coarse-grained structure")
    parser.add_argument("output",  help="provide a name of an all-atom structure")

    # Optional arguments
    parser.add_argument("--ref",   help="reference structure")
    parser.add_argument("--mapping",
        nargs   = "*",
        default = [],
        help    = "cis/trans/chiral information of molecules. defaults: $mstool/mapping/martini.protein.c36m.dat and $mstool/mapping/martini.lipid.c36.dat")

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

    # True / False
    parser.add_argument('--backbone',    action='store_true',  dest='backbone',
        help = 'Pre-construct protein backbone atoms using CG2AA algorithm by Lombardi1, Marti, Capece (default)')
    parser.add_argument('--no-backbone', action='store_false', dest='backbone', 
        help='turn off pre-constructing protein backbone atoms using CG2AA algorithm')
    parser.set_defaults(backbone=True)

    return parser.parse_args()


def main():
    args = parse_args()
    mstool.Backmap(structure=args.input, out=args.output, refstructure=args.ref, 
        mapping=args.mapping, mapping_add=args.mapping_add,
        ff=args.ff, ff_add=args.ff_add, backbone=args.backbone)

if __name__ == '__main__':
    main()
    

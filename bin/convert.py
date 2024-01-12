import argparse
import mstool

description = "Change a file format between pdb and dms"

def parse_args():
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)

    # Positional arguments
    parser.add_argument("input",  help="provide an input structure")
    parser.add_argument("output", help="provide an output structure")

    # Optional arguments
    parser.add_argument("--ff",
        nargs   = "*",
        default = [],
        help    = "all-atom forcefield. defaults: $mstool/FF/charmm36/charmm36.xml and $mstool/FF/charmm36/water.xml")

    parser.add_argument("--ff_add",
        nargs   = "*",
        default = [],
        help    = "additional all-atom forcefield")

    # True / False
    parser.add_argument("--create-bonds",    action='store_true',  dest='bonds', help='create bonds before making a DMS file')
    parser.add_argument("--no-create-bonds", action='store_false', dest='bonds', help='do not create bonds (default)') 
    parser.set_defaults(bonds=True)

    return parser.parse_args()

def main():
    args = parse_args()
    u = mstool.Universe(args.input, ff=args.ff, ff_add=args.ff_add, create_bonds=args.bonds)
    u.write(args.output)

if __name__ == '__main__':
    main()


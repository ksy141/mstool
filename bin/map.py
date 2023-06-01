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

description = """This script maps an atomistic structure into a coarse-grained structure"""

def parse_args():
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)

    # Positional arguments
    parser.add_argument("input",  help="provide an all-atom structure")
    parser.add_argument("output", help="provide a name of coarse-grained structure")

    # Optional arguments
    parser.add_argument("--mapping",
        nargs   = "*",
        default = [],
        help    = "mapping information of molecules. defaults: $mstool/mapping/martini.protein.c36m.dat and $mstool/mapping/martini.lipid.c36.dat")

    parser.add_argument("--mapping_add",
        nargs   = "*",
        default = [],
        help    = "additional mapping information of molecules")

    # True / False
    parser.add_argument('--BB2CA',    action='store_true',  dest='BB2CA')
    parser.add_argument('--no-BB2CA', action='store_false', dest='BB2CA', 
        help='The tool uses the positions of CA atoms for BB atoms by default. Use this option to turn it off')
    parser.set_defaults(BB2CA=True)

    return parser.parse_args()

def main():
    args = parse_args()
    mstool.Map(structure=args.input, out=args.output, mapping=args.mapping, mapping_add=args.mapping_add, BB2CA=args.BB2CA)


if __name__ == '__main__':
    main()

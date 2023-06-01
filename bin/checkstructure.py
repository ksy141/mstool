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

description = """This script checks whether a structure has a flipped isomer"""

def parse_args():
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)

    # Positional arguments
    parser.add_argument("input", help="provide an all-atom structure")
    
    # Optional arguments
    parser.add_argument("--mapping",
        nargs   = "*",
        default = [],
        help    = "cis/trans/chiral information of molecules. defaults:  $mstool/mapping/martini.protein.c36m.dat and $mstool/mapping/martini.lipid.c36.dat")

    parser.add_argument("--mapping_add",
        nargs   = "*",
        default = [],
        help    = "additional mapping information of molecules")

    parser.add_argument("--log",
        default=None,
        help="provide the name of a log file. If not provided, the file will not be created.")

    return parser.parse_args()

def main():
    args = parse_args()
    mstool.CheckStructure(structure=args.input, mapping=args.mapping, mapping_add=args.mapping_add, log=args.log)

if __name__ == '__main__':
    main()

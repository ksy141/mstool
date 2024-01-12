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

description = """This script minimizes a structure using reduced energy minimization (REM)"""

def parse_args():
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    
    # Positional arguments
    parser.add_argument("input",   help="provide a coarse-grained structure")
    parser.add_argument("output",  help="provide a name of an all-atom structure")

    # Optional arguments
    parser.add_argument("--protein",       default=None,   help="protein structure (openMM will add bonds)")
    parser.add_arugment("--rock",          default=None,   help="molecule structure that should not move during REM (called ROCK)")
    parser.add_arugment("--rockrcut",      default=1.2,    type=float, help="cutoff distance (nm) for elastic network model for ROCK molecules")
    parser.add_arugment("--rockKbond",     default=None,   type=float, help="bond strenght (kJ/mol/nm^2) for elastic network model for ROCK molecules")

    parser.add_argument("--rcut",          default=1.2,    type=float, help="cutoff distance (nm)")
    parser.add_argument("--A",             default=100.0,  type=float, help="repulsion parameter (kJ/mol)")
    parser.add_argument("--C",             default=50.0,   type=float, help="Coulomb parameter (kJ/mol)")
    parser.add_argument("--Kchiral",       default=300.0,  type=float, help="Chiral torsion parameter (kJ/mol)")
    parser.add_argument("--Kpeptide",      default=300.0,  type=float, help="Peptide torsion parameter (kJ/mol)")
    parser.add_argument("--Kcistrans",     default=300.0,  type=float, help="Peptide torsion parameter (kJ/mol)")
    
    parser.add_argument("--bfactor_posre", default=0.5,    type=float, help="atoms whose bfactors are larger than this will have positional restraints")
    parser.add_argument("--fcx",           default=1000.0, type=float, help="positional restraints for x (kJ/mol/nm^2)")
    parser.add_argument("--fcy",           default=1000.0, type=float, help="positional restraints for y (kJ/mol/nm^2)")
    parser.add_argument("--fcz",           default=1000.0, type=float, help="positional restraints for z (kJ/mol/nm^2)")
    parser.add_argument("--refposre",      default=None,   type=float, help="positional restraints if atoms exist in this structure")


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
    parser.add_argument('--pbc',    action='store_true',  dest='pbc',
        help = 'apply PBC (default)')
    parser.add_argument('--no-pbc', action='store_false', dest='pbc',
        help = 'PBC is not taken into account (useful for loop prediction)')
    parser.set_defaults(pbc=True)

    return parser.parse_args()

def main():
    args = parse_args()
    mstool.REM(structure=args.input, out=args.output, protein=args.protein, refposre=args.refposre,
        rock=args.rock, rockrcut=args.rockrcut, rockKbond=args.rockKbond, rockresname='ROCK', 
        rcut=args.rcut, pbc=args.pbc,
        A=args.A, C=args.C,
        mapping=args.mapping, mapping_add=args.mapping_add,
        ff=args.ff, ff_add=args.ff_add,
        Kchiral=args.Kchiral, Kpeptide=args.Kpeptide, Kcistrans=args.Kcistrans,
        fcx=args.fcx, fcy=args.fcy, fcz=args.fcz, bfactor_posre=args.bfactor_posre)

    mstool.CheckStructure(structure=args.output, mapping=args.mapping, mapping_add=args.mapping_add)

if __name__ == '__main__':
    main()

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
    parser.add_argument('--ref',     nargs='?', default='None',    help='provide a reference structure (this is for a test purpose)')
    parser.add_argument('--workdir', nargs='?', default='workdir', help='provide a workdir path')
    parser.add_argument('--A',       nargs='?', default=100.0,     help='(soft) repulsion parameter')
    parser.add_argument('--C',       nargs='?', default=50.0,      help='(soft) Coulomb parameter')
    
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
    
    return parser.parse_args()


def rmsd(a, b):
    N = len(a)
    return np.sqrt(np.sum((a - b) ** 2) / N)

def compare_rmsd(refcg, u, resids):
    bA1ref = refcg.atoms.resid.isin(resids)
    bA2ref = refcg.atoms.name == 'BB'

    bA1cg  = u.atoms.resid.isin(resids)
    bA2cg  = u.atoms.name == 'BB'

    refxyz   = refcg.atoms[ bA1ref ][['x','y','z']].to_numpy()
    modelxyz = u.atoms[     bA1cg  ][['x','y','z']].to_numpy()
    print('RMSD of modeled BB+SC*: {:.2f}'.format(rmsd(refxyz, modelxyz)))

    refxyz   = refcg.atoms[ bA1ref & bA2ref ][['x','y','z']].to_numpy()
    modelxyz = u.atoms[     bA1cg  & bA2cg  ][['x','y','z']].to_numpy()
    print('RMSD of modeled BB+SC*: {:.2f}'.format(rmsd(refxyz, modelxyz)))


def main():
    args = parse_args()

    ### Make a workdir
    if os.path.exists(args.workdir):
        shutil.rmtree(args.workdir)
        #raise Exception(args.workdir + ' already exists')
    os.mkdir(args.workdir)
    
    ### Read and save a protein structure
    u = mstool.Universe(args.protein)
    u.atoms.bfactor = 1.0
    u.write(args.workdir + '/step1_input.dms')
    u.write(args.workdir + '/step1_input.pdb')
    
    u   = mstool.Universe(args.workdir + '/step1_input.dms')
    seq = mstool.Seq(fasta=args.fasta)
    
    ### Fix mutations
    if args.mutate:
        u = mstool.Mutate(u, seq)
        u.write(args.workdir + '/step1_mutate.dms')
        u.write(args.workdir + '/step1_mutate.pdb')
    else:
        u.write(args.workdir + '/step1_mutate.dms')
        u.write(args.workdir + '/step1_mutate.pdb')
    u = mstool.Universe(args.workdir + '/step1_mutate.dms')
    
    ### Map 
    cg = mstool.Map(structure   = args.workdir + '/step1_mutate.dms', 
                    mapping     = args.mapping, 
                    mapping_add = args.mapping_add,
                    BB2CA       = True,
                    out         = args.workdir + '/step2.dms')
    
    
    ### Map - reference
    if args.ref != 'None':
        refcg = mstool.Map(structure = args.ref,
                           mapping = args.mapping,
                           mapping_add = args.mapping_add,
                           BB2CA = True,
                           out = args.workdir + '/ref.cg.dms')
        resids = set(refcg.atoms.resid) - set(cg.atoms.resid)
        print('-------------------------------')
        print('modeling the following residues:')
        print(sorted(resids))
    

    ### Fill
    newcg = mstool.Fill(structure   = args.workdir + '/step2.dms',
                        sequence    = args.fasta,
                        mapping     = args.mapping,
                        mapping_add = args.mapping_add,
                        out         = args.workdir + '/step3.dms')

    
    ### Martini
    martini = mstool.ReadMartini()
    mstool.MartinizeDMS(args.workdir + '/step3.dms', martini=martini, output=args.workdir + '/step3.ff.dms')
    mstool.dumpsql(args.workdir + '/step3.ff.dms')
    
    ### From now on, step3.ff.dms is the reference dms
    u = mstool.Universe(args.workdir + '/step3.ff.dms')
    
    
    # --------------------------------------------------------------------------------------------------#
    ### openMM system
    system, dms = mstool.DMS2openmm(
            dms = args.workdir + '/step3.ff.dms',
            nonbondedMethod = 'CutoffNonPeriodic',
            soft = args.soft,
            A    = args.A,
            C    = args.C).make()
    
    ### RUN SOFT SIMS
    integrator = LangevinMiddleIntegrator(310*kelvin, 1/picosecond, 0.002*picoseconds)
    simulation = Simulation(dms.topology, system, integrator)
    simulation.context.setPositions(dms.positions)
    
    print('-------------------------------')
    if args.soft:
        print('Soft interactions are turned on')
    else:
        print('Soft interactions are turned off')
    
    print('Initial: {:+.2f} kJ/mol'.format(mstool.getEnergy(simulation)))
    simulation.minimizeEnergy()
    print('EM:      {:+.2f} kJ/mol'.format(mstool.getEnergy(simulation)))
    
    numpypositions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True) * 10
    u.setpositions(numpypositions)
    u.write(args.workdir + '/step4.dms')
    u.write(args.workdir + '/step4.pdb')
    mstool.dumpsql(args.workdir + '/step4.dms')
    
    if args.ref != 'None':
    
        bA1ref = refcg.atoms.resid.isin(resids)
        bA2ref = refcg.atoms.name == 'BB'
    
        bA1cg  = u.atoms.resid.isin(resids)
        bA2cg  = u.atoms.name == 'BB'
    
        refxyz   = refcg.atoms[ bA1ref ][['x','y','z']].to_numpy()
        modelxyz = u.atoms[     bA1cg  ][['x','y','z']].to_numpy()
        print('RMSD of modeled BB+SC*: {:.2f}'.format(rmsd(refxyz, modelxyz)))
    
        refxyz   = refcg.atoms[ bA1ref & bA2ref ][['x','y','z']].to_numpy()
        modelxyz = u.atoms[     bA1cg  & bA2cg  ][['x','y','z']].to_numpy()
        print('RMSD of modeled BB+SC*: {:.2f}'.format(rmsd(refxyz, modelxyz)))
    
    
    # --------------------------------------------------------------------------------------------------#
    ### openMM system
    system, dms = mstool.DMS2openmm(
            dms = args.workdir + '/step3.ff.dms',
            nonbondedMethod = 'CutoffNonPeriodic',
            soft = False).make()
    
    ### Use previous simulation's positions
    integrator2 = LangevinMiddleIntegrator(310*kelvin, 1/picosecond, 0.002*picoseconds)
    simulation2 = Simulation(dms.topology, system, integrator2)
    simulation2.context.setPositions(simulation.context.getState(getPositions=True).getPositions())
    
    print('-------------------------------')
    print('Regular Energy Minimization')
    print('Initial: {:+.2f} kJ/mol'.format(mstool.getEnergy(simulation2)))
    simulation2.minimizeEnergy()
    print('EM:      {:+.2f} kJ/mol'.format(mstool.getEnergy(simulation2)))
    numpypositions = simulation2.context.getState(getPositions=True).getPositions(asNumpy=True) * 10
    u.setpositions(numpypositions)
    u.write(args.workdir + '/step5.dms')
    u.write(args.workdir + '/step5.pdb')
    if args.ref != 'None':
        compare_rmsd(refcg, u, resids)
    
    
    print('-------------------------------')
    print('Short NVT')
    print('Initial: {:+.2f} kJ/mol'.format(mstool.getEnergy(simulation2)))
    simulation2.reporters.append(PDBReporter(args.workdir + '/step6.trj.pdb', 10))
    simulation2.reporters.append(StateDataReporter(args.workdir + '/step6.csv', 10, 
            step=True, potentialEnergy=True, temperature=True))
    simulation2.step(1000)
    print('EM:      {:+.2f} kJ/mol'.format(mstool.getEnergy(simulation2)))
    
    numpypositions = simulation2.context.getState(getPositions=True).getPositions(asNumpy=True) * 10
    u.setpositions(numpypositions)
    u.write(args.workdir + '/step6.dms')
    u.write(args.workdir + '/step6.pdb')
    u.write(args.workdir + '/cg.final.dms')
    u.write(args.workdir + '/cg.final.pdb')
    if args.ref != 'None':
        compare_rmsd(refcg, u, resids)



    # --------------------------------------------------------------------------------------------------#
    ### Backmap


if __name__ == '__main__':
    main()


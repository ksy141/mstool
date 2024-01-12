import os
import argparse
import shutil
import numpy  as np
import pandas as pd

from   .mutate                import Mutate
from   .fill                  import Fill
from   .readmappings          import ReadMappings
from   .universe              import Universe
from   .map                   import Map
from   .readmartini           import ReadMartini
from   .martinizedms          import MartinizeDMS
from   .dms2openmm            import DMS2openmm
from   .rem                   import REM
from   .checkstructure        import CheckStructure
from   .backmap               import Backmap

from   ..utils.protein_sel    import three2one
from   ..utils.dump           import dumpsql
from   ..utils.openmmutils    import getEnergy

from   openmm.app             import *
from   openmm                 import *
from   openmm.unit            import *

class LoopModeler:
    def __init__(self, protein, fasta, ref=None, workdir='workdir', 
        A=100, C=50, soft=True, mutate=True, 
        mapping=[], mapping_add=[], ff=[], ff_add=[],
        fc1 = 50, fc2 = 1000, Kchiral=300.0, Kpeptide=300.0):

        self.protein     = protein
        self.fasta       = fasta
        self.workdir     = workdir
        self.A           = A
        self.C           = C
        self.soft        = soft
        self.mutate      = mutate
        self.mapping     = mapping
        self.mapping_add = mapping_add
        self.fc1         = fc1


        ### step1: workdir
        self.create_workdir()


        ### step2: mutate
        self.fix_mutation(workdir + '/step1_input.dms',
                          workdir + '/step2_mutate.dms')


        ### step2.1: remove dangling CONH
        self.remove_dangling_CONH(workdir + '/step2_mutate.dms',
                                  workdir + '/step2_mutate_new.dms')


        ### step3: map
        Map(structure   = workdir + '/step2_mutate.dms',
            out         = workdir + '/step3_cg.dms',
            mapping     = mapping,
            mapping_add = mapping_add,
            BB2CA       = True)


        ### step4: fill martini loops
        Fill(structure   = workdir + '/step3_cg.dms',
             out         = workdir + '/step4_cgfilled.dms',
             sequence    = fasta,
             mapping     = mapping,
             mapping_add = mapping_add)


        ### step5: create a martini FF dms
        self.makemartini()


        ### steo6: run Martini EM
        self.runCG(input  = workdir + '/step5_ff.dms', 
                   output = workdir + '/step6_minimized.dms',
                   soft   = soft)


        ### step7: backmap (output must be a pdb so that openMM recognizes protein residues)
        Backmap(structure   = workdir + '/step6_minimized.dms',
                out         = workdir + '/step7_backmapped.pdb',
                mapping     = mapping,
                mapping_add = mapping_add, 
                ff          = ff,
                ff_add      = ff_add,
                backbone    = True)


        ### step8: REM + steer MD
        REM(structure   = workdir + '/step7_backmapped.pdb',
            out         = workdir + '/step8_minimized.pdb',
            refposre    = workdir + '/step2_mutate_new.dms',
            pbc         = False,
            A           = A,
            C           = C,
            mapping     = mapping,
            mapping_add = mapping_add,
            ff          = ff,
            ff_add      = ff_add,
            fcx         = fc2,
            fcy         = fc2,
            fcz         = fc2,
            Kchiral     = Kchiral,
            Kpeptide    = Kpeptide,
            Kcistrans   = 0.0)


        ### check structure
        CheckStructure(structure   = workdir + '/step8_minimized.pdb',
                       mapping     = mapping,
                       mapping_add = mapping_add)


        ### step9: combine
        u = self.combine_two(xtal    = workdir + '/step1_input.dms', 
                             backmap = workdir + '/step8_minimized.pdb')
        u.write(workdir + '/step9_final.dms')
        u.write(workdir + '/step9_final.pdb')


        ### vis
        self.vis()



    def create_workdir(self):
        if os.path.exists(self.workdir):
            #shutil.rmtree(self.workdir)
            raise Exception(self.workdir + ' already exists')
        os.mkdir(self.workdir)

        ### Read and save a protein structure
        u = Universe(self.protein)
        u.atoms.bfactor = 1.0
        u.write(self.workdir + '/step1_input.dms')
        u.write(self.workdir + '/step1_input.pdb')


    def fix_mutation(self, structure, out):
        if self.mutate:
            Mutate(structure = structure,
                   fasta     = self.fasta,
                   out       = out)
        else:
            u = Universe(structure).write(out)


    def map2martini(self):
        Map(structure   = self.workdir + '/step2_mutate.dms',
            mapping     = self.mapping,
            mapping_add = self.mapping_add,
            BB2CA       = True,
            out         = self.workdir + '/step3_cg.dms')


    def fillmartini(self):
        Fill(structure   = self.workdir + '/step3_cg.dms',
             sequence    = self.fasta,
             mapping     = self.mapping,
             mapping_add = self.mapping_add,
             out         = self.workdir + '/step4_cgfilled.dms')


    def makemartini(self):
        martini = ReadMartini()
        # dms default unit: kcal/mol, A^
        # 50 kcal/mol/A^2 -> 0.5 * 50 * 4.184 * 100 kJ/mol/nm^2
        MartinizeDMS(dms     = self.workdir + '/step4_cgfilled.dms',
                     martini = martini, 
                     output  = self.workdir + '/step5_ff.dms',
                     fcx     = self.fc1,
                     fcy     = self.fc1,
                     fcz     = self.fc1)

        dumpsql(      self.workdir + '/step5_ff.dms')
        u = Universe( self.workdir + '/step5_ff.dms')
        u.write(      self.workdir + '/step5_ff.pdb')


    def runCG(self, input, output, soft):

        ### step5_ff.dms is a reference and contains all the FF information
        ### please note that you steal positions from the input
        system, dms = DMS2openmm(
                dms_in          = self.workdir + '/step5_ff.dms',
                nonbondedMethod = 'CutoffNonPeriodic',
                soft            = soft,
                A               = self.A,
                C               = self.C).make()
        

        ### RUN SIMS
        integrator = LangevinMiddleIntegrator(310*kelvin, 1/picosecond, 0.002*picoseconds)
        simulation = Simulation(dms.topology, system, integrator)
        simulation.context.setPositions(DesmondDMSFile(input).positions)
        
        print('-------------------------------')
        if soft:
            print('Soft interactions are turned on')
        else:
            print('Soft interactions are turned off')
        
        print('E0: %.3e kJ/mol' %getEnergy(simulation))
        simulation.minimizeEnergy()
        print('E1: %.3e kJ/mol' %getEnergy(simulation))
        print('-------------------------------')
        
        u = Universe(self.workdir + '/step5_ff.dms')
        numpypositions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True) * 10
        u.setpositions(numpypositions)
        u.write(output)


    def check_structure(self, structure):
        CheckStructure(structure   = structure,
                       mapping     = self.mapping,
                       mapping_add = self.mapping_add)


    def combine_two(self, xtal, backmap):
        exclude = ['ACE', 'NMA', 'NME'] + list(three2one.keys())
        xtal    = Universe(xtal)
        backmap = Universe(backmap)
        
        ### select non-protein residues
        df = xtal.atoms[~xtal.atoms['resname'].isin(exclude)]

        ### combine
        backmap.atoms = pd.concat([backmap.atoms, df], ignore_index=True)
        backmap.sort()

        ### copy dimension
        backmap.cell       = xtal.cell
        backmap.dimensions = xtal.dimensions

        ### atoms exist in xtal have bfactor of 1.0
        bfactors = []
        for index, atom in backmap.atoms.iterrows():
            bA1 = atom['name']  == xtal.atoms.name
            bA2 = atom['resid'] == xtal.atoms.resid
            bA3 = atom['chain'] == xtal.atoms.chain

            xtalatoms = xtal.atoms[bA1 & bA2 & bA3]
            if len(xtalatoms) == 0:
                bfactors.append(0.0)
                #atom.bfactor = 0.0

            elif len(xtalatoms) == 1:
                bfactors.append(1.0)
                #atom.bfactor = 1.0

            else:
                assert 0 == 1, 'more than one atom with the same name, resid, chain?'

        backmap.atoms.bfactor = bfactors
        return backmap


    def remove_dangling_CONH(self, input, output):
        u = Universe(input)

        saved = []
        for index, atom in u.atoms.iterrows():
            resid   = atom.resid
            chain   = atom.chain
            resname = atom.resname
            record  = chain + ':' + str(resid)

            if resname not in three2one.keys():
                continue
                
            if record in saved:
                continue

            bA1 = u.atoms.chain == chain
            bA2 = u.atoms.resid == resid + 1
            bA3 = u.atoms.resid == resid - 1
            
            if len(u.atoms[bA1 & bA2]) > 0 and len(u.atoms[bA1 & bA3]) > 0:
                saved.append(record)
                
            elif len(u.atoms[bA1 & bA2]) == 0 and len(u.atoms[bA1 & bA3]) > 0:
                bA4 = u.atoms.resid == resid
                bA5 = u.atoms.name.isin(['C', 'O'])
                #u.atoms.loc[bA1 & bA4 & bA5, 'bfactor'] = 0.0
                u.atoms.loc[bA1 & bA4, 'bfactor'] = 0.0
                saved.append(record)
                
            elif len(u.atoms[bA1 & bA2]) > 0 and len(u.atoms[bA1 & bA3]) == 0:
                bA4 = u.atoms.resid == resid
                bA5 = u.atoms.name.isin(['N', 'H', 'HN'])
                #u.atoms.loc[bA1 & bA4 & bA5, 'bfactor'] = 0.0
                u.atoms.loc[bA1 & bA4, 'bfactor'] = 0.0
                saved.append(record)
                
            elif len(u.atoms[bA1 & bA2]) == 0 and len(u.atoms[bA1 & bA3]) == 0:
                bA4 = u.atoms.resid == resid
                bA5 = u.atoms.name.isin(['C', 'O', 'N', 'H', 'HN'])
                #u.atoms.loc[bA1 & bA4 & bA5, 'bfactor'] = 0.0
                u.atoms.loc[bA1 & bA4, 'bfactor'] = 0.0
                saved.append(record)

        df = u.atoms[u.atoms.bfactor > 0.5]
        new = Universe(data=df)
        new.dimensions = u.dimensions
        new.cell = new.cell
        new.write(output)


    def vis(self):
        with open(self.workdir + '/vis.cxc', 'w') as W:
            W.write('''
open step1_input.pdb
open step9_final.pdb

camera ortho
set bgColor white
lighting shadow true
lighting soft
graphics silhouettes true
cartoon suppress False

### Loop residues
select (#2 & @@bfactor < 0.5 & @CA) residues true

# cartoon/ribbon
color sel black target cr

# atom
color sel &  C black target a
show  sel target a

~select
''')
        with open(self.workdir + '/vis.tcl', 'w') as W:
            W.write('''
# vmd -e vis.tcl
display projection Orthographic
display shadows on
display ambientocclusion on
# color Display Background white
axes location Off

mol new step1_input.pdb
mol modstyle 0 0 NewCartoon 0.300000 30.000000 4.100000 0
mol modcolor 0 0 ColorID 1 # Red
mol modmaterial 0 0 AOShiny

mol new step9_final.pdb
mol modstyle 0 1 NewCartoon 0.300000 30.000000 4.100000 0
mol modmaterial 0 1 AOShiny

mol addrep 1
mol modstyle 1 1 Licorice 0.300000 30.000000 30.000000
mol modmaterial 1 1 AOShiny
mol modselect 1 1 "same residue as (beta < 0.5 and name CA)"
''')


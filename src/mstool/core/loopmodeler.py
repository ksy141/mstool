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
from   .dmsfile               import DMSFile
from   .rem                   import REM
from   .checkstructure        import CheckStructure
from   .ungroup               import Ungroup
from   .solvate_martini       import SolvateMartini

from   ..utils.protein_sel    import three2one
from   ..utils.dump           import dumpsql
from   ..utils.openmmutils    import getEnergy, runMartiniEM, runMartiniNPT, runMartiniEMNPT


from   openmm.app             import *
from   openmm                 import *
from   openmm.unit            import *

class LoopModeler:
    def __init__(self, protein, fasta, ref=None, workdir='workdir', 
        A=100, C=50, soft=True, mutate=True, t=15.0, extend_termini={},
        mapping=[], mapping_add=[], ff=[], ff_add=[],
        fc1 = 50, fc2 = 2000, Kchiral=300.0, Kpeptide=300.0, 
        nsteps=0, dt=0.002, dcdfreq=1000, csvfreq=1000, helix=False,
        addnbtype=['ZCARBON', 0.34, 1.51]):

        self.protein     = protein
        self.fasta       = fasta
        self.workdir     = workdir
        self.A           = A
        self.C           = C
        self.soft        = soft
        self.mutate      = mutate
        self.t           = t
        self.mapping     = mapping
        self.mapping_add = mapping_add
        self.fc1         = fc1
        self.helix       = helix
        self.addnbtype   = addnbtype

        ### step1: workdir
        self.create_workdir(protein, 
                            workdir + '/step1_input.pdb')

        ### step1.1: translate
        self.translate(workdir + '/step1_input.pdb',
                       workdir + '/step1_translated.pdb')


        ### step2: mutate
        self.fix_mutation(workdir + '/step1_translated.pdb',
                          workdir + '/step2_mutate.pdb')


        ### step2.1: remove dangling CONH
        self.remove_dangling_CONH(workdir + '/step2_mutate.pdb',
                                  workdir + '/step2_mutate_new.pdb')


        ### step3: map
        Map(structure   = workdir + '/step2_mutate.pdb',
            out         = workdir + '/step3_cg.pdb',
            mapping     = mapping,
            mapping_add = mapping_add,
            BB2CA       = True,
            add_notavailableAAAtoms = True)


        ### step4: Fill martini loops
        Fill(structure      = workdir + '/step3_cg.pdb',
             out            = workdir + '/step4_cgfilled.pdb',
             sequence       = fasta,
             extend_termini = extend_termini,
             mapping        = mapping,
             mapping_add    = mapping_add)


        ### step4.2: Solvate
        if t == 0 or t is None:
            shutil.copy(workdir + '/step4_cgfilled.pdb',
                        workdir + '/step4_solvated.pdb')
            nonbondedMethod = 'CutoffNonPeriodic'

        else:
            SolvateMartini(structure = workdir + '/step4_cgfilled.pdb',
                           out       = workdir + '/step4_solvated.pdb',
                           removedr  = 4.0,
                           center    = False,
                           conc      = 0.0,
                           pbc       = True,
                           t         = t)
            nonbondedMethod = 'CutoffPeriodic'


        ### step5: create a martini FF dms
        self.makemartini(workdir + '/step4_solvated.pdb',
                         workdir + '/step5_ff.dms')


        ### step6: run Martini simulation
        dms = DMSFile(workdir + '/step5_ff.dms')
        dms.createSystem(REM=True, tapering='shift', martini=True, nonbondedMethod=nonbondedMethod)
        dms.runEMNPT(dt=dt, nsteps=nsteps, dcdfreq=dcdfreq, csvfreq=csvfreq, semiisotropic=False, out=workdir + '/step6_minimized.dms')
        Universe(workdir + '/step6_minimized.dms').write(workdir + '/step6_minimized.pdb')

        #runMartiniEM(dms_in = workdir + '/step5_ff.dms',
        #             out    = workdir + '/step6_minimized.pdb',
        #             soft   = soft,
        #             nonbondedMethod = nonbondedMethod,
        #             nonbondedCutoff = 1.1)

        if t == 0 or t is None or nsteps == 0:
            shutil.copy(workdir + '/step6_minimized.pdb',
                        workdir + '/step6_NPT.pdb')
        #else:
        #    runMartiniEMNPT(dms_in = workdir + '/step5_ff.dms',
        #                    pos_in = workdir + '/step6_minimized.pdb',
        #                    out    = workdir + '/step6_NPT.pdb',
        #                    soft   = False,
        #                    nonbondedMethod = nonbondedMethod,
        #                    nonbondedCutoff = 1.1,
        #                    dt      = dt,
        #                    nsteps  = nsteps,
        #                    dcdfreq = dcdfreq,
        #                    csvfreq = csvfreq)


        ### step7: ungroup (output must be a pdb so that openMM recognizes protein residues)
        Ungroup(structure   = workdir + '/step6_NPT.pdb',
                out         = workdir + '/step7_backmapped.pdb',
                mapping     = mapping,
                mapping_add = mapping_add, 
                ff          = ff,
                ff_add      = ff_add,
                backbone    = True,
                water_resname = 'DONOTBACKMAPWATER')


        ### step8: REM + steer MD
        ### Kpeptide is reduced because a crystal structure can have cis peptide
        REM(structure   = workdir + '/step7_backmapped.pdb',
            out         = workdir + '/step8_minimized.pdb',
            refposre    = workdir + '/step2_mutate_new.pdb',
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
            Kchiral     = Kchiral * 0.25,
            Kpeptide    = Kpeptide * 0.25,
            Kcistrans   = 0.0,
            turn_off_EMNVT= True)


        ### check structure
        CheckStructure(structure   = workdir + '/step8_minimized.pdb',
                       mapping     = mapping,
                       mapping_add = mapping_add)


        ### translate back
        u = Universe(workdir + '/step8_minimized.pdb')
        u.atoms[['x','y','z']] -= self.dr
        u.write(workdir + '/step8_minimized.pdb')


        ### step9: combine
        u = self.combine_two(xtal    = workdir + '/step1_input.pdb', 
                             backmap = workdir + '/step8_minimized.pdb')
        u.write(workdir + '/step9_final.pdb')

        ### vis
        self.vis()



    def create_workdir(self, structure, out):
        if os.path.exists(self.workdir):
            #shutil.rmtree(self.workdir)
            raise Exception(self.workdir + ' already exists')
        os.mkdir(self.workdir)

        ### Read and save a protein structure
        u = Universe(structure)
        u.atoms.bfactor = 1.0
        u.atoms.segname = ''
        u.write(out)


    def translate(self, structure, out):
        u = Universe(structure)
        cog = u.atoms[['x','y','z']].mean(axis=0)
        u.atoms[['x','y','z']] -= cog
        dx = max(u.atoms['x']) - min(u.atoms['x'])
        dy = max(u.atoms['y']) - min(u.atoms['y'])
        dz = max(u.atoms['z']) - min(u.atoms['z'])
        pbc = max(dx, dy, dz) + self.t * 2
 
        u.dimensions = [pbc, pbc, pbc, 90, 90, 90]
        u.cell = [[pbc, 0, 0], [0, pbc, 0], [0, 0, pbc]]

        u.atoms[['x','y','z']] += np.array([pbc/2, pbc/2, pbc/2])
        self.dr = -cog + np.array([pbc/2, pbc/2, pbc/2])

        ### HSE to HIE
        ### This is the only histidine residue that openMM cannot recognize...
        ### HIS is automatically changed to HSD
        bA = u.atoms.resname == 'HSE'
        u.atoms.loc[bA, 'resname'] = 'HIE'

        bA2 = u.atoms.resname == 'ACE'
        bA3 = u.atoms.resname == 'NMA'
        bA4 = u.atoms.resname == 'NME'

        u = Universe(data = u.atoms[~(bA2 | bA3 | bA4)])
        u.dimensions = [pbc, pbc, pbc, 90, 90, 90]
        u.cell = [[pbc, 0, 0], [0, pbc, 0], [0, 0, pbc]]
        u.write(out)


    def fix_mutation(self, structure, out):
        if self.mutate:
            Mutate(structure = structure,
                   fasta     = self.fasta,
                   out       = out)
        else:
            shutil.copy(structure, out)
            #u = Universe(structure).write(out)


    def makemartini(self, structure, out):
        martini = ReadMartini(addnbtype=self.addnbtype)
        # dms default unit: kcal/mol/A^2
        # 50 kcal/mol/A^2 -> 0.5 * 50 * 4.184 * 100 kJ/mol/nm^2
        MartinizeDMS(dms_in  = structure,
                     martini = martini, 
                     out     = out,
                     fcx     = self.fc1,
                     fcy     = self.fc1,
                     fcz     = self.fc1,
                     helix   = self.helix,
                     addnbtype = self.addnbtype[0])
        dumpsql(out)


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
            bA4 = xtal.atoms.resname.isin(['ACE', 'NMA', 'NME'])

            xtalatoms = xtal.atoms[bA1 & bA2 & bA3 & ~bA4]
            if len(xtalatoms) == 0:
                bfactors.append(0.0)
                #atom.bfactor = 0.0

            elif len(xtalatoms) == 1:
                bfactors.append(1.0)
                #atom.bfactor = 1.0

            else:
                assert 0 == 1, print(xtalatoms)

        backmap.atoms.bfactor = bfactors
        return backmap


    def remove_dangling_CONH(self, structure, out):
        u = Universe(structure)

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
        new.write(out)


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


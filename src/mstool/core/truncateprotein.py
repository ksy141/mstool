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
from   .ungroup               import Ungroup
from   .solvate_martini       import SolvateMartini

from   ..utils.protein_sel    import three2one
from   ..utils.dump           import dumpsql
from   ..utils.openmmutils    import getEnergy, runMartiniEM, runMartiniNPT, runMartiniEMNPT
from   ..utils.amberselection import amberSelection
from   ..utils.add_termini_atoms import addTerminiAtoms

from   ..utils.openmmutils import runEM, addPosre

class TruncateProtein:
    def __init__(self, protein, truncate, workdir='workdir', 
        A=100, C=50, soft=True, mutate=True, extend_termini={},
        mapping=[], mapping_add=[], ff=[], ff_add=[],
        fc1 = 50, fc2 =5000, Kchiral=300.0, Kpeptide=300.0, 
        add_termini_atoms=True):

        self.protein     = protein
        self.truncate    = truncate
        self.workdir     = workdir
        self.A           = A
        self.C           = C
        self.soft        = soft
        self.mutate      = mutate
        self.mapping     = mapping
        self.mapping_add = mapping_add
        self.fc1         = fc1
        self.fc2         = fc2
        self.add_termini_atoms = add_termini_atoms

        ### step1: workdir
        self.create_workdir(protein, 
                            workdir + '/step1_input.pdb')

        ### step1: add termini atoms
        if add_termini_atoms:
            addTerminiAtoms(workdir + '/step1_input.pdb',
                            workdir + '/step1_input_termini.pdb')
        else:
            shutil.copy(workdir + '/step1_input.pdb',
                        workdir + '/step1_input_termini.pdb')

        ### step1: translate
        self.translate(workdir + '/step1_input_termini.pdb',
                       workdir + '/step1_translated.pdb')


        ### step1: read truncate
        self.read_truncate()


        ### step1: truncate
        self.truncate_protein(workdir + '/step1_translated.pdb',
                              workdir + '/step1_truncate.pdb')

        ### step2: reassign resid
        self.reassignresid(workdir + '/step1_truncate.pdb',
                           workdir + '/step2_new.pdb', 
                           workdir + '/seq.fasta')

        ### step3: map
        Map(structure   = workdir + '/step2_new.pdb',
            out         = workdir + '/step3_cg.pdb',
            mapping     = mapping,
            mapping_add = mapping_add,
            BB2CA       = True)


        ### step4: Fill martini loops
        Fill(structure      = workdir + '/step3_cg.pdb',
             out            = workdir + '/step4_cgfilled.pdb',
             sequence       = workdir + '/seq.fasta',
             extend_termini = extend_termini,
             mapping        = mapping,
             mapping_add    = mapping_add)


        ### step5: create a martini FF dms
        self.makemartini(workdir + '/step4_cgfilled.pdb',
                         workdir + '/step5_ff.dms')


        ### step6: run Martini simulation
        runMartiniEM(dms_in = workdir + '/step5_ff.dms',
                     out    = workdir + '/step6_minimized.pdb',
                     soft   = soft,
                     nonbondedMethod = 'CutoffNonPeriodic',
                     nonbondedCutoff = 1.1)

        ### step7: ungroup (output must be a pdb so that openMM recognizes protein residues)
        Ungroup(structure   = workdir + '/step6_minimized.pdb',
                out         = workdir + '/step7_backmapped.pdb',
                mapping     = mapping,
                mapping_add = mapping_add, 
                ff          = ff,
                ff_add      = ff_add,
                backbone    = True)


        ### step8: REM + steer MD
        ### Kpeptide is reduced because a crystal structure can have cis peptide
        REM(structure   = workdir + '/step7_backmapped.pdb',
            out         = workdir + '/step8_minimized.pdb',
            refposre    = workdir + '/step2_new.pdb',
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
            turn_off_EMNVT = True)


        ### check structure
        CheckStructure(structure   = workdir + '/step8_minimized.pdb',
                       mapping     = mapping,
                       mapping_add = mapping_add)


        ### step9: combine
        u = self.combine_two(xtal    = workdir + '/step1_input_termini.pdb', 
                             backmap = workdir + '/step8_minimized.pdb',
                             out     = workdir + '/step9_final.pdb')

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
        u   = Universe(structure)
        cog = u.atoms[['x','y','z']].mean(axis=0)
        dx  = max(u.atoms['x']) - min(u.atoms['x'])
        dy  = max(u.atoms['y']) - min(u.atoms['y'])
        dz  = max(u.atoms['z']) - min(u.atoms['z'])
        pbc = max(dx, dy, dz) + 15.0 * 2
 
        u.dimensions = [pbc, pbc, pbc, 90, 90, 90]
        u.cell = [[pbc, 0, 0], [0, pbc, 0], [0, 0, pbc]]
        self.pbc = pbc

        ### HSE to HIE
        ### This is the only histidine residue that openMM cannot recognize...
        ### HIS is automatically changed to HSD
        bA = u.atoms.resname == 'HSE'
        u.atoms.loc[bA, 'resname'] = 'HIE'
        u.write(out)


    def read_truncate(self):
        self.truncate_data = {}
        with open(self.truncate) as fp:
            for line in fp:
                line = line.strip()
                if not line: continue
                one, two, three = line.split()
                start = amberSelection(one)
                end   = amberSelection(three)
                chain = start['chain']
                assert start['chain'] == end['chain'],\
                    'chain is different %s %s' %(start['chain'], end['chain'])

                if chain not in self.truncate_data.keys():
                    self.truncate_data[chain] = []
                self.truncate_data[chain].append([start['resid'], two, end['resid']])

        for key in self.truncate_data.keys():
            self.truncate_data[key].sort()

        print(self.truncate_data)


    def truncate_protein(self, structure, out):
        u  = Universe(structure)
        bA = u.atoms.resname.isin(three2one.keys())
        proteinu = Universe(data=u.atoms[bA])
        proteinu.dimensions = u.dimensions
        proteinu.cell       = u.cell

        removes = []
        for key, value in self.truncate_data.items():
            bAchain  = proteinu.atoms.chain == key

            for v in value:
                bAresid = proteinu.atoms.resid.isin(np.arange(v[0]+1, v[2]))
                removes.extend(list(proteinu.atoms[bAchain & bAresid].index))

        proteinu.atoms.drop(removes, inplace=True)
        proteinu.write(out)
        return out


    def reassignresid(self, structure, out, fasta):
        u = Universe(structure)
        chains = list(dict.fromkeys(u.atoms.chain))
        self.truncate_xtal = {}
        self.truncate_new  = {}

        seq = ''
        for c in chains:
            self.truncate_xtal[c] = {}
            self.truncate_new[c]  = {}

            seq  += f'> {c}\n'
            bAchain = u.atoms.chain == c

            minresid = min(u.atoms[bAchain].resid)
            maxresid = max(u.atoms[bAchain].resid)

            newresid = 1
            for r in np.arange(minresid, maxresid+1):
                bAresid = u.atoms.resid == r
                residue = u.atoms[bAchain & bAresid]

                # residue that doesn't exist in truncate
                if len(residue) == 0: continue

                resname = residue.resname.values[0]
                resid   = residue.resid.values[0]
                u.atoms.loc[bAchain & bAresid, 'resid'] = newresid
                seq    += three2one[resname]

                # chain that does not have truncations
                if c not in self.truncate_data.keys():
                    self.truncate_xtal[c][resid] = newresid

                # chain that has truncations
                if c in self.truncate_data.keys():
                    added = False

                    for v in self.truncate_data[c]:
                        if r == v[0]:
                            # flexible end
                            self.truncate_new[c][resid] = newresid

                            for ii, addseq in enumerate(v[1], 1):
                                self.truncate_new[c][resid + ii] = newresid + ii
                            newresid += len(v[1])
                            seq += v[1]
                            added = True

                        elif r == v[2]:
                            # flexible end
                            self.truncate_new[c][resid] = newresid
                            added = True

                    if not added:
                        self.truncate_xtal[c][resid] = newresid

                newresid += 1

            seq += '\n\n'

        u.write(out)
        with open(fasta, 'w') as fp:
            fp.write(seq)
            

    def makemartini(self, structure, out):
        martini = ReadMartini()
        # dms default unit: kcal/mol/A^2
        # 50 kcal/mol/A^2 -> 0.5 * 50 * 4.184 * 100 kJ/mol/nm^2
        MartinizeDMS(dms_in  = structure,
                     martini = martini, 
                     out     = out,
                     fcx     = self.fc1,
                     fcy     = self.fc1,
                     fcz     = self.fc1)
        dumpsql(out)


    def combine_two(self, xtal, backmap, out):
        xtal     = Universe(xtal)
        xtal.atoms.bfactor = 1.0

        backmap  = Universe(backmap)
        backmap.atoms.bfactor = 0.0

        bA       = xtal.atoms.resname.isin(three2one.keys())
        proteinu = Universe(data=xtal.atoms[bA])
        otheru   = Universe(data=xtal.atoms[~bA])

        chains   = list(dict.fromkeys(proteinu.atoms.chain)) 

        saves = []
        for c in chains:
            bA1chain = proteinu.atoms.chain == c
            bA1resid = proteinu.atoms.resid.isin(list(self.truncate_xtal[c].keys()))
            saves.append(proteinu.atoms[bA1chain & bA1resid])

            bA2chain = backmap.atoms.chain == c
            bA2resid = backmap.atoms.resid.isin(list(self.truncate_new[c].values()))
            atoms    = pd.concat([backmap.atoms[bA2chain & bA2resid]], ignore_index=True)

            for key, value in self.truncate_new[c].items():
                bA3 = atoms.resid == value
                atoms.loc[bA3, 'resid'] = key
            saves.append(atoms)
        

        ### optional EM - comment out this session
        u = Universe(data=pd.concat(saves, ignore_index=True))
        u.atoms.sort_values(by=['chain', 'resid', 'name'], inplace=True)
        u.dimensions = [self.pbc, self.pbc, self.pbc, 90, 90, 90]
        u.cell       = [[self.pbc, 0, 0], [0, self.pbc, 0], [0, 0, self.pbc]]
        u.write(self.workdir + '/step9_protein.pdb')

        cef = addPosre(Universe(self.workdir + '/step9_protein.pdb'), 
                       bfactor_posre=0.5, fcx=500, fcy=500, fcz=500)

        runEM(structure  = self.workdir + '/step9_protein.pdb', 
             forcefield = 'charmm36.xml',
             out = self.workdir + '/step9_protein_EM.pdb',
             addForces = [cef])

        u = Universe(self.workdir + '/step9_protein_EM.pdb')
        newu = Universe(data=pd.concat([u.atoms, otheru.atoms], ignore_index=True))
        ###
        

        newu = Universe(data=pd.concat([*saves, otheru.atoms], ignore_index=True))
        newu.dimensions = xtal.dimensions
        newu.cell       = xtal.cell
        newu.atoms.sort_values(by=['chain', 'resid', 'name'], inplace=True)
        newu.write(out)


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


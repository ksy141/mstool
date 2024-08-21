import os
import shutil
import numpy as np
from  openmm.app import *
from .readmappings     import ReadMappings
from .readxml          import ReadXML
from .ungroup          import Ungroup
from .rem              import REM
from .universe         import Universe, Merge
from .checkstructure   import CheckStructure
from .checktetrahedron import CheckTetrahedron

class Backmap:

    def __init__(self, structure, workdir='workdir', AA=None, protein=None, backbone=True,
                 mapping = [], mapping_add = [],
                 ff      = [], ff_add = [],
                 Kchiral=300, Kpeptide=300, Kcistrans=300, Kdihedral=300,
                 fcx=1000.0, fcy=1000.0, fcz=1000.0,
                 rock=None, rockCtype='CTL3', rockHtype='HAL3',
                 rcut=1.2, pbc=True, A=100, C=50,
                 add_bonds = True, remversion='v4',
                 water_resname='W', water_chain=None, water_number=4, water_fibor=2.0, water_chain_dms=True, 
                 use_AA_structure=True, AA_structure=[], AA_structure_add=[], AA_shrink_factor=0.8,
                 use_existing_workdir=False, fileindex=1, pdbsave=True, cospower=2,
                 nsteps=10000, rem_nsteps=0, turn_off_EMNVT=False, T=310, sanitizeMartini=True,
                 changename={':CHOL':':CHL1',':ION@NA':':SOD@SOD',':NA@NA':':SOD@SOD',
                             ':ION@CL':':CLA@CLA',':CL@CL':':CLA@CLA',':ION@CA':':CAL@CAL'},
                 changename_add={}):

        ### workdir
        if not use_existing_workdir: os.mkdir(workdir)
        if use_existing_workdir and not os.path.exists(workdir): os.mkdir(workdir)

        ### save args
        args = locals()
        with open(workdir + '/args_backmap.txt', 'w') as W:
            for key, value in args.items():
                if key == 'self': continue
                W.write(f'{key:30} = {value}\n')

        ### Read Mapping and XML and compare these two
        if not isinstance(ff_add, list): ff_add = [ff_add]
        if not isinstance(mapping_add, list): mapping_add = [mapping_add]
        map = ReadMappings(mapping=mapping, mapping_add=mapping_add)
        xml = ReadXML(ff=ff, ff_add=ff_add)
        self.checkMappingXML(map, xml)

        ### Read AA and determine whether one should use protein or rock
        if AA:
            try:
                pdb = PDBFile(AA)
                forcefield = ForceField('charmm36.xml', 'charmm36/water.xml')
                system = forcefield.createSystem(pdb.topology)
                protein = AA
                rock    = None
                print('All of the residues in the AA argument is recognizable in openMM')
                print(f'Using protein="{AA}"')

            except:
                rock = AA
                protein = None
                print('Some of the residues in the AA argument is not recognizable in openMM')
                print(f'Using rock="{AA}"')

        ### Ungroup
        if sanitizeMartini:
            structure = Universe(structure)
            for key, value in changename_add.items():
                changename[key] = value
            structure.changeName(changename)

        Ungroup(structure, out=workdir + f'/step{fileindex}_ungroup.dms', 
                mapping=mapping, mapping_add=mapping_add, backbone=backbone,
                water_resname=water_resname,
                water_chain=water_chain, water_number=water_number,
                water_fibor=water_fibor, water_chain_dms=water_chain_dms,
                sort=True, use_AA_structure=use_AA_structure,
                AA_structure=AA_structure, AA_structure_add=AA_structure_add, AA_shrink_factor=0.8)
        
        REM(structure   = workdir + f'/step{fileindex}_ungroup.dms', 
            outrem      = workdir + f'/step{fileindex+1}_rem.dms',
            out         = workdir + f'/step{fileindex+2}_em.dms',
            rockout     = workdir + f'/step{fileindex+2}_rock.dms',
            nonrockout  = workdir + f'/step{fileindex+2}_nonrock.dms',
            rockprefix  = workdir + '/ROCK',
            protein     = protein,
            rock        = rock,
            rockCtype   = rockCtype,
            rockHtype   = rockHtype,
            mapping     = mapping, 
            mapping_add = mapping_add,
            ff          = ff, 
            ff_add      = ff_add,
            A           = A,
            C           = C,
            cospower    = cospower,
            rcut        = rcut,
            pbc         = pbc,
            nsteps      = nsteps,
            rem_nsteps  = rem_nsteps,
            T           = T,
            version     = remversion,
            Kchiral     = Kchiral,
            Kpeptide    = Kpeptide, 
            Kcistrans   = Kcistrans, 
            Kdihedral   = Kdihedral,
            turn_off_EMNVT = turn_off_EMNVT)
        

        ### CheckStructure
        CheckStructure(structure   = workdir + f'/step{fileindex+2}_em.dms', 
                       log         = workdir + '/log.txt',
                       mapping     = mapping,
                       mapping_add = mapping_add)


        ### Combine
        if rock:
            if os.path.exists('ROCK.dms'): os.rename('ROCK.dms', workdir + '/ROCK.dms')
            if os.path.exists('ROCK.xml'): os.rename('ROCK.xml', workdir + '/ROCK.xml')

            u1 = Universe(AA)
            #memory-intensive job
            #if len(u1.bonds) == 0: u1.addBondFromDistance()

            shutil.copyfile(workdir + f'/step{fileindex+2}_nonrock.dms', 
                            workdir + f'/step{fileindex+3}_nonprotein.dms')

            u2 = Universe(workdir + f'/step{fileindex+2}_nonrock.dms')
            u = Merge(u1.atoms, u2.atoms)
            u.bonds = list(u1.bonds) + list(len(u1.atoms) + np.array(u2.bonds))
            u.dimensions = u2.dimensions
            u.cell       = u2.cell
            u.write(workdir + f'/step{fileindex+3}_final.dms')
            u.write(workdir + f'/step{fileindex+3}_final.pdb')
            
            u1.dimensions = u2.dimensions
            u1.cell       = u2.cell
            u1.write(workdir + f'/step{fileindex+3}_protein.pdb')

        else:
            shutil.copyfile(workdir + f'/step{fileindex+2}_em.dms',
                            workdir + f'/step{fileindex+3}_final.dms')
            Universe(workdir + f'/step{fileindex+3}_final.dms').write(
                     workdir + f'/step{fileindex+3}_final.pdb')
        
        ### PDB
        if pdbsave:
            step1file = workdir + f'/step{fileindex}_ungroup'
            Universe(step1file + '.dms').write(step1file + '.pdb')

            step2file = workdir + f'/step{fileindex+1}_rem'
            Universe(step2file + '.dms').write(step2file + '.pdb')

            step3file = workdir + f'/step{fileindex+2}_em'
            Universe(step3file + '.dms').write(step3file + '.pdb')

            #step3file = workdir + f'/step{fileindex+3}_final'
            #Universe(step3file + '.dms').write(step3file + '.pdb')
            #CheckTetrahedron(step3file + '.dms', ff=ff, ff_add=ff_add)


    def checkMappingXML(self, map, xml):
        skips = ['ASH', 'GLH', 'HIS', 'HID', 'HIE', 'HIP', 'CHOL', 'ION', 'NA', 'CL']
        for resname in map.RESI.keys():
            if resname in skips: continue

            if resname not in xml.RESI.keys():
                print(f'WARNING: {resname} is defined in mapping files but not in forcefield files. This will likely cause a probelm if your system contains this residue.')
                continue
            
            aa_map = map.RESI[resname]['AAAtoms']
            aa_xml = xml.RESI[resname]['names']
            if set(aa_map) == set(aa_xml) and len(aa_map) == len(aa_xml): continue

            for atom in set(aa_map) - set(aa_xml):
                print(f'WARNING: {resname} has {atom} in mapping files but not in forcefield files. This will likely cause a probelm if your system contains this atom.')
            
            for atom in set(aa_xml) - set(aa_map):
                print(f'WARNING: {resname} has {atom} in forcefield files but not in mapping files. This will likely cause a probelm if your system contains this atom.')


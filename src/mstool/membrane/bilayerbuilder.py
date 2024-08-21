import os
import shutil
import numpy as np
from  .bilayer               import Bilayer, make_rect2
from  ..core.map             import Map
from  ..core.readmartini     import ReadMartini
from  ..core.solvate_martini import SolvateMartini, ionize
from  ..core.martinizedms    import MartinizeDMS
from  ..core.dmsfile         import DMSFile
from  ..core.universe        import Universe
from  ..core.backmap         import Backmap
from  ..utils.dump           import dumpsql
from  ..utils.openmmutils    import addPosre, addPosrePeriodicZ, addFlatBottomZ
from  ..utils.protein_sel    import one2three, three2one
pwd = os.path.dirname(os.path.realpath(__file__))

# martini3 ALA SC1: addnbtype=['ZCARBON', 0.34, 1.51]
# CHARMM36 HA3: addnbtype=['ZCARBON', 0.23876, 0.1]

class BilayerBuilder:
    def __init__(self, workdir='workdir', protein=None, upper={}, lower={}, between={}, sep=0.0, dx=8.0, area=None, water=50.0, 
                 rcut=4.0, rcut_use_enclosed_protein=False, mode='shift', 
                 dN=10, rockCtype='CTL3', rockHtype='HAL3',
                 martini=[], martini_add=[], lipidpath=pwd+'/../../../FF/martini2.2/structures/',
                 mapping=[], mapping_add=[],
                 ff=[], ff_add=[],
                 removedr=6.0, aa_nsteps=0, fc=10.0, fc2=10.0,
                 dt_rem=0.025, cg_nsteps_rem=100000,
                 dt=0.020, cg_nsteps=500000,
                 frictionCoeff=10.0, barfreq=10, nonbondedCutoff=1.1, 
                 dcdfreq=1000, csvfreq=1000, tension=0.0, P=1.0, T=310, 
                 improper_prefactor=0.99, use_existing_folder=False,
                 hydrophobic_thickness=30, ionconc=0.15,
                 use_AA_structure=True, AA_shrink_factor=0.8,
                 remove_solvent=False,
                 solvate=True,
                 tapering='shift',
                 protein_dict_pbc_shrink_factor=0.9,
                 changename={}, addnbtype=['ZCARBON', 0.34, 1.51]):

        '''Bilayer builder.
        Parameters
        ----------
        workdir: str=workdir
            name for the working directory (All the results will be saved here).
            If protein is not provided, the final structure will be step6_final.dms, inside the working directory.
            If protein is provided, the protein structure will be protein.dms; lipids will be step6_nonprotein.dms; combined will be step6_final.dms.
        protein: str=None
            provide a protein input structure file (.pdb or .dms). None (default) for protein is not needed in your membrane.
            Protein should be prealigned and oriented before providing a structure to the workflow 
            with respect to a membrane whose center is (0,0,0) and whose normal is [0, 0, 1].
        upper: dict={}
            type (key) and number (value) of each lipid in the upper leaflet.
            e.g., upper={'POPC': 100, 'DOPC': 100} will place 100 POPC molecules and 100 DOPC molecules in the upper leaflet.
        lower: dict={}
            type and number of each lipid in the lower leaflet.
            e.g., lower={'POPC': 100, 'DOPC': 100} will place 100 POPC molecules and 100 DOPC molecules in the lower leaflet.
        dx: float=8.0
            distance between the lipids at the initial structure.
        water: float=50.0
            water thickness in Z dimension (A).
        rcut: float=3.0
            lipids that are closer than rcut (A) to protein will be shifted (mode=shift) or removed (mode=remove)
        mode: str=shift
            shift will move lipids that overlap with protein. (Therefore, the total number of lipids will be as given.)
            remove will remove lipids that overlap with protein. (Therefore, the total number of lipids will be reduced than given.)
        dN: int=5
            Number of additional layers in XY dimensions that will be made if you shift overlapped lipids.
            dN = 5 (default) is usually fine, but if you have a protein-crowded membrane structure,
            you should increase this value.
        martini: list=[]
            Martini force field. If not provided, $mstool/FF/martini2.2 will be read.
        martini_add: list=[]
            Additional Martini force field. Useful when you add new molecules.
        lipidpath: $mstool/FF/martini2.2/structures/
            Path to a folder that contains the structures of lipids.
            Phospholipids that have names of GL*/C*/D* can be internally constructed.
            However, cholesterols and other molecules are not internally constructed.
            Therefore, users should put lipid structures that cannot be internally constructed to this path.
            The default is $mstool/FF/martini2.2/structures/
            The filename should be ``name_of_moleculetype.pdb`` or ``name_of_moleculetype.dms``.
            The lipid should have an orientation of +z, and approximately centered at the origin.
            Imagine that this lipid is located at the upper leaflet of a bilayer whose center is 0.
        mapping: list=[]
            Mapping files. If not provided, $mstool/mapping/martini.lipid.c36.dat and $mstool/mapping/martini.protein.c36m.dat will be read.
        mapping_add: list=[]
            Additional mapping files. Useful when you add new molecules.
        ff: list=[]
            openMM Charmm36 XML files. If not provided, $mstool/FF/charmm36/* will be read.
        ff_add: list = []
            Additional openMM Charmm36 XML files. Useful when you add new molecules.
        removedr: float=8.5
            Remove water beads that are closer than removedr to solutes (protein + lipids).
        dt: float=0.03
            Integration time (ps) for coarse-grained simulations.
            If you encounter numerical issue (Energy is Nan, or particles are all over the place), reduce dt while increasing nsteps.
        cg_nsteps: int=100000
            Number of coarse-grained simulation steps.
            If your membrane does not look equilibrated (e.g., water not equally dispersed in solution), increase this number.
        aa_nsteps: int=10000
            Number of all-atom simulation steps after backmapping
        frictionCoeff: float=5.0
            Friction coefficient for coarse-grained simulations.
        barfreq: int=1
            Frequency of MonteCarloMembraneBarostat.
        nonbondedCutoff: float=1.1
            Cutoff distance in coarse-grained simulations.
        improper_prefactor: float=0.99
            k*(theta-theta0)^2 causes discontinuity at +-180 degree.
            Therefore, in mstool, it is replaced with k*(acos(cos(theta-theta0)))^2.
            However, this causes numerical stability issue in Linux (no problem in Mac).
            Therefore, I have to introduce a prefactor, so that cosine values do not exceed 1.
            k*(acos(improper_prefactor * cos(theta-theta0)))^2.
            It works fine for most of the case. However, when you compare potential energy calculated with gromacs and openMM, 
            make sure you change it to 1.0.
        hydrophobic_thickness: float=30
            Hydrophobic thickness in A.
        ionconc: float=0.15
            Ion concentration in M.
        '''

        
        ### workdir
        if (not use_existing_folder) or (use_existing_folder and not os.path.exists(workdir)):
            os.mkdir(workdir)

        ### save args
        args = locals()
        with open(workdir + '/args_builder.txt', 'w') as W:
            for key, value in args.items():
                if key == 'self': continue
                W.write(f'{key:30} = {value}\n')

        ### area
        if area:
            ntotal  = int(area / dx**2)
            nupper0 = sum(upper.values())
            nlower0 = sum(lower.values())
            for key, value in upper.items():
                upper[key] = int(ntotal * value / nupper0)
            for key, value in lower.items():
                lower[key] = int(ntotal * value / nlower0)
            print(upper)
            print(lower)

            assert sum(upper.values()) > 0, 'the number of upper lipids = 0; increase area'
            assert sum(lower.values()) > 0, 'the number of upper lipids = 0; increase area'

        ### protein
        if protein:
            if isinstance(protein, str):
                proteinU = Universe(protein)
 
            elif isinstance(protein, Universe):
                proteinU = protein
 
            elif isinstance(protein, dict):
                # protein = {'GPCR.pdb': 3, 'ompf.pdb': 2}
                # proteinlayer_keys = ['GPCR.pdb', 'ompf.pdb']
                # proteinlayer_list = [0, 0, 0, 1, 1]
                proteinlayer_key  = [key for key in protein.keys()]
                proteinlayer_list = []
                for idx, number in enumerate(protein.values()):
                    proteinlayer_list += [idx] * number
                np.random.shuffle(proteinlayer_list)
                proteindata = {'name': [], 'resname': [], 'resid': [], 'chain': [], 'x': [], 'y': [], 'z': []}
                
                proteinN = len(proteinlayer_list)
                pbcx = np.sqrt(max([sum(upper.values()), sum(lower.values())])) * dx * protein_dict_pbc_shrink_factor
                proteinP, _ = make_rect2(proteinN, dx=pbcx/int(np.sqrt(proteinN)), dN=1)
                for i in range(len(proteinP)):
                    key = proteinlayer_key[proteinlayer_list[i]]
                    if isinstance(key, str):
                        # structure information exists.
                        univ         = Universe(key)
                        positions    = univ.atoms[['x','y','z']].to_numpy()
                        save_resname = univ.atoms['resname'].tolist()
                        save_name    = univ.atoms['name'].tolist()
                        save_resid   = univ.atoms['resid'].tolist()
                        save_chain   = np.char.add(univ.atoms['chain'].tolist(), [str(i)] * len(positions))
        
                    elif isinstance(key, Universe):
                        # a structure e.g., protein or lipid with a structure
                        positions     = key.atoms[['x','y','z']].to_numpy()
                        save_resname  = key.atoms['resname'].tolist()
                        save_name     = key.atoms['name'].tolist()
                        save_resid    = key.atoms['resid'].tolist()
                        save_chain    = np.char.add(key.atoms['chain'].tolist(), [str(i)] * len(positions))
        
                    else:
                        raise IOError('input protein structure is not found')
        
                    finalpositions = positions + proteinP[i]
                    proteindata['chain'].extend(save_chain)
                    proteindata['resname'].extend(save_resname)
                    proteindata['resid'].extend(save_resid)
                    proteindata['name'].extend(save_name)
                    proteindata['x'].extend(finalpositions[:,0])
                    proteindata['y'].extend(finalpositions[:,1])
                    proteindata['z'].extend(finalpositions[:,2])
        
                proteinU = Universe(data=proteindata)
                proteinU.atoms = proteinU.atoms.astype({'resid': int})

            else:
                raise IOError('input protein structure is not found')

            Hatoms   = proteinU.atoms.name.str.startswith('H')
            proteinU.atoms.bfactor = 0.0
            proteinU.atoms.loc[~Hatoms, 'bfactor'] = 1.0
            #this is a memory-intensive job
            #if len(proteinU.bonds) == 0: proteinU.addBondFromDistance()
            proteinU.write(workdir + '/protein.dms', wrap=False)
            proteinU.write(workdir + '/protein.pdb', wrap=False)
  
       
        ### Read Martini
        #martiniff = ReadMartini(ff=martini, ff_add=martini_add, define={'FLEXIBLE': 'True'})
        martiniff = ReadMartini(ff=martini, ff_add=martini_add, constraint_to_bond=True, Kc2b=10000.0, addnbtype=addnbtype)

        ### Mapping & Translate
        if protein: Map(workdir + '/protein.dms', workdir + '/protein_CG.dms', add_notavailableAAAtoms=True, mapping=mapping, mapping_add=mapping_add)
        
        ### Construct a bilayer
        instance = Bilayer(protein=workdir + '/protein_CG.dms' if protein else None, 
                           upper=upper, lower=lower, between=between, sep=sep,
                           waterz=water, rcut=rcut, rcut_use_enclosed_protein=rcut_use_enclosed_protein,
                           out=workdir + '/step1.bilayer.dms', 
                           martini=martiniff, 
                           lipidpath=lipidpath,
                           hydrophobic_thickness=hydrophobic_thickness,
                           mode=mode, dN=dN, dx=dx)

        ### Solvate
        if solvate:
            usol = SolvateMartini(workdir + '/step1.bilayer.dms', removedr=removedr, conc=0.0, waterchain='6')
            cell = usol.cell
            dim  = usol.dimensions
            bA1  = usol.atoms.name == 'W'
            bA2  = usol.atoms.z    <  hydrophobic_thickness / 2 + sep/2 + 10
            bA3  = usol.atoms.z    > -hydrophobic_thickness / 2 - sep/2 - 10
            usol = Universe(data=usol.atoms[~(bA1 & bA2 & bA3)])
        
            qtot = 0
            for resn, value in usol.atoms.groupby('resn'):
                resname = value.resname.values[0]
                if resname == 'W': continue
                if resname in martiniff.martini['molecules']:
                    qtot += np.sum(martiniff.martini['molecules'][value.resname.values[0]]['atoms']['q'])
                else:
                    qtot += 0

            usol = ionize(usol, conc=ionconc, posionchain='4', negionchain='5', qtot=qtot)
            usol.dimensions = dim
            usol.cell       = cell
 
            #usol = SolvateMartini(workdir + '/step1.bilayer.dms', removedr=removedr, conc=ionconc)
            #if protein:
            #    usol = SolvateMartini(workdir + '/step1.bilayer.dms', removedr=removedr, conc=ionconc)
            #else:
            #    usol = SolvateMartini(workdir + '/step1.bilayer.dms', removedr=removedr, conc=ionconc, membrane=True)

        else:
            usol = instance.universe
            barfreq=0 # make it NVT
            remove_solvent=False # remove_solvent is redundant

        ### Translate
        shift = usol.dimensions[0:3] / 2
        usol.atoms[['x','y','z']] += shift
        usol.write(workdir + '/step1.sol.dms')
        usol.write(workdir + '/step1.sol.pdb')

        ### Martinize
        MartinizeDMS(workdir + '/step1.sol.dms', out = workdir + '/step1.martini.dms', martini=martiniff, addnbtype=addnbtype[0])
        dumpsql(workdir + '/step1.martini.dms')

        ### Create system & add posz to GL1 and GL2
        # dms = DMSFile(workdir + '/step1.martini.dms')
        # martiniU = Universe(workdir + '/step1.martini.dms')
        # martiniU.atoms.loc[((martiniU.atoms.name == 'GL1')), 'bfactor'] = 1.0
        # positional restraints and barostat do not like each other....
        # https://github.com/openmm/openmm/issues/1854
        #posre = addPosrePeriodicZ(martiniU, bfactor_posre=0.5, k=fc)
        #dms.createSystem(REM=False, tapering=tapering, martini=True, nonbondedCutoff=nonbondedCutoff, nonbondedMethod='CutoffPeriodic', improper_prefactor=improper_prefactor, removeCMMotion=False)
        #dms.system.addForce(posre)

        #dms = DMSFile(workdir + '/step1.martini.dms')
        #dms.createSystem(REM=REM, tapering=tapering, martini=True, nonbondedCutoff=nonbondedCutoff, nonbondedMethod='CutoffPeriodic', improper_prefactor=improper_prefactor, removeCMMotion=True)
        #dms.runEMNPT(dt=0.002, nsteps=cg_nsteps_2fs, frictionCoeff=frictionCoeff, barfreq=barfreq, T=T, semiisotropic=True)
        #dms.runEMNPT(dt=0.005, nsteps=cg_nsteps_5fs, frictionCoeff=frictionCoeff, barfreq=0,       T=T, EM=False)
        #dms.runEMNPT(dt=dt,    nsteps=cg_nsteps,     frictionCoeff=frictionCoeff, barfreq=0,       T=T, EM=False, out=workdir + '/step2.dms')

        #dms = DMSFile(workdir + '/step1.martini.dms')
        #dms.createSystem(REM=REM, tapering=tapering, martini=True, nonbondedCutoff=nonbondedCutoff, nonbondedMethod='CutoffPeriodic', improper_prefactor=improper_prefactor, removeCMMotion=True)
        #dms.runEMNPT(dt=0.002, nsteps=cg_nsteps_2fs, frictionCoeff=frictionCoeff, barfreq=barfreq, T=T, semiisotropic=True)
        #dms.runEMNPT(dt=0.005, nsteps=cg_nsteps_5fs, frictionCoeff=frictionCoeff, barfreq=0,       T=T, EM=False)
        #dms.runEMNPT(dt=dt,    nsteps=cg_nsteps,     frictionCoeff=frictionCoeff, barfreq=0,       T=T, EM=False, out=workdir + '/step2.dms')

        ### REM BEFORE NPT
        #dms.createSystem(REM=True, tapering=tapering, martini=True, nonbondedCutoff=nonbondedCutoff, nonbondedMethod='CutoffPeriodic', improper_prefactor=improper_prefactor, removeCMMotion=False)
        #dms.system.addForce(posre)
        #dms.runEMNPT(workdir + '/step2.rem.dms', emout=None, dt=dt, nsteps=0, frictionCoeff=frictionCoeff, barfreq=0, T=T)
        
        martiniU = Universe(workdir + '/step1.martini.dms')
        #martiniU.atoms.loc[((martiniU.atoms.name == 'GL1')), 'bfactor'] = 2.0
        martiniU.atoms.loc[((martiniU.atoms.name == 'PO4')), 'bfactor'] = 2.0

        dms = DMSFile(workdir + '/step1.martini.dms')
        dms.createSystem(REM=True,  tapering=tapering, martini=True, nonbondedCutoff=nonbondedCutoff, nonbondedMethod='CutoffPeriodic', improper_prefactor=improper_prefactor, removeCMMotion=True)
        dms.system.addForce(addPosre(martiniU, bfactor_posre=1.5, fcx=0.0, fcy=0.0, fcz=fc))
        #if protein: dms.system.addForce(addPosre(martiniU, bfactor_posre=1.5, fcx=0.0, fcy=0.0, fcz=fc))
        #dms.system.addForce(addFlatBottomZ(martiniU, bfactor_posre=0.5, radius=hydrophobic_thickness/2, rfb=0.1, R0z=shift[2], fc=fc, chain='UPPER'))
        #dms.system.addForce(addFlatBottomZ(martiniU, bfactor_posre=0.5, radius=hydrophobic_thickness/2, rfb=0.1, R0z=shift[2], fc=fc, chain='LOWER'))
        dms.runEMNPT(dt=dt_rem, nsteps=int(cg_nsteps_rem), frictionCoeff=frictionCoeff, barfreq=barfreq, dcdfreq=dcdfreq, csvfreq=csvfreq, tension=tension, P=P, T=T, semiisotropic=True, out=workdir + '/step2.rem.dms')

        #dms = DMSFile(workdir + '/step2.rem.dms')
        dms.createSystem(REM=False, tapering=tapering, martini=True, nonbondedCutoff=nonbondedCutoff, nonbondedMethod='CutoffPeriodic', improper_prefactor=improper_prefactor, removeCMMotion=True)
        if protein: dms.system.addForce(addPosre(martiniU, bfactor_posre=1.5, fcx=0.0, fcy=0.0, fcz=fc2))
        dms.runEMNPT(dt=dt, nsteps=int(cg_nsteps), frictionCoeff=frictionCoeff, barfreq=barfreq, dcdfreq=dcdfreq, csvfreq=csvfreq, tension=tension, P=P, T=T, semiisotropic=True, out=workdir + '/step2.dms')
        
        ### Translate Back
        u = Universe(workdir + '/step2.dms')

        bA1  = u.atoms.resname == 'W'
        bA2  = u.atoms.chain   == '4'
        bA3  = u.atoms.chain   == '5'
        bA   = bA1 | bA2 | bA3
        NtotalW = len(u.atoms[bA1])

        if protein:
            membranecenter = shift
        else:
            membranecenter = u.atoms[~bA][['x','y','z']].mean(axis=0).to_numpy()

        u.atoms[['x','y','z']] -= membranecenter


        ### APL
        Ntotal = instance.upperN + instance.lowerN
        APL = u.dimensions[0] * u.dimensions[1] / (Ntotal / 2)
        WL  = NtotalW * 4 / Ntotal
        print(f'APL: {APL:.2f} A^2')
        print(f'W/L: {WL:.2f}')

        if remove_solvent:
            dimensions = u.dimensions
            cell = u.cell
            bA1  = u.atoms.resname == 'W'
            bA2  = u.atoms.chain   == '4'
            bA3  = u.atoms.chain   == '5'
            bA   = bA1 | bA2 | bA3
            u    = Universe(data=u.atoms[~bA])
            u.dimensions = dimensions
            u.cell = cell
        
        u.sort()
        u.write(workdir + '/step3.dms', wrap=True)
        u.write(workdir + '/step3.pdb', wrap=True)
        #u.atoms[['x','y','z']] -= u.dimensions[0:3] / 2
        #u.write(workdir + '/step3.dms')
        #u.write(workdir + '/step3.pdb')

        if protein:
            #proteinshifted = Universe(workdir + '/protein.pdb')
            #proteinshifted.atoms[['x','y','z']] += shift - u.dimensions[0:3]/2
            #proteinshifted.write(workdir + '/step8_protein.pdb')

            noprotein = Universe(data=u.atoms[~u.atoms.resname.isin(three2one.keys())])
            noprotein.dimensions = u.dimensions
            noprotein.cell = u.cell
            noprotein.write(workdir + '/step3.noprotein.dms')
            Backmap(workdir + '/step3.noprotein.dms', workdir=workdir, use_existing_workdir=True, nsteps=int(aa_nsteps),
                    AA=workdir + '/protein.dms', fileindex=4, mapping=mapping, mapping_add=mapping_add, ff=ff, ff_add=ff_add,
                    use_AA_structure=use_AA_structure, AA_shrink_factor=AA_shrink_factor, rockCtype=rockCtype, rockHtype=rockHtype, T=T, water_chain='6789',
                    changename=changename)
        else:
            Backmap(workdir + '/step3.dms', workdir=workdir, use_existing_workdir=True, nsteps=int(aa_nsteps),
                    AA=None, fileindex=4, mapping=mapping, mapping_add=mapping_add, ff=ff, ff_add=ff_add,
                    use_AA_structure=use_AA_structure, AA_shrink_factor=AA_shrink_factor, rockCtype=rockCtype, rockHtype=rockHtype, T=T, water_chain='6789',
                    changename=changename)


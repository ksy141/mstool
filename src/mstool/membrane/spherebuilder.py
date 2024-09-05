import os
import shutil
import numpy as np
from  .sphere                import Sphere
from  .sphereprotein         import SphereProtein
from  ..core.map             import Map
from  ..core.backmap         import Backmap
from  ..core.readmartini     import ReadMartini
from  ..core.solvate_martini import SolvateMartini, ionize
from  ..core.martinizedms    import MartinizeDMS
from  ..core.dmsfile         import DMSFile
from  ..core.universe        import Universe
from  ..core.backmap         import Backmap
from  ..utils.dump           import dumpsql
from  ..utils.openmmutils    import addSpherePosre
from  ..utils.protein_sel    import one2three, three2one
pwd = os.path.dirname(os.path.realpath(__file__))

class SphereBuilder:
    def __init__(self, radius, workdir='workdir', protein=None, upper={}, lower={}, between={}, sep=0.0,
                 water=50.0, rcut=4.0, rockCtype='CTL3', rockHtype='HAL3',
                 martini=[], martini_add=[], lipidpath=pwd+'/../../../FF/martini2.2/structures/',
                 mapping=[], mapping_add=[],
                 ff=[], ff_add=[],
                 removedr=4.5, aa_nsteps=0, fc=10.0,
                 dt_rem=0.025, cg_nsteps_rem=100000,
                 dt=0.020, cg_nsteps=100000,
                 frictionCoeff=10.0, barfreq=10, nonbondedCutoff=1.1, 
                 dcdfreq=1000, csvfreq=1000, P=1.0, T=310,
                 improper_prefactor=0.99, use_existing_folder=False,
                hydrophobic_thickness=30, ionconc=0.15,
                 use_AA_structure=True,
                 alpha=0.0, beta=0.0, gamma=0.0,
                 remove_solvent=False,
                 solvate=True,
                 tapering='shift',
                 changename={}, addnbtype=['ZCARBON', 0.34, 1.51]):

        '''Sphere builder.
        Parameters
        ----------
        radius: float
            Provide the radius of the sphere
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
        water: float=50.0
            water thickness (A).
        rcut: float=3.0
            lipids that are closer than rcut (A) to protein will be removed (mode=remove)
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
        cg_nsteps: int=1000
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
        
        ### protein
        if protein:
            if isinstance(protein, str):
                proteinU = Universe(protein)

            elif isinstance(protein, dict):
                SphereProtein(radius=radius, protein=protein, out=workdir + '/protein.dms')
                proteinU = Universe(workdir + '/protein.dms')

            Hatoms   = proteinU.atoms.name.str.startswith('H')
            proteinU.atoms.bfactor = 0.0
            proteinU.atoms.loc[~Hatoms, 'bfactor'] = 1.0
            if len(proteinU.bonds) == 0: proteinU.addBondFromDistance()
            proteinU.write(workdir + '/protein.dms', wrap=False)
            proteinU.write(workdir + '/protein.pdb', wrap=False)
       
        ### Read Martini
        #martiniff = ReadMartini(ff=martini, ff_add=martini_add, define={'FLEXIBLE': 'True'})
        martiniff = ReadMartini(ff=martini, ff_add=martini_add, constraint_to_bond=True, Kc2b=10000.0, addnbtype=addnbtype)

        ### Mapping & Translate
        if protein: Map(workdir + '/protein.dms', workdir + '/protein_CG.dms', add_notavailableAAAtoms=True, mapping=mapping, mapping_add=mapping_add)

        ### Construct a bilayer
        instance = Sphere( protein=workdir + '/protein_CG.dms' if protein else None, 
                           upper=upper, lower=lower, between=between, sep=sep,
                           water=water, rcut=rcut,
                           out=workdir + '/step1.bilayer.dms', 
                           martini=martiniff, 
                           lipidpath=lipidpath,
                           hydrophobic_thickness=hydrophobic_thickness,
                           alpha=alpha, beta=beta, gamma=gamma, r=radius)

        ### Solvate
        if solvate:
            #usol = SolvateMartini(workdir + '/step1.bilayer.dms', removedr=removedr, conc=ionconc, 
            #                      posionchain='4', negionchain='5', waterchain='6')

            usol = SolvateMartini(workdir + '/step1.bilayer.dms', removedr=removedr, conc=0.0, waterchain='6')
            cell = usol.cell
            dim  = usol.dimensions
            bA1  = usol.atoms.name == 'W'
            bA2  = usol.atoms.x ** 2 + usol.atoms.y ** 2 + usol.atoms.z ** 2 < (radius + hydrophobic_thickness/2 + sep/2) ** 2
            bA3  = usol.atoms.x ** 2 + usol.atoms.y ** 2 + usol.atoms.z ** 2 > (radius - sep/2) ** 2
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
        
        else:
            usol = instance.universe
            # make it NVT
            barfreq=0
            # remove_solvent is redundant
            remove_solvent=False

        ### Translate
        shift = usol.dimensions[0:3] / 2
        usol.atoms[['x','y','z']] += shift
        usol.write(workdir + '/step1.sol.dms')
        usol.write(workdir + '/step1.sol.pdb')

        ### Martinize
        MartinizeDMS(workdir + '/step1.sol.dms', out = workdir + '/step1.martini.dms', martini=martiniff, addnbtype=addnbtype[0])
        dumpsql(workdir + '/step1.martini.dms')

        ### Create system & add posz to GL1 and GL2
        #dms = DMSFile(workdir + '/step1.martini.dms')
        #try:
        #    dms.createSystem(REM=False, tapering=tapering, martini=True, nonbondedCutoff=nonbondedCutoff, nonbondedMethod='CutoffPeriodic', improper_prefactor=improper_prefactor, removeCMMotion=False)
        #    martiniU = Universe(workdir + '/step1.martini.dms')
        #    martiniU.atoms.loc[((martiniU.atoms.name == 'GL1')), 'bfactor'] = 1.0
        #    dms.system.addForce(addSpherePosre(martiniU, bfactor_posre=0.5, 
        #                                       radius=radius + hydrophobic_thickness/2,
        #                                       rfb=0.1,
        #                                       R0=shift,
        #                                       fc=fc,
        #                                       chain='UPPER'))

        #    dms.system.addForce(addSpherePosre(martiniU, bfactor_posre=0.5, 
        #                                       radius=radius - hydrophobic_thickness/2,
        #                                       rfb=0.1,
        #                                       R0=shift,
        #                                       fc=fc,
        #                                       chain='LOWER'))
        #    dms.runEMNPT(workdir + '/step2.dms', emout=workdir + '/step2.em.dms', dt=dt, nsteps=cg_nsteps, frictionCoeff=frictionCoeff, barfreq=barfreq, semiisotropic=False, T=T)
        #except:
        #### REM
        #dms.createSystem(REM=True, tapering=tapering, martini=True, nonbondedCutoff=nonbondedCutoff, nonbondedMethod='CutoffPeriodic', improper_prefactor=improper_prefactor, removeCMMotion=False)
        #martiniU = Universe(workdir + '/step1.martini.dms')
        #martiniU.atoms.loc[((martiniU.atoms.name == 'GL1')), 'bfactor'] = 1.0
        #dms.system.addForce(addSpherePosre(martiniU, bfactor_posre=0.5, 
        #                                   radius=radius + hydrophobic_thickness/2,
        #                                   rfb=0.1,
        #                                   R0=shift,
        #                                   fc=fc,
        #                                   chain='UPPER'))

        #dms.system.addForce(addSpherePosre(martiniU, bfactor_posre=0.5, 
        #                                   radius=radius - hydrophobic_thickness/2,
        #                                   rfb=0.1,
        #                                   R0=shift,
        #                                   fc=fc,
        #                                   chain='LOWER'))
        #dms.runEMNPT(workdir + '/step2.tmp.dms', dt=dt, nsteps=0, frictionCoeff=frictionCoeff, T=T)
        
        ### NPT
        martiniU = Universe(workdir + '/step1.martini.dms')
        #martiniU.atoms.loc[((martiniU.atoms.name == 'GL1')), 'bfactor'] = 2.0
        martiniU.atoms.loc[((martiniU.atoms.name == 'PO4')), 'bfactor'] = 2.0

        dms = DMSFile(workdir + '/step1.martini.dms')
        dms.createSystem(REM=True,  tapering=tapering, martini=True, nonbondedCutoff=nonbondedCutoff, nonbondedMethod='CutoffPeriodic', improper_prefactor=improper_prefactor, removeCMMotion=True)
        dms.system.addForce(addSpherePosre(martiniU, bfactor_posre=1.5, radius=radius + hydrophobic_thickness/2 + sep/2, rfb=0.1, R0=shift, fc=fc, chain='0'))
        dms.system.addForce(addSpherePosre(martiniU, bfactor_posre=1.5, radius=radius - hydrophobic_thickness/2 - sep/2, rfb=0.1, R0=shift, fc=fc, chain='1'))
        dms.runEMNPT(dt=dt_rem, nsteps=int(cg_nsteps_rem), frictionCoeff=frictionCoeff, barfreq=barfreq, dcdfreq=dcdfreq, csvfreq=csvfreq, P=P, T=T, semiisotropic=False, out=workdir + '/step2.rem.dms')

        dms.createSystem(REM=False, tapering=tapering, martini=True, nonbondedCutoff=nonbondedCutoff, nonbondedMethod='CutoffPeriodic', improper_prefactor=improper_prefactor, removeCMMotion=True)
        dms.system.addForce(addSpherePosre(martiniU, bfactor_posre=1.5, radius=radius + hydrophobic_thickness/2 + sep/2, rfb=0.1, R0=shift, fc=fc, chain='0'))
        dms.system.addForce(addSpherePosre(martiniU, bfactor_posre=1.5, radius=radius - hydrophobic_thickness/2 - sep/2, rfb=0.1, R0=shift, fc=fc,chain='1'))
        dms.runEMNPT(dt=dt, nsteps=int(cg_nsteps), frictionCoeff=frictionCoeff, barfreq=barfreq, dcdfreq=dcdfreq, csvfreq=csvfreq, P=P, T=T, semiisotropic=False, out=workdir + '/step2.dms')


        #dms.runEMNPT(dt=0.002, nsteps=cg_nsteps_2fs, frictionCoeff=frictionCoeff, barfreq=barfreq, T=T, semiisotropic=False)
        #dms.runEMNPT(dt=0.005, nsteps=cg_nsteps_5fs, frictionCoeff=frictionCoeff, barfreq=0,       T=T, EM=False)
        #dms.runEMNPT(dt=dt,    nsteps=cg_nsteps,     frictionCoeff=frictionCoeff, barfreq=0,       T=T, EM=False, out=workdir + '/step2.dms')
        
        ### Translate Back
        u = Universe(workdir + '/step2.dms')
        u.atoms[['x','y','z']] -= shift

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

        if protein:
            noprotein = Universe(data=u.atoms[~u.atoms.resname.isin(three2one.keys())])
            noprotein.dimensions = u.dimensions
            noprotein.cell = u.cell
            noprotein.write(workdir + '/step3.noprotein.dms')
            Backmap(workdir + '/step3.noprotein.dms', workdir=workdir, use_existing_workdir=True, nsteps=int(aa_nsteps),
                    AA=workdir + '/protein.dms', fileindex=4, mapping=mapping, mapping_add=mapping_add, ff=ff, ff_add=ff_add,
                    use_AA_structure=use_AA_structure, rockCtype=rockCtype, rockHtype=rockHtype, T=T, water_chain='6789',
                    changename=changename)
        else:
            Backmap(workdir + '/step3.dms', workdir=workdir, use_existing_workdir=True, nsteps=int(aa_nsteps),
                    AA=None, fileindex=4, mapping=mapping, mapping_add=mapping_add, ff=ff, ff_add=ff_add,
                    use_AA_structure=use_AA_structure, rockCtype=rockCtype, rockHtype=rockHtype, T=T, water_chain='6789',
                    changename=changename)


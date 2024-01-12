import os
import shutil
from  .bilayer               import Bilayer
from  ..core.map             import Map
from  ..core.backmap         import Backmap
from  ..core.readmartini     import ReadMartini
from  ..core.solvate_martini import SolvateMartini
from  ..core.martinizedms    import MartinizeDMS
from  ..core.dmsfile         import DMSFile
from  ..core.universe        import Universe
from  ..core.backmap         import Backmap
from  ..utils.dump           import dumpsql
from  ..utils.openmmutils    import addPosre
pwd = os.path.dirname(os.path.realpath(__file__))

class BilayerBuilder:
    def __init__(self, workdir='workdir', protein=None, upper={}, lower={}, dx=8.0, waterz=30.0, rcut=3.0, 
                 mode='shift', dN=5, 
                 martini=[], martini_add=[], lipidpath=pwd+'/../../../FF/martini2.2/structures/',
                 mapping=[], mapping_add=[],
                 ff=[], ff_add=[],
                 removedr=8.5, dt=0.007, nsteps=50000, frictionCoeff=5.0, barfreq=100, nonbondedCutoff=1.1, 
                 improper_prefactor=0.99, use_existing_folder=False):
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
        waterz: float=30.0
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
        dt: float=0.007
            Integration time for coarse-grained simulations.
            dt=0.020 works fine for no virtual sites including systems.
            dt=0.007 is the largest value for virtual sites containing systems (e.g., cholesterol, CHL1)
            If you encounter numerical issue (Energy is Nan, or particles are all over the place), reduce dt while increasing nsteps.
        nsteps: int=50000
            Number of coarse-grained simulation steps.
            If your membrane does not look equilibrated (e.g., water not equally dispersed in solution), increase this number.
        frictionCoeff: float=5.0
            Friction coefficient for coarse-grained simulations.
        barfreq: int=100
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
        '''

        
        ### workdir
        if (not use_existing_folder) or (use_existing_folder and not os.path.exists(workdir)):
            os.mkdir(workdir)

        if protein: Universe(protein).write(workdir + '/protein.dms', wrap=False)
       
        ### Read Martini
        martiniff = ReadMartini(ff=martini, ff_add=martini_add, define={'FLEXIBLE': 'True'})

        ### Mapping & Translate
        #if protein: Map(workdir + '/protein.dms', workdir + '/protein_CG.dms', add_notavailableAAAtoms=True)
        if protein: Map(workdir + '/protein.dms', workdir + '/protein_CG.dms', add_notavailableAAAtoms=False)

        ### Construct a bilayer
        u = Bilayer(protein=workdir + '/protein_CG.dms', 
                    upper=upper, lower=lower, 
                    dx=dx, waterz=waterz, rcut=rcut,
                    out=workdir + '/step1.bilayer.dms', 
                    mode=mode, dN=dN, martini=martiniff, 
                    lipidpath=lipidpath)

        ### Solvate
        usol = SolvateMartini(workdir + '/step1.bilayer.dms', removedr=removedr)

        ### Translate
        shift = usol.dimensions[0:3] / 2
        usol.atoms[['x','y','z']] += shift
        usol.atoms.loc[(usol.atoms.name == 'GL1') | (usol.atoms.name == 'GL2'), 'bfactor'] = 1.0
        usol.write(workdir + '/step1.sol.dms')

        ### Martinize
        MartinizeDMS(workdir + '/step1.sol.dms', out = workdir + '/step1.martini.dms', martini=martiniff)
        dumpsql(workdir + '/step1.martini.dms')

        ### Run
        dms = DMSFile(workdir + '/step1.martini.dms')
        dms.createSystem(REM=False, martini=True, nonbondedCutoff=nonbondedCutoff, nonbondedMethod='CutoffPeriodic', improper_prefactor=improper_prefactor, removeCMMotion=False)
        #dms.system.addForce(addPosre(Universe(workdir + '/step1.martini.dms'), bfactor_posre=0.5, fcx=0.0, fcy=0.0, fcz=100.0))
        dms.runEMNPT(workdir + '/step2.dms', dt=dt, nsteps=nsteps, frictionCoeff=frictionCoeff, barfreq=barfreq)
        
        ### Translate Back
        u = Universe(workdir + '/step2.dms')
        u.atoms[['x','y','z']] -= shift
        u.write(workdir + '/step3.dms', wrap=True)

        ### Backmapping
        Backmap(workdir + '/step3.dms', workdir=workdir, use_existing_workdir=True,
                AA=protein, fileindex=4, mapping=mapping, mapping_add=mapping_add, ff=ff, ff_add=ff_add)


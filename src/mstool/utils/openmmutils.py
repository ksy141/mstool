import sqlite3
import shutil
import copy
import numpy as np
from openmm.unit       import *
from openmm.app        import *
from openmm            import *
from .protein_sel      import three2one
from ..core.universe   import Universe, Merge
from ..core.dms2openmm import DMS2openmm
from ..core.readxml    import ReadXML

def getEnergy(simulation):
    return simulation.context.getState(getEnergy=True).getPotentialEnergy()._value

def getCell(simulation):
    return simulation.context.getState().getPeriodicBoxVectors(asNumpy=True)._value * 10

def getPositions(simulation):
    return simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions(asNumpy=True)._value * 10

def Pressure(P=1.0, T=310.0, barfreq=100):
    return MonteCarloBarostat(P*bar, T*kelvin)

def MembranePressure(P=1.0, T=310.0, r=0.0, barfreq=100):
    return MonteCarloMembraneBarostat(P*bar, r*bar*nanometer, T*kelvin,
             MonteCarloMembraneBarostat.XYIsotropic, MonteCarloMembraneBarostat.ZFree)

def getDMSEnergy(dms_in, nonbondedCutoff=1.1, nonbondedMethod='CutoffPeriodic'):

    system, dms = DMS2openmm(
            dms_in          = dms_in,
            nonbondedMethod = nonbondedMethod,
            nonbondedCutoff = nonbondedCutoff,
            soft            = False).make()

    integrator = LangevinMiddleIntegrator(310*kelvin, 1/picosecond, 0.002*picoseconds)
    simulation = Simulation(dms.topology, system, integrator)
    simulation.context.setPositions(DesmondDMSFile(dms_in).positions)

    PE = getEnergy(simulation)
    print('{:7s}: {:10.3f} kJ/mol'.format('openMM', PE))
    return PE

def runEM(structure, forcefield, out=None, nonbondedCutoff=1.1, nonbondedMethod='CutoffNonPeriodic', addForces=[]):
    '''
    forcefield:
    charmm36.xml
    amber99sbildn.xml
    '''
    if structure.split('.')[-1] == 'dms':
        pdb = DesmondDMSFile(structure)
    if structure.split('.')[-1] == 'pdb':
        pdb = PDBFile(structure)
    
    if isinstance(forcefield, list):
        forcefield = ForceField(*forcefield)
    else:
        forcefield = ForceField(forcefield)
    system = forcefield.createSystem(pdb.topology)

    ### ADD additional forces
    for addf in addForces:
        system.addForce(copy.deepcopy(addf))

    integrator = LangevinMiddleIntegrator(310*kelvin, 1/picosecond, 0.002*picoseconds)
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    platform = simulation.context.getPlatform().getName()
    print('-------------------------------')
    print("Platform: ", platform)
    print('E0: %.3e kJ/mol' %getEnergy(simulation))
    simulation.minimizeEnergy()
    print('E1: %.3e kJ/mol' %getEnergy(simulation))
    print('-------------------------------')
    ### SAVE THE LAST FRAME
    #pos = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)._value * 10
    pos = getPositions(simulation)
    
    if out:
        u = Universe(structure)
        u.atoms[['x','y','z']] = pos
        u.write(out)

#def runEMNPT(structure, forcefield, out=None, nonbondedCutoff=1.1, nonbondedMethod='CutoffPeriodic',
#             nsteps=10000, dcdfreq=1000, csvfreq=1000, dt=0.002, P=1.0, T=310, semiisotropic=False,
#             addForces=[], barfreq=100, frictionCoeff=1.0):
#    '''
#    forcefield:
#    charmm36.xml
#    amber99sbildn.xml
#    '''
#    if structure.split('.')[-1] == 'dms':
#        pdb = DesmondDMSFile(structure)
#    if structure.split('.')[-1] == 'pdb':
#        pdb = PDBFile(structure)
#    
#    if isinstance(forcefield, list):
#        forcefield = ForceField(*forcefield)
#    else:
#        forcefield = ForceField(forcefield)
#
#    if nonbondedMethod == 'CutoffPeriodic':
#        nonbondedMethod = CutoffPeriodic
#    else:
#        nonbondedMethod = CutoffNonPeriodic
#
#    system = forcefield.createSystem(pdb.topology, nonbondedMethod=nonbondedMethod)
#
#    ### ADD BAROSTAT
#    if barfreq > 0:
#        if semiisotropic:
#            addForces.append(MembranePressure(P=P, T=T, barfreq=barfreq))
#        else:
#            addForces.append(Pressure(P=P, T=T, barfreq=barfreq))
#
#    ### ADD additional forces
#    for addf in addForces:
#        system.addForce(copy.deepcopy(addf))
#    
#    ### PREPARE SIMS
#    integrator = LangevinMiddleIntegrator(T*kelvin, frictionCoeff/picosecond, dt*picoseconds)
#    simulation = Simulation(pdb.topology, system, integrator)
#    simulation.context.setPositions(pdb.positions)
#    platform = simulation.context.getPlatform().getName()
#    print('-------------------------------')
#    print("Platform: ", platform)
#    print('E0: %.3e kJ/mol' %getEnergy(simulation))
#    simulation.minimizeEnergy()
#    print('E1: %.3e kJ/mol' %getEnergy(simulation))
#    print('-------------------------------')
#
#    ### RUN NPT
#    prefix = '.'.join(out.split('.')[:-1])
#    simulation.reporters.append(DCDReporter(      prefix + '.dcd', dcdfreq))
#    simulation.reporters.append(StateDataReporter(prefix + '.csv', csvfreq, step=True, potentialEnergy=True, temperature=True))
#    simulation.step(nsteps)
#
#    cell = getCell(simulation)
#    dimensions = [cell[0][0], cell[1][1], cell[2][2], 90, 90, 90]
#
#    ### SAVE THE LAST FRAME
#    pos = getPositions(simulation)
#    
#    if out:
#        u = Universe(structure)
#        u.atoms[['x','y','z']] = pos
#        u.cell = cell
#        u.dimensions = dimensions
#        u.write(out)

def runEMNPT(structure, ff=[], ff_add=[], out=None, nonbondedCutoff=1.2, nonbondedMethod='CutoffPeriodic',
             nsteps=10000, dcdfreq=1000, csvfreq=1000, dt=0.002, P=1.0, T=310, tension=0.0, semiisotropic=False,
             addForces=[], barfreq=100, frictionCoeff=1.0, EM=True, shift=False):
    
    ### Structure -> Modeller
    if not isinstance(structure, list):
        structure = [structure]
    
    modeller_combined = []
    universe_combined = []
    for struct in structure:
        if struct.split('.')[-1] == 'dms':
            pdb = DesmondDMSFile(struct)
        if struct.split('.')[-1] == 'pdb':
            pdb = PDBFile(struct)
        modeller_combined.append([pdb.topology, pdb.positions.in_units_of(nanometer)])
        realpbc = pdb.topology.getPeriodicBoxVectors().in_units_of(nanometer)
        universe_combined.append(Universe(struct).atoms)
    u = Merge(*universe_combined)
    
    # pdb.positions: Quantity([Vec3(0,0,0), Vec3(1,1,1)], unit=nanometer)
    # shiftpbc:      Quantity(Vec3(10, 10, 10), unit=nanometer)

    if shift:
        shiftpbc = realpbc[0]/2 + realpbc[1]/2 + realpbc[2]/2
    else:
        shiftpbc = Quantity(value=Vec3(x=0, y=0, z=0), unit=nanometer)

    modeller = Modeller(modeller_combined[0][0], 
                        Quantity([pos.value_in_unit(nanometer) + shiftpbc.value_in_unit(nanometer) \
                                  for pos in modeller_combined[0][1]], unit=nanometer))

    for i in range(1, len(modeller_combined)):
        modeller.add(modeller_combined[i][0],
                     Quantity([pos.value_in_unit(nanometer) + shiftpbc.value_in_unit(nanometer) \
                               for pos in modeller_combined[i][1]], unit=nanometer))

    modeller.topology.setPeriodicBoxVectors(realpbc)

    unique_bonds = set()
    for bond in modeller.topology.bonds():
        i0, i1 = bond[0].index, bond[1].index
        if i0 > i1:
            i0, i1 = i1, i0
        unique_bonds.add((i0, i1))
    u.bonds = [list(bond) for bond in unique_bonds]

    ### FF -> CreateSystem
    if not isinstance(ff_add, list):
        ff_add = [ff_add]
    
    xml = ReadXML(ff=ff, ff_add=ff_add)
    forcefield = ForceField(*xml.ff)

    if nonbondedMethod == 'CutoffPeriodic':
        nonbondedMethod = CutoffPeriodic
    else:
        nonbondedMethod = CutoffNonPeriodic
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=nonbondedMethod)

    ### ADD BAROSTAT
    if barfreq > 0:
        if semiisotropic:
            addForces.append(MembranePressure(P=P, T=T, r=tension, barfreq=barfreq))
        else:
            addForces.append(Pressure(P=P, T=T, barfreq=barfreq))

    ### ADD additional forces
    for addf in addForces:
        system.addForce(copy.deepcopy(addf))
    
    ### PREPARE SIMS
    integrator = LangevinMiddleIntegrator(T*kelvin, frictionCoeff/picosecond, dt*picoseconds)
    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)
    platform = simulation.context.getPlatform().getName()
    print('-------------------------------')
    print("Platform: ", platform)
    print('E0: %.3e kJ/mol' %getEnergy(simulation))
    if EM: simulation.minimizeEnergy()
    print('E1: %.3e kJ/mol' %getEnergy(simulation))
    print('-------------------------------')

    ### RUN NPT
    prefix = '.'.join(out.split('.')[:-1])
    simulation.reporters.append(DCDReporter(      prefix + '.dcd', dcdfreq))
    simulation.reporters.append(StateDataReporter(prefix + '.csv', csvfreq, step=True, potentialEnergy=True, temperature=True))
    simulation.step(nsteps)

    cell = getCell(simulation)
    dimensions = [cell[0][0], cell[1][1], cell[2][2], 90, 90, 90]

    ### SAVE THE LAST FRAME
    pos = getPositions(simulation)
    
    if out:
        u.atoms[['x','y','z']] = pos
        u.cell = cell
        u.dimensions = dimensions
        u.write(out)


def runMartiniEM(dms_in, out, pos_in=None, soft=False, A=200, C=50,
    nonbondedCutoff=1.1, nonbondedMethod='CutoffPeriodic', T=310, dt=0.002,
    addForces=[]):
    '''
    CutoffNonPeriodic
    CutoffPeriodic
    NoCutoff
    '''
    system, dms = DMS2openmm(
            dms_in          = dms_in,
            nonbondedMethod = nonbondedMethod,
            nonbondedCutoff = nonbondedCutoff,
            soft            = soft,
            A               = A,
            C               = C).make()
    
    ### ADD additional forces
    for addf in addForces:
        system.addForce(copy.deepcopy(addf))

    ### PREPARE SIMS
    integrator = LangevinMiddleIntegrator(T*kelvin, 1/picosecond, dt*picoseconds)
    simulation = Simulation(dms.topology, system, integrator)

    ### SET POSITIONS
    if pos_in:
        if pos_in.split('.')[-1] == 'dms':
            simulation.context.setPositions(DesmondDMSFile(pos_in).positions)
        elif pos_in.split('.')[-1] == 'pdb':
            simulation.context.setPositions(PDBFile(pos_in).positions)
    
    else:
        simulation.context.setPositions(DesmondDMSFile(dms_in).positions)
    
    ### RUN SIMS
    print('-------------------------------')
    if soft:
        print('Soft interactions are turned on')
    else:
        print('Soft interactions are turned off')

    platform = simulation.context.getPlatform().getName()
    print("Platform: ", platform)
    print('E0: %.3e kJ/mol' %getEnergy(simulation))
    simulation.minimizeEnergy()
    print('E1: %.3e kJ/mol' %getEnergy(simulation))
    print('-------------------------------')

    ### SAVE THE LAST FRAME
    #pos = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)._value * 10
    pos = getPositions(simulation)

    if out.split('.')[-1] == 'dms':
        shutil.copy(dms_in, out)

        conn    = sqlite3.connect(out)
        cursor  = conn.cursor()
        n_atoms = cursor.execute('SELECT COUNT(DISTINCT id) FROM particle;').fetchall()[0][0]

        assert n_atoms == len(pos), 'n_atoms is different'

        for index in range(n_atoms):
            cursor.execute('UPDATE particle SET x = ?, y = ?, z = ? WHERE id = ?', 
                (*pos[index][0:3], index))

        conn.commit()
        conn.close()

    else:
        u = Universe(dms_in)
        u.atoms[['x','y','z']] = pos
        u.write(out)


def runMartiniNPT(dms_in, out, pos_in=None, soft=False, A=200, C=50,
    nonbondedCutoff=1.1, nonbondedMethod='CutoffPeriodic',
    nsteps=10000, dcdfreq=1000, csvfreq=1000, dt=0.02, P=1.0, T=310, semiisotropic=False,
    addForces=[], barfreq=100, frictionCoeff=5.0):
    '''
    CutoffNonPeriodic
    CutoffPeriodic
    NoCutoff
    '''
    system, dms = DMS2openmm(
            dms_in          = dms_in,
            nonbondedMethod = nonbondedMethod,
            nonbondedCutoff = nonbondedCutoff,
            soft            = soft,
            A               = A,
            C               = C).make()

    ### ADD BAROSTAT
    if semiisotropic:
        addForces.append(MembranePressure(P=P, T=T, barfreq=barfreq))
    else:
        addForces.append(Pressure(P=P, T=T, barfreq=barfreq))

    ### ADD additional forces
    for addf in addForces:
        system.addForce(copy.deepcopy(addf))
    
    ### PREPARE SIMS
    integrator = LangevinMiddleIntegrator(T*kelvin, frictionCoeff/picosecond, dt*picoseconds)
    simulation = Simulation(dms.topology, system, integrator)

    ### SET POSITIONS
    if pos_in:
        if pos_in.split('.')[-1] == 'dms':
            simulation.context.setPositions(DesmondDMSFile(pos_in).positions)
        elif pos_in.split('.')[-1] == 'pdb':
            simulation.context.setPositions(PDBFile(pos_in).positions)
    
    else:
        simulation.context.setPositions(DesmondDMSFile(dms_in).positions)

    ### RUN NPT
    prefix = '.'.join(out.split('.')[:-1])
    simulation.reporters.append(DCDReporter(      prefix + '.dcd', dcdfreq))
    simulation.reporters.append(StateDataReporter(prefix + '.csv', csvfreq, step=True, potentialEnergy=True, temperature=True))
    simulation.step(nsteps)
    
    cell = getCell(simulation)
    dimensions = [cell[0][0], cell[1][1], cell[2][2], 90, 90, 90]

    ### SAVE THE LAST FRAME
    #pos = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)._value * 10
    pos = getPositions(simulation)

    if out.split('.')[-1] == 'dms':
        shutil.copy(dms_in, out)

        conn    = sqlite3.connect(out)
        cursor  = conn.cursor()
        n_atoms = cursor.execute('SELECT COUNT(DISTINCT id) FROM particle;').fetchall()[0][0]

        assert n_atoms == len(pos), 'n_atoms is different'

        for index in range(n_atoms):
            cursor.execute('UPDATE particle SET x = ?, y = ?, z = ? WHERE id = ?', 
                (*pos[index][0:3], index))
        
        # global_cell id starts from 1 while particle id starts from 0 (why?)
        ids = sorted([tmp[0] for tmp in cursor.execute('SELECT id FROM global_cell;').fetchall()])
        assert len(set(ids)) == 3, f'global_cell has {len(set(ids))} entries (it should have been 3)'
        cursor.execute('UPDATE global_cell SET x=?, y=?, z=? WHERE id=?', (*cell[0], ids[0]))
        cursor.execute('UPDATE global_cell SET x=?, y=?, z=? WHERE id=?', (*cell[1], ids[1]))
        cursor.execute('UPDATE global_cell SET x=?, y=?, z=? WHERE id=?', (*cell[2], ids[2]))

        conn.commit()
        conn.close()

    else:
        u = Universe(dms_in)
        u.atoms[['x','y','z']] = pos
        u.cell = cell
        u.dimensions = dimensions
        u.write(out)



def runMartiniEMNPT(dms_in, out, pos_in=None, soft=False, A=200, C=50,
    nonbondedCutoff=1.1, nonbondedMethod='CutoffPeriodic',
    nsteps=10000, dcdfreq=1000, csvfreq=1000, dt=0.02, P=1.0, T=310, 
    semiisotropic=False, barfreq=100, addForces=[], frictionCoeff=5.0):
    '''
    CutoffNonPeriodic
    CutoffPeriodic
    NoCutoff
    '''
    system, dms = DMS2openmm(
            dms_in          = dms_in,
            nonbondedMethod = nonbondedMethod,
            nonbondedCutoff = nonbondedCutoff,
            soft            = soft,
            A               = A,
            C               = C).make()

    ### ADD BAROSTAT
    if semiisotropic:
        addForces.append(MembranePressure(P=P, T=T, barfreq=barfreq))
    else:
        addForces.append(Pressure(P=P, T=T, barfreq=barfreq))

    ### ADD additional forces
    for addf in addForces:
        system.addForce(copy.deepcopy(addf))
    
    ### PREPARE SIMS
    integrator = LangevinMiddleIntegrator(T*kelvin, frictionCoeff/picosecond, dt*picoseconds)
    simulation = Simulation(dms.topology, system, integrator)

    ### SET POSITIONS
    if pos_in:
        if pos_in.split('.')[-1] == 'dms':
            simulation.context.setPositions(DesmondDMSFile(pos_in).positions)
        elif pos_in.split('.')[-1] == 'pdb':
            simulation.context.setPositions(PDBFile(pos_in).positions)
    
    else:
        simulation.context.setPositions(DesmondDMSFile(dms_in).positions)

    
    ### RUN EM
    print('-------------------------------')
    if soft:
        print('Soft interactions are turned on')
    else:
        print('Soft interactions are turned off')

    platform = simulation.context.getPlatform().getName()
    print("Platform: ", platform)
    print('E0: %.3e kJ/mol' %getEnergy(simulation))
    simulation.minimizeEnergy()
    print('E1: %.3e kJ/mol' %getEnergy(simulation))
    print('-------------------------------')


    ### RUN NPT
    prefix = '.'.join(out.split('.')[:-1])
    simulation.reporters.append(DCDReporter(      prefix + '.dcd', dcdfreq))
    simulation.reporters.append(StateDataReporter(prefix + '.csv', csvfreq, step=True, potentialEnergy=True, temperature=True))
    simulation.step(nsteps)

    cell = getCell(simulation)
    dimensions = [cell[0][0], cell[1][1], cell[2][2], 90, 90, 90]

    ### SAVE THE LAST FRAME
    #pos = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)._value * 10
    pos = getPositions(simulation)

    if out.split('.')[-1] == 'dms':
        shutil.copy(dms_in, out)

        conn    = sqlite3.connect(out)
        cursor  = conn.cursor()
        n_atoms = cursor.execute('SELECT COUNT(DISTINCT id) FROM particle;').fetchall()[0][0]

        assert n_atoms == len(pos), 'n_atoms is different'

        for index in range(n_atoms):
            cursor.execute('UPDATE particle SET x = ?, y = ?, z = ? WHERE id = ?', 
                (*pos[index][0:3], index))
        
        # global_cell id starts from 1 while particle id starts from 0 (why?)
        ids = sorted([tmp[0] for tmp in cursor.execute('SELECT id FROM global_cell;').fetchall()])
        assert len(set(ids)) == 3, f'global_cell has {len(set(ids))} entries (it should have been 3)'
        cursor.execute('UPDATE global_cell SET x=?, y=?, z=? WHERE id=?', (*cell[0], ids[0]))
        cursor.execute('UPDATE global_cell SET x=?, y=?, z=? WHERE id=?', (*cell[1], ids[1]))
        cursor.execute('UPDATE global_cell SET x=?, y=?, z=? WHERE id=?', (*cell[2], ids[2]))

        conn.commit()
        conn.close()

    else:
        u = Universe(dms_in)
        u.atoms[['x','y','z']] = pos
        u.cell = cell
        u.dimensions = dimensions
        u.write(out)



# On Mac, all of the functions below work; On Linux, energy minimization fails.
# ctf = CustomTorsionForce('k * (acos(cos(theta-theta0)))^2')
# ctf = CustomTorsionForce('k * (cos(theta) - cos(theta0))^2')
# ctf = CustomTorsionForce('k * (theta-theta0)^2')
# Instead, using the below
# ctf = k*min(dtheta, 2*pi-dtheta)^2; dtheta = abs(theta-theta0); pi = 3.1415926535


def addPeptideTorsions(u, Kpeptide):
    ctf = CustomTorsionForce("k * (acos(0.999 * cos(theta-theta0)))^2")
    ctf.setName('PeptideTorsion')
    ctf.addPerTorsionParameter("k")
    ctf.addPerTorsionParameter("theta0")

    protein_resnames = three2one.keys()
    bA  = u.atoms.resname.isin(protein_resnames)
    bAN = u.atoms.name == 'N'
    bAH = (u.atoms.name.isin(['HN', 'H'])) | ((u.atoms.resname == 'PRO') & (u.atoms.name == 'CD'))
    bAC = u.atoms.name == 'C'
    bAO = u.atoms.name == 'O'
    bACA = u.atoms.name == 'CA'

    N = 0
    Natoms = u.atoms[bA & bAN]
    Hatoms = u.atoms[bA & bAH]
    Catoms = u.atoms[bA & bAC]
    Oatoms = u.atoms[bA & bAO]
    CAatoms = u.atoms[bA & bACA]

    for i in range(len(Catoms)):
        # Catom is a series object; other atoms are DataFrame objects
        Catom = Catoms.iloc[i]
        resid = Catom.resid
        chain = Catom.chain

        Oatom = Oatoms[ (Oatoms.resid == resid)   & (Oatoms.chain == chain)]
        Natom = Natoms[ (Natoms.resid == resid+1) & (Natoms.chain == chain)]
        Hatom = Hatoms[ (Hatoms.resid == resid+1) & (Hatoms.chain == chain)]
        CAatom0 = CAatoms[ (CAatoms.resid == resid) & (CAatoms.chain == chain) ]
        CAatom1 = CAatoms[ (CAatoms.resid == resid+1) & (CAatoms.chain == chain) ]
        
        if len(Oatom) *   len(Natom) * len(Hatom) != 1: continue
        if len(CAatom0) * len(Natom) * len(CAatom1) != 1: continue

        OatomIndex = Oatom.index.values[0]
        NatomIndex = Natom.index.values[0]
        HatomIndex = Hatom.index.values[0]
        CatomIndex = Catom.name
        # Catom['name'] -> CA
        # Catom.name    -> its index

        CAatom0Index = CAatom0.index.values[0]
        CAatom1Index = CAatom1.index.values[0]

        ctf.addTorsion(OatomIndex, CatomIndex, NatomIndex, HatomIndex, [Kpeptide, 3.141592])
        N += 1

        ctf.addTorsion(CAatom0Index, CatomIndex, NatomIndex, CAatom1Index, [Kpeptide, 3.141592])
        N += 1

    print(f'Adding PeptideTorsion for {N:d} isomers with K={Kpeptide:.2f}')
    return ctf


def addCisTransTorsions(u, Kcistrans, mapping, exclude=[], turn_off_torsion_warning=False):
    if not isinstance(exclude, list):
        exclude = [exclude]

    ctf = CustomTorsionForce("k * (acos(0.999 * cos(theta-theta0)))^2")
    ctf.setName('CisTransTorsion')
    ctf.addPerTorsionParameter("k")
    ctf.addPerTorsionParameter("theta0")

    N = 0
    for resname in set(u.atoms.resname):
        if resname in exclude: continue
        if resname not in mapping.RESI.keys():
            if resname == 'TIP3':
                continue
            elif resname.startswith('ROCK'):
                continue
            else:
                print(f'Warning: {resname} not in the mapping scheme - skipping CisTransTorsions for this residue')
                continue
        for isomer in ['cis', 'trans']:
            atomset = mapping.RESI[resname][isomer]
            for atoms in atomset:
                bA  = u.atoms.resname == resname
                bA0 = u.atoms.name    == atoms[0]
                bA1 = u.atoms.name    == atoms[1]
                bA2 = u.atoms.name    == atoms[2]
                bA3 = u.atoms.name    == atoms[3]

                atomA = u.atoms[bA & bA0].index
                atomB = u.atoms[bA & bA1].index
                atomC = u.atoms[bA & bA2].index
                atomD = u.atoms[bA & bA3].index
                
                checkbA = not (len(atomA) == len(atomB) == len(atomC) == len(atomD))
                if checkbA and turn_off_torsion_warning:
                    continue
                
                assert len(atomA) == len(atomB) == \
                    len(atomC) == len(atomD), \
                    "the length of atoms for cistrans torsions is different"

                for a, b, c, d in zip(atomA, atomB, atomC, atomD):
                    N += 1
                    if isomer == 'cis':
                        ctf.addTorsion(a, b, c, d, [Kcistrans, 0.000])

                    elif isomer == 'trans':
                        ctf.addTorsion(a, b, c, d, [Kcistrans, 3.141592])

    print(f'Adding CisTransTorsion for {N:d} isomers with K={Kcistrans:.2f}')
    return ctf


def addDihedralTorsions(u, Kdihedral, mapping, exclude=[], turn_off_torsion_warning=False):
    if not isinstance(exclude, list):
        exclude = [exclude]

    ctf = CustomTorsionForce("k * (acos(0.999 * cos(theta-theta0)))^2")
    ctf.setName('DihedralTorsion')
    ctf.addPerTorsionParameter("k")
    ctf.addPerTorsionParameter("theta0")


    N = 0
    for resname in set(u.atoms.resname):
        if resname in exclude: continue
        if resname not in mapping.RESI.keys():
            if resname == 'TIP3':
                continue
            elif resname.startswith('ROCK'):
                continue
            else:
                print(f'Warning: {resname} not in the mapping scheme - skipping DihedralTorsions for this residue')
                continue

        dihedrals = mapping.RESI[resname]['dihedral']
        for dihedral in dihedrals:
            for ii in range(len(dihedral) // 5):
                dihe = dihedral[ii * 5 : (ii+1) * 5]

                bA  = u.atoms.resname == resname
                bA0 = u.atoms.name    == dihe[0]
                bA1 = u.atoms.name    == dihe[1]
                bA2 = u.atoms.name    == dihe[2]
                bA3 = u.atoms.name    == dihe[3]

                atomA = u.atoms[bA & bA0].index
                atomB = u.atoms[bA & bA1].index
                atomC = u.atoms[bA & bA2].index
                atomD = u.atoms[bA & bA3].index

                checkbA = not (len(atomA) == len(atomB) == len(atomC) == len(atomD))
                if checkbA and turn_off_torsion_warning: continue

                assert len(atomA) == len(atomB) == \
                    len(atomC) == len(atomD), \
                    "the length of atoms for dihedral torsions is different"

                for a, b, c, d in zip(atomA, atomB, atomC, atomD):
                    N += 1
                    ctf.addTorsion(a, b, c, d, [Kdihedral, float(dihe[4]) * 3.141592 / 180])

    print(f'Adding DihedralTorsion for {N:d} isomers with K={Kdihedral:.2f}')
    return ctf


def addAntiDihedralTorsions(u, Kdihedral, mapping, exclude=[], turn_off_torsion_warning=False):
    if not isinstance(exclude, list):
        exclude = [exclude]

    ctf = CustomTorsionForce("k * (acos(0.999 * cos(theta-theta0)))^2")
    ctf.setName('AntiDihedralTorsion')
    ctf.addPerTorsionParameter("k")
    ctf.addPerTorsionParameter("theta0")


    N = 0
    for resname in set(u.atoms.resname):
        if resname in exclude: continue
        if resname not in mapping.RESI.keys():
            if resname == 'TIP3':
                continue
            elif resname.startswith('ROCK'):
                continue
            else:
                print(f'Warning: {resname} not in the mapping scheme - skipping AntiDihedralTorsions for this residue')
                continue

        dihedrals = mapping.RESI[resname]['antidihedral']
        for dihedral in dihedrals:
            for ii in range(len(dihedral) // 5):
                dihe = dihedral[ii * 5 : (ii+1) * 5]

                bA  = u.atoms.resname == resname
                bA0 = u.atoms.name    == dihe[0]
                bA1 = u.atoms.name    == dihe[1]
                bA2 = u.atoms.name    == dihe[2]
                bA3 = u.atoms.name    == dihe[3]

                atomA = u.atoms[bA & bA0].index
                atomB = u.atoms[bA & bA1].index
                atomC = u.atoms[bA & bA2].index
                atomD = u.atoms[bA & bA3].index

                checkbA = not (len(atomA) == len(atomB) == len(atomC) == len(atomD))
                if checkbA and turn_off_torsion_warning: continue

                assert len(atomA) == len(atomB) == \
                    len(atomC) == len(atomD), \
                    "the length of atoms for dihedral torsions is different"

                for a, b, c, d in zip(atomA, atomB, atomC, atomD):
                    N += 1
                    ctf.addTorsion(a, b, c, d, [-Kdihedral, float(dihe[4]) * 3.141592 / 180])

    print(f'Adding AntiDihedralTorsion for {N:d} isomers with K={Kdihedral:.2f}')
    return ctf


def addChiralTorsions(u, Kchiral, mapping, exclude=[], turn_off_torsion_warning=False):
    if not isinstance(exclude, list):
        exclude = [exclude]

    ctf = CustomTorsionForce('k * (acos(0.999 * cos(theta-theta0)))^2')
    ctf.setName('ChiralTorsion')
    ctf.addPerTorsionParameter("k")
    ctf.addPerTorsionParameter("theta0")

    N = 0
    for resname in set(u.atoms.resname):
        if resname in exclude: continue
        if resname not in mapping.RESI.keys():
            if resname == 'TIP3':
                continue
            elif resname.startswith('ROCK'):
                continue
            else:
                print(f'Warning: {resname} not in the mapping scheme - skipping ChiralTorsions for this residue')
                continue
        chirals = mapping.RESI[resname]['chiral']
        for chiral in chirals:
            bA  = u.atoms.resname == resname
            bA0 = u.atoms.name    == chiral[0]
            bA1 = u.atoms.name    == chiral[1]
            bA2 = u.atoms.name    == chiral[2]
            bA3 = u.atoms.name    == chiral[3]              
            bA4 = u.atoms.name    == chiral[4]

            target = u.atoms[bA & bA0].index
            center = u.atoms[bA & bA1].index
            atomC  = u.atoms[bA & bA2].index
            atomD  = u.atoms[bA & bA3].index
            atomE  = u.atoms[bA & bA4].index
            
            checkbA = not (len(target) == len(center) == len(atomC) == len(atomD) == len(atomE))
            if checkbA and turn_off_torsion_warning: continue
                    
            assert len(target) == len(center) == \
                len(atomC) == len(atomD) == len(atomE), \
                "the length of atoms for chiral torsions is different"

            for a, b, c, d, e in zip(target, center, atomC, atomD, atomE):
                ctf.addTorsion(b, c, d, e, [Kchiral, -35.26 * 3.141592/180])
                ctf.addTorsion(b, d, e, c, [Kchiral, -35.26 * 3.141592/180])
                ctf.addTorsion(b, e, c, d, [Kchiral, -35.26 * 3.141592/180])
    
                ctf.addTorsion(a, c, d, e, [Kchiral, -70.53 * 3.141592/180])
                ctf.addTorsion(a, d, e, c, [Kchiral, -70.53 * 3.141592/180])
                ctf.addTorsion(a, e, c, d, [Kchiral, -70.53 * 3.141592/180])
                N += 1

    print('Adding ChiralTorsion for {:d} chirals'.format(N))
    return ctf


def addPosre(u, bfactor_posre, fcx, fcy, fcz):
    '''Apply positional restraints on atoms whose bfactors are larger than bfactor_posre'''
    cef = CustomExternalForce("hkx*(x-x0)^2+hky*(y-y0)^2+hkz*(z-z0)^2")
    cef.addPerParticleParameter("x0")
    cef.addPerParticleParameter("y0")
    cef.addPerParticleParameter("z0")
    cef.addPerParticleParameter("hkx")
    cef.addPerParticleParameter("hky")
    cef.addPerParticleParameter("hkz")

    if 'bfactor' not in u.atoms:
        print('bfactor not in atoms; skipping Positional Restraints')
        return cef

    bA = u.atoms.bfactor > bfactor_posre
    df = u.atoms[bA]

    for index, row in df.iterrows():
        x0d = (row.x * angstrom).value_in_unit(nanometer)
        y0d = (row.y * angstrom).value_in_unit(nanometer)
        z0d = (row.z * angstrom).value_in_unit(nanometer)
        hfcxd = fcx*kilojoule_per_mole/nanometer**2
        hfcyd = fcy*kilojoule_per_mole/nanometer**2
        hfczd = fcz*kilojoule_per_mole/nanometer**2
        cef.addParticle(index,[ x0d, y0d, z0d, hfcxd,  hfcyd,  hfczd])

    print(f'Adding Posre of ({fcx:.1f}, {fcy:.1f}, {fcz:.1f}) kJ/mol/nm^2 for {len(df):d} atoms whose bfactor > {bfactor_posre:.2f})')
    return cef

def addPosrePeriodic(u, bfactor_posre, k):
    '''Apply positional restraints on atoms whose bfactors are larger than bfactor_posre'''
    cef = CustomExternalForce("k * periodicdistance(x, y, z, x0, y0, z0)^2")
    cef.addPerParticleParameter("x0")
    cef.addPerParticleParameter("y0")
    cef.addPerParticleParameter("z0")
    cef.addPerParticleParameter("k")

    if 'bfactor' not in u.atoms:
        print('bfactor not in atoms; skipping Positional Restraints')
        return cef

    bA = u.atoms.bfactor > bfactor_posre
    df = u.atoms[bA]

    for index, row in df.iterrows():
        x0d = (row.x * angstrom).value_in_unit(nanometer)
        y0d = (row.y * angstrom).value_in_unit(nanometer)
        z0d = (row.z * angstrom).value_in_unit(nanometer)
        fc  = k * kilojoule_per_mole/nanometer**2
        cef.addParticle(index,[ x0d, y0d, z0d, fc])

    print('Adding Periodic Posre for {:d} atoms whose bfactor > {:.2f})'.format(len(df), bfactor_posre))
    return cef


def addPosrePeriodicZ(u, bfactor_posre, k):
    '''Apply positional restraints on atoms whose bfactors are larger than bfactor_posre'''
    cef = CustomExternalForce("k * periodicdistance(x, y, z, x, y, z0)^2")
    cef.addPerParticleParameter("z0")
    cef.addPerParticleParameter("k")

    if 'bfactor' not in u.atoms:
        print('bfactor not in atoms; skipping Positional Restraints')
        return cef

    bA = u.atoms.bfactor > bfactor_posre
    df = u.atoms[bA]

    for index, row in df.iterrows():
        z0d = (row.z * angstrom).value_in_unit(nanometer)
        fc  = k * kilojoule_per_mole/nanometer**2
        cef.addParticle(index,[ z0d, fc])

    print('Adding Periodic Posre Z for {:d} atoms whose bfactor > {:.2f})'.format(len(df), bfactor_posre))
    return cef


def addRefPosre(u, refstructure, fcx, fcy, fcz):
    ref = Universe(refstructure)

    cef = CustomExternalForce("hkx*(x-x0)^2+hky*(y-y0)^2+hkz*(z-z0)^2")
    cef.addPerParticleParameter("x0")
    cef.addPerParticleParameter("y0")
    cef.addPerParticleParameter("z0")
    cef.addPerParticleParameter("hkx")
    cef.addPerParticleParameter("hky")
    cef.addPerParticleParameter("hkz")

    N = 0
    for index, atom in u.atoms.iterrows():
        bA1 = atom['name']  == ref.atoms.name
        bA2 = atom['resid'] == ref.atoms.resid
        bA3 = atom['chain'] == ref.atoms.chain

        refatoms = ref.atoms[bA1 & bA2 & bA3]
        if len(refatoms) == 0:
            continue
        elif len(refatoms) == 1:
            refatom = refatoms.iloc[0]
            #print(refatom.resid, refatom.chain, refatom['name'], atom.x, refatom.x, atom.y, refatom.y, atom.z, refatom.z)
            x0d = (refatom.x * angstrom).value_in_unit(nanometer)
            y0d = (refatom.y * angstrom).value_in_unit(nanometer)
            z0d = (refatom.z * angstrom).value_in_unit(nanometer)
            hfcxd = fcx*kilojoule_per_mole/nanometer**2
            hfcyd = fcy*kilojoule_per_mole/nanometer**2
            hfczd = fcz*kilojoule_per_mole/nanometer**2
            cef.addParticle(index,[ x0d, y0d, z0d, hfcxd,  hfcyd,  hfczd])
            N += 1

        else:
            assert 0 == 1, '/{:s}:{:d}@{:s} '.format(atom['chain'], atom['resid'], atom['name']) + \
            'more than one atom with the same name, resid, chain?'

    print('Adding RefPosre for {:d} atoms that exist in {:s}'.format(N, refstructure))
    return cef


def addRefPosrePeriodic(u, refstructure, k):
    ref = Universe(refstructure)

    cef = CustomExternalForce("k * periodicdistance(x, y, z, x0, y0, z0)^2")
    cef.addPerParticleParameter("x0")
    cef.addPerParticleParameter("y0")
    cef.addPerParticleParameter("z0")
    cef.addPerParticleParameter("k")

    N = 0
    for index, atom in u.atoms.iterrows():
        bA1 = atom['name']  == ref.atoms.name
        bA2 = atom['resid'] == ref.atoms.resid
        bA3 = atom['chain'] == ref.atoms.chain

        refatoms = ref.atoms[bA1 & bA2 & bA3]
        if len(refatoms) == 0:
            continue
        elif len(refatoms) == 1:
            refatom = refatoms.iloc[0]
            #print(refatom.resid, refatom.chain, refatom['name'], atom.x, refatom.x, atom.y, refatom.y, atom.z, refatom.z)
            x0d = (refatom.x * angstrom).value_in_unit(nanometer)
            y0d = (refatom.y * angstrom).value_in_unit(nanometer)
            z0d = (refatom.z * angstrom).value_in_unit(nanometer)
            fc  = k * kilojoule_per_mole/nanometer**2
            cef.addParticle(index,[ x0d, y0d, z0d, k])
            N += 1

        else:
            assert 0 == 1, '/{:s}:{:d}@{:s} '.format(atom['chain'], atom['resid'], atom['name']) + \
            'more than one atom with the same name, resid, chain?'

    print('Adding RefPosre for {:d} atoms that exist in {:s}'.format(N, refstructure))
    return cef

def addBonds(u, xml, pdb=None):
    print("Adding bonds for non-protein residues - started")
    bond_records = []

    # add bonds except for protein residues
    protein_resnames = three2one.keys()
    for resname in set(u.atoms.resname):
        # openMM PDBFile takes care of protein bonds
        if resname in protein_resnames: continue

        if resname not in xml.RESI.keys():
            print(f'Warning: openMM xml does not have {resname}')
            continue

        bonds = xml.RESI[resname]['bonds']

        for bond in bonds:
            bA  = u.atoms.resname == resname
            bA0 = u.atoms.name == bond[0]
            bA1 = u.atoms.name == bond[1]

            atomA = u.atoms[bA & bA0].index
            atomB = u.atoms[bA & bA1].index

            assert len(atomA) == len(atomB), \
                    f"{resname}: {bond[0]} and {bond[1]} - the length of atoms for bonds is different"

            for a, b in zip(atomA, atomB):
                bond_records.append([a, b])


    print("Adding bonds for non-protein residues - finished")
    if pdb:
        pdbatoms = [atom for atom in pdb.topology.atoms()]
        for bond in bond_records:
            a = bond[0]
            b = bond[1]
            pdb.topology.addBond(pdbatoms[a], pdbatoms[b])

        return bond_records, pdb

    else:
        return bond_records


def getTopPDB(structure, ff=[], ff_add=[]):
    xml = ReadXML(ff=ff, ff_add=ff_add)
    u   = Universe(structure)
    if structure.split('.')[-1] == 'pdb':
        pdb = PDBFile(structure)
    elif structure.split('.')[-1] == 'dms':
        pdb = DesmondDMSFile(structure)
    pdbatoms = [atom for atom in pdb.topology.atoms()]

    print("Adding bonds for non-protein residues - started")
    # add bonds except for protein residues
    protein_resnames = three2one.keys()
    for resname in set(u.atoms.resname):
        # openMM PDBFile takes care of protein bonds
        if resname in protein_resnames: continue

        if resname not in xml.RESI.keys():
            print(f'Warning: openMM xml does not have {resname}')
            continue

        bonds = xml.RESI[resname]['bonds']

        for bond in bonds:
            bA  = u.atoms.resname == resname
            bA0 = u.atoms.name == bond[0]
            bA1 = u.atoms.name == bond[1]

            atomA = u.atoms[bA & bA0].index
            atomB = u.atoms[bA & bA1].index

            assert len(atomA) == len(atomB), "the length of atoms for bonds is different"

            for a, b in zip(atomA, atomB):
                pdb.topology.addBond(pdbatoms[a], pdbatoms[b])

    print("Adding bonds for non-protein residues - finished")
    return pdb


def getBonds(structure, ff=[], ff_add=[]):
    xml  = ReadXML(ff=ff, ff_add=ff_add)
    data = []
    
    if isinstance(structure, str):
        u = Universe(structure)
    else:
        u = structure

    print("Adding bonds - started")
    for resname in set(u.atoms.resname):
        if resname not in xml.RESI.keys():
            print(f'Warning: openMM xml does not have {resname}')
            continue

        bonds = xml.RESI[resname]['bonds']

        for bond in bonds:
            bA0 = u.atoms.resname == resname
            bA1 = u.atoms.name.isin(bond)
            new = u.atoms[bA0 & bA1]

            for i in range(len(new.index) - 1):
                if new.iloc[i].resn == new.iloc[i+1].resn:
                    data.append([new.iloc[i].id, new.iloc[i+1].id])

            #bA0 = u.atoms.name == bond[0]
            #bA1 = u.atoms.name == bond[1]

            #atomA = u.atoms[bA & bA0].index
            #atomB = u.atoms[bA & bA1].index
            #
            #assert len(atomA) == len(atomB), f"bond: /{resname} @{bond[0]} @{bond[1]}"

            #for a, b in zip(atomA, atomB):
            #    data.append([a, b])
    

    ### PROTEIN BACKBONE
    protein_bA     = u.atoms.resname.isin(three2one.keys())
    protein_chains = u.atoms[protein_bA].chain.unique()

    for chain in protein_chains:
        chain_bA    = u.atoms.chain == chain
        chain_atoms = u.atoms[protein_bA & chain_bA]
        resid_min   = chain_atoms.resid.min()
        resid_max   = chain_atoms.resid.max()
        
        ### C : +N bonds
        bAN = chain_atoms['name'] == 'N'
        bAC = chain_atoms['name'] == 'C'
        assert np.sum(bAN) == np.sum(bAC), 'protein: len(N) != len(C)'
            
        tmp = np.array([chain_atoms[bAC].index[:-1], chain_atoms[bAN].index[1:]], dtype=np.int64).T
        data.extend(list(tmp))
        

        ### 0N - HT1 HT2 HT3
        residue_atoms = chain_atoms[chain_atoms.resid == resid_min]
        bAN = residue_atoms['name'] == 'N'
        bAH = residue_atoms['name'].isin(['HT1', 'HT2', 'HT3'])
        assert np.sum(bAN) == 1, f'{np.sum(bAN)} N termini in chain {chain}?'

        for i in residue_atoms[bAH].index:
            data.append([residue_atoms[bAN].index[0], i])

        ### -1C - OT1 OT2
        residue_atoms = chain_atoms[chain_atoms.resid == resid_max]
        bAC = residue_atoms['name'] == 'C'
        bAO = residue_atoms['name'].isin(['OT1', 'OT2', 'OXT1', 'OXT2', 'OXT'])
        assert np.sum(bAC) == 1, f'{np.sum(bAC)} C termini in chain {chain}?'

        for i in residue_atoms[bAO].index:
            data.append([residue_atoms[bAC].index[0], i])

    print("Adding bonds - finished")
    return data

def addFlatBottomZ(u, bfactor_posre, radius, rfb, R0z=0.0, fc=1000.0, chain=False):
    '''Apply Flat-bottomed position restraints for sphere simulations

    | d - radius | < rfb  ---> pot = 0
      d - radius   < -rfb ---> pot = (radius - rfb - x)**2
      d - radius   > +rfb ---> pot = (radius + rfb - x)**2
    
    Parameters
    ----------
    u : Universe
    bfactor_posre : float
        atoms whose bfactors are larger than bfactor_posre will have this restraints
    radius : float
        radius in unit of Angstrom
    rfb : float
        distance from the center with a flat potential in unit of Angstrom
    fc : float
        force constant in unit of kilojoule_per_mole/naometer**2
    R0 : list or array
        a reference position with a shape of ``(3,)``.
    '''

    pot  = "fc * (radius + rfb - d)^2 * step( d - radius - rfb) + "
    pot += "fc * (radius - rfb - d)^2 * step(-d + radius - rfb); "
    pot += "d =  ((z-R0z)^2)^0.5;"
    
    cef = CustomExternalForce(pot)
    cef.addGlobalParameter("fc",  fc * kilojoule_per_mole/nanometer**2)
    cef.addGlobalParameter("R0z",    R0z    * 0.1 * nanometer)
    cef.addGlobalParameter("rfb",    rfb    * 0.1 * nanometer)
    cef.addPerParticleParameter("radius")

    bA1 = u.atoms.bfactor > bfactor_posre
    if chain:
        bA2 = u.atoms.chain == chain
        df  = u.atoms[bA1 & bA2]
    else:
        df = u.atoms[bA1]

    for index, row in df.iterrows():
        cef.addParticle(index, [radius * 0.1 * nanometer])

    if chain:
        print(f'Adding FlatBottomZ for {len(df):d} atoms whose bfactor > {bfactor_posre:.2f} and chain is {chain})')
    else:
        print(f'Adding FlatBottomZ for {len(df):d} atoms whose bfactor > {bfactor_posre:.2f}')

    return cef



def addFlatBottomSphere(u, bfactor_posre, radius, rfb, R0=[0,0,0], fc=1000.0):
    '''Apply Flat-bottomed position restraints for sphere simulations

    | d - radius | < rfb  ---> pot = 0
      d - radius   < -rfb ---> pot = (radius - rfb - x)**2
      d - radius   > +rfb ---> pot = (radius + rfb - x)**2
    
    Parameters
    ----------
    u : Universe
    bfactor_posre : float
        atoms whose bfactors are larger than bfactor_posre will have this restraints
    radius : float
        radius in unit of Angstrom
    rfb : float
        distance from the center with a flat potential in unit of Angstrom
    fc : float
        force constant in unit of kilojoule_per_mole/naometer**2
    R0 : list or array
        a reference position with a shape of ``(3,)``.
    '''

    pot  = "fc * (radius + rfb - d)^2 * step(+d - radius - rfb) + "
    pot += "fc * (radius - rfb - d)^2 * step(-d + radius - rfb); "
    pot += "d = (  (x-R0x)^2 + (y-R0y)^2 + (z-R0z)^2  )^0.5;"
    
    cef = CustomExternalForce(pot)
    cef.addGlobalParameter("fc",  fc * kilojoule_per_mole/nanometer**2)
    cef.addGlobalParameter("R0x",    R0[0]  * 0.1 * nanometer)
    cef.addGlobalParameter("R0y",    R0[1]  * 0.1 * nanometer)
    cef.addGlobalParameter("R0z",    R0[2]  * 0.1 * nanometer)
    cef.addGlobalParameter("rfb",    rfb    * 0.1 * nanometer)
    cef.addGlobalParameter("radius", radius * 0.1 * nanometer)

    bA = u.atoms.bfactor > bfactor_posre
    df = u.atoms[bA]

    for index, row in df.iterrows():
        cef.addParticle(index)

    print('Adding FlatBottomSphere for {:d} atoms whose bfactor > {:.2f})'.format(len(df), bfactor_posre))
    return cef


def verifyFlatBottomSphere(radius=300, rfb=50, R0=[0,0,0], fc=1.0):
    def step(x):
        return 0 if x < 0 else 1
    
    def pot1(d):
        if d < radius - rfb:
            return (d - (radius - rfb))**2
        elif d > radius + rfb:
            return (d - (radius + rfb))**2
        else:
            return 0
    
    def pot2(d):
        term1 = fc * (radius + rfb - d) ** 2 * step(+d - radius - rfb)
        term2 = fc * (radius - rfb - d) ** 2 * step(-d + radius - rfb)
        return term1 + term2
        
    ds = np.linspace(0, 400, 100); y1 = []; y2 = [];
    for d in ds:
        y1.append(pot1(d))
        y2.append(pot2(d))
    
    print(y1)
    print(y2)


def addSpherePosre(u, bfactor_posre, radius, rfb, R0=[0,0,0], fc=100.0, chain=None):
    '''Apply Flat-bottomed position restraints for sphere simulations

    | d - radius | < rfb  ---> pot = 0
      d - radius   < -rfb ---> pot = (radius - rfb - x)**2
      d - radius   > +rfb ---> pot = (radius + rfb - x)**2
    
    Parameters
    ----------
    u : Universe
    bfactor_posre : float
        atoms whose bfactors are larger than bfactor_posre will have this restraints
    radius : float
        radius in unit of Angstrom
    rfb : float
        distance from the center with a flat potential in unit of Angstrom
    fc : float
        force constant in unit of kilojoule_per_mole/naometer**2
    R0 : list or array
        a reference position with a shape of ``(3,)``.
    '''

    pot  = "fc * (radius + rfb - d)^2 * step( d - radius - rfb) + "
    pot += "fc * (radius - rfb - d)^2 * step(-d + radius - rfb); "
    pot += "d = (  (x-R0x)^2 + (y-R0y)^2 + (z-R0z)^2  )^0.5;"
    
    cef = CustomExternalForce(pot)
    cef.addGlobalParameter("fc",     fc * kilojoule_per_mole/nanometer**2)
    cef.addGlobalParameter("R0x",    R0[0]  * 0.1 * nanometer)
    cef.addGlobalParameter("R0y",    R0[1]  * 0.1 * nanometer)
    cef.addGlobalParameter("R0z",    R0[2]  * 0.1 * nanometer)
    cef.addGlobalParameter("rfb",    rfb    * 0.1 * nanometer)
    cef.addPerParticleParameter("radius")

    bA1 = u.atoms.bfactor > bfactor_posre
    if chain:
        bA2 = u.atoms.chain == chain
        df  = u.atoms[bA1 & bA2]
    else:
        df = u.atoms[bA1]

    for index, row in df.iterrows():
        cef.addParticle(index, [radius * 0.1 * nanometer])

    if chain:
        print(f'Adding Sphere Posre for {len(df):d} atoms whose bfactor > {bfactor_posre:.2f} and chain is {chain})')
    else:
        print(f'Adding Sphere Posre for {len(df):d} atoms whose bfactor > {bfactor_posre:.2f}')
    return cef


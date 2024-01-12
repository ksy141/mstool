from openmm.unit import *
from openmm.app  import *
from openmm      import *

def getEnergy(simulation):
    return simulation.context.getState(getEnergy=True).getPotentialEnergy()._value

def Pressure(system, P=1.0, T=310.0):
    system.addForce(MonteCarloBarostat(P*bar, T*kelvin))
    return system

def MembranePressure(system, P=1.0, T=310.0, r=0.0):
    system.addForce(MonteCarloMembraneBarostat(P*bar, r*bar*nanometer, T*kelvin,
            MonteCarloMembraneBarostat.XYIsotropic, MonteCarloMembraneBarostat.ZFree))
    return system

def getDMSEnergy(dms_in, nonbondedCutoff=1.1, nonbondedMethod='CutoffPeriodic'):
    from ..core.dms2openmm import DMS2openmm

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


def runMartiniEM(dms_in, out, pos_in=None, soft=False, A=200, C=50,
    nonbondedCutoff=1.1, nonbondedMethod='CutoffPeriodic', T=310, dt=0.002,
    addForces=[]):
    '''
    CutoffNonPeriodic
    CutoffPeriodic
    NoCutoff
    '''

    from ..core.dms2openmm import DMS2openmm
    from ..core.universe   import Universe
    import sqlite3

    system, dms = DMS2openmm(
            dms_in          = dms_in,
            nonbondedMethod = nonbondedMethod,
            nonbondedCutoff = nonbondedCutoff,
            soft            = soft,
            A               = A,
            C               = C).make()
    
    ### ADD additional forces
    for Force in addForces:
        system.addForce(Force)

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
    
    print('E0: %.3e kJ/mol' %getEnergy(simulation))
    simulation.minimizeEnergy()
    print('E1: %.3e kJ/mol' %getEnergy(simulation))
    print('-------------------------------')

    ### SAVE THE LAST FRAME
    pos = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)._value * 10

    if out.split('.')[-1] == 'dms':
        import shutil
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
    addForces=[]):
    '''
    CutoffNonPeriodic
    CutoffPeriodic
    NoCutoff
    '''
    
    from ..core.dms2openmm import DMS2openmm
    from ..core.universe   import Universe

    system, dms = DMS2openmm(
            dms_in          = dms_in,
            nonbondedMethod = nonbondedMethod,
            nonbondedCutoff = nonbondedCutoff,
            soft            = soft,
            A               = A,
            C               = C).make()

    ### ADD BAROSTAT
    if semiisotropic:
        system = MembranePressure(system, P=P, T=T)
    else:
        system = Pressure(system, P=P, T=T)

    ### ADD additional forces
    for Force in addForces:
        system.addForce(Force)
    
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
    prefix = '.'.join(out.split('.')[:-1])
    simulation.reporters.append(DCDReporter(      prefix + '.dcd', dcdfreq))
    simulation.reporters.append(StateDataReporter(prefix + '.csv', csvfreq, step=True, potentialEnergy=True, temperature=True))
    simulation.step(nsteps)

    ### SAVE THE LAST FRAME
    u = Universe(dms_in)
    numpypositions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True) * 10
    u.setpositions(numpypositions)
    u.write(out)



# On Mac, all of the functions below work; On Linux, energy minimization fails.
# ctf = CustomTorsionForce('k * (acos(cos(theta-theta0)))^2')
# ctf = CustomTorsionForce('k * (cos(theta) - cos(theta0))^2')
# ctf = CustomTorsionForce('k * (theta-theta0)^2')
# Instead, using the below
# ctf = k*min(dtheta, 2*pi-dtheta)^2; dtheta = abs(theta-theta0); pi = 3.1415926535


def addPeptideTorsions(u, Kpeptide):
    ctf = CustomTorsionForce("k*min(dtheta, 2*pi-dtheta)^2; dtheta = abs(theta-theta0); pi = 3.1415926535")
    ctf.setName('PeptideTorsion')
    ctf.addPerTorsionParameter("k")
    ctf.addPerTorsionParameter("theta0")

    protein_resnames = three2one.keys()
    bA  = u.atoms.resname.isin(protein_resnames)
    bAN = u.atoms.name == 'N'
    bAH = (u.atoms.name.isin(['HN', 'H'])) | ((u.atoms.resname == 'PRO') & (u.atoms.name == 'CD'))
    bAC = u.atoms.name == 'C'
    bAO = u.atoms.name == 'O'

    N = 0
    Natoms = u.atoms[bA & bAN]
    Hatoms = u.atoms[bA & bAH]
    Catoms = u.atoms[bA & bAC]
    Oatoms = u.atoms[bA & bAO]

    for i in range(len(Catoms)):
        # Catom is a series object; other atoms are DataFrame objects
        Catom = Catoms.iloc[i]
        resid = Catom.resid
        chain = Catom.chain

        Oatom = Oatoms[ (Oatoms.resid == resid)   & (Oatoms.chain == chain)]
        Natom = Natoms[ (Natoms.resid == resid+1) & (Natoms.chain == chain)]
        Hatom = Hatoms[ (Hatoms.resid == resid+1) & (Hatoms.chain == chain)]
        if len(Oatom) * len(Natom) * len(Hatom) != 1: continue

        OatomIndex = Oatom.index.values[0]
        NatomIndex = Natom.index.values[0]
        HatomIndex = Hatom.index.values[0]
        CatomIndex = Catom.name
        # Catom['name'] -> CA
        # Catom.name    -> its index
        ctf.addTorsion(OatomIndex, CatomIndex, NatomIndex, HatomIndex, [Kpeptide, 3.141592])
        N += 1

    print('Adding PeptideTorsion for {:d} isomers'.format(N))
    return ctf


def addCisTransTorsions(u, Kcistrans, mapping, exclude=[]):
    if not isinstance(exclude, list):
        exclude = [exclude]

    ctf = CustomTorsionForce("k*min(dtheta, 2*pi-dtheta)^2; dtheta = abs(theta-theta0); pi = 3.1415926535")
    ctf.setName('CisTransTorsion')
    ctf.addPerTorsionParameter("k")
    ctf.addPerTorsionParameter("theta0")


    N = 0
    for resname in set(u.atoms.resname):
        if resname in exclude: continue
        if resname not in mapping.RESI.keys():
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

                assert len(atomA) == len(atomB) == \
                    len(atomC) == len(atomD), \
                    "the length of atoms for cistrans torsions is different"

                for a, b, c, d in zip(atomA, atomB, atomC, atomD):
                    N += 1
                    if isomer == 'cis':
                        ctf.addTorsion(a, b, c, d, [Kcistrans, 0.000])

                    elif isomer == 'trans':
                        ctf.addTorsion(a, b, c, d, [Kcistrans, 3.141592])

    print('Adding CisTransTorsion for {:d} isomers'.format(N))
    return ctf


def addChiralTorsions(u, Kchiral, mapping, exclude=[]):
    if not isinstance(exclude, list):
        exclude = [exclude]

    ctf = CustomTorsionForce('k * (theta-theta0)^2')
    ctf.setName('ChiralTorsion')
    ctf.addPerTorsionParameter("k")
    ctf.addPerTorsionParameter("theta0")

    N = 0
    for resname in set(u.atoms.resname):
        if resname in exclude: continue
        if resname not in mapping.RESI.keys():
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

    print('Adding Posre for {:d} atoms whose bfactor > {:.2f})'.format(len(df), bfactor_posre))
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

            assert len(atomA) == len(atomB), "the length of atoms for bonds is different"

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

    print('Adding Posre for {:d} atoms whose bfactor > {:.2f})'.format(len(df), bfactor_posre))
    return cef



from   openmm        import *
from   openmm.app    import *
from   openmm.unit   import *

from   .universe     import Universe
from   .readmappings import ReadMappings
from   .readxml      import ReadXML

from   ..utils.protein_sel import three2one
from   ..utils.rock        import Rock
from   ..utils.openmmutils import *

import numpy  as np
import pandas as pd
import sys
sys.setrecursionlimit(1000000)

# On Mac, all of the functions below work; On Linux, energy minimization fails.
# ctf = CustomTorsionForce('k * (acos(cos(theta-theta0)))^2')
# ctf = CustomTorsionForce('k * (cos(theta) - cos(theta0))^2')
# ctf = CustomTorsionForce('k * (theta-theta0)^2')


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






class REM:
    # pdb1 = PDBFile('villin.pdb')
    # pdb2 = PDBFile('ala_ala_ala.pdb')
    # modeller = Modeller(pdb1.topology, pdb1.positions)
    # modeller.add(pdb2.topology, pdb2.positions)
    # mergedTopology = modeller.topology
    # mergedPositions = modeller.positions

    def __init__(self, structure=None, out=None, protein=None, refposre=None,
        rock=None, rockrcut=1.2, rockKbond=5000.0, rockresname='ROCK',
        rcut=1.2, pbc=True, 
        A=100, C=50,
        mapping = [], mapping_add = [], 
        ff      = [], ff_add = [],
        Kchiral=300, Kpeptide=300, Kcistrans=300,
        fcx = 1000.0, fcy = 1000.0, fcz = 1000.0,
        bfactor_posre = 0.5, add_bonds=True):


        ### NonbondedMethod
        if pbc:
            self.nonbondedMethod = CutoffPeriodic
        else:
            self.nonbondedMethod = CutoffNonPeriodic


        ### Parameters
        self.A     = A
        self.C     = C
        self.rcut  = rcut
        u          = Universe(structure)


        if not isinstance(ff_add, list):
            ff_add = [ff_add]


        ### Make a rock file
        if rock:
            rr         = Rock(structure=rock, out='ROCK', rcut=rockrcut, Kbond=rockKbond, resname=rockresname)
            rrdms      = DesmondDMSFile(rr.dms)
            ff_add    += ['ROCK.xml']


        ### Read XML
        mapping    = ReadMappings(mapping=mapping, mapping_add=mapping_add)
        xml        = ReadXML(ff=ff, ff_add=ff_add)
        forcefield = ForceField(*xml.ff)


        ### Read a structure file
        ext = structure.split('.')[-1]
        if ext == 'pdb' or ext == 'PDB':
            pdb = PDBFile(structure)
        elif ext == 'dms' or ext == 'DMS':
            pdb = DesmondDMSFile(structure)
        else:
            assert 0 == 1, 'Please provide a pdb or dms file'



        ### Add bonds for non-protein residues in u
        if add_bonds:
            bonds, pdb = addBonds(u, xml, pdb)

        ### Combine systems (rock should be the first because of the bonds added later)
        modeller_combined = []
        universe_combined = []

        if rock:
            modeller_combined.append([rrdms.topology, rrdms.positions])
            universe_combined.append(Universe(rr.dms).atoms)

        if protein:
            proteinpdb = PDBFile(protein)
            modeller_combined.append([proteinpdb.topology, proteinpdb.positions])
            universe_combined.append(Universe(protein).atoms)

        modeller_combined.append([pdb.topology, pdb.positions])
        universe_combined.append(u.atoms)


        ### Make a modeller
        modeller = Modeller(modeller_combined[0][0], modeller_combined[0][1])
        for i in range(1, len(modeller_combined)):
            modeller.add(modeller_combined[i][0], modeller_combined[i][1])

        if pbc:
            modeller.topology.setPeriodicBoxVectors(pdb.topology.getPeriodicBoxVectors())

        self.final = modeller
        print(self.final.topology)


        ### Make a universe
        u.atoms = pd.concat(universe_combined, ignore_index=True)


        ### Create a system
        self.system = forcefield.createSystem(self.final.topology, nonbondedMethod=self.nonbondedMethod, nonbondedCutoff=self.rcut*nanometers)


        ### Add rock bonds (need to be defined after the system is defined)
        if rock:
            for i, force in enumerate(self.system.getForces()):
                if force.getName() == 'HarmonicBondForce':
                    for bond in rr.bonds:
                        force.addBond(*bond)



        ### Add posre
        if refposre:
            self.system.addForce(addRefPosre(u, refposre, fcx, fcy, fcz))
        else:
            # rock molecules have bfactor of 0.0 -> will have no restraints
            self.system.addForce(addPosre(u, bfactor_posre, fcx, fcy, fcz))



        ### Update Nonbonded + Add IsomerTorsions
        self.system.addForce(self.updateCustomNonbondedForce())
        #self.system.addForce(self.updateCustomNonbondedForceOld())
        self.removeForces(['LennardJones', 'LennardJones14', 'NonbondedForce'])

        print("Adding Isomer Torsions - started")
        self.system.addForce(addPeptideTorsions(  u, Kpeptide))
        self.system.addForce(addCisTransTorsions( u, Kcistrans, mapping, exclude=[rockresname]))
        self.system.addForce(addChiralTorsions(   u, Kchiral,   mapping, exclude=[rockresname]))
        print("Adding Isomer Torsions - finished")


        ### Run EM with Additional Torsions
        self.positions = self.final.positions
        print("-------------------------------")
        print("Running REM with additional torsions")
        self.runREM()

        ### Run EM without Additional Torsions
        print("Running REM without additional torsions")
        self.removeForces(['PeptideTorsion', 'CisTransTorsion', 'ChiralTorsion'])
        self.runREM()


        ### Save
        u.atoms[['x','y','z']] = self.numpypositions

        if rock:
            rockstruct = Universe(data = u.atoms[u.atoms.resname == 'ROCK'])
            rockstruct.dimensions = u.dimensions
            rockstruct.cell       = u.cell
            rockstruct.write('ROCK_rem.pdb')

            nonrockstruct = Universe(data = u.atoms[u.atoms.resname != 'ROCK'])
            new = Universe(data = pd.concat([Universe(rock).atoms, nonrockstruct.atoms], ignore_index=True))
            new.dimensions = u.dimensions
            new.cell       = u.cell
            if out: new.write(out)

        else:
            u.sort()
            if out: u.write(out)
        

        self.forces = { force.__class__.__name__ : force for force in self.system.getForces() }



    def removeForces(self, removes):
        for remove in removes:
            for i, force in enumerate(self.system.getForces()):
                if force.getName() == remove:
                    self.system.removeForce(i)


    def updateCustomNonbondedForce(self, excl=2):
        forces = { force.__class__.__name__ : force for force in self.system.getForces() }

        # Original Nonbonded Force (only carries charges for charmm36)
        onf = forces['NonbondedForce']

        # Custom Nonbonded Force
        # cnf.getEnergyFunction() -> 'acoef(type1, type2)/r^12 - bcoef(type1, type2)/r^6;'
        cnf = forces['CustomNonbondedForce']
        nonbondedmethod = cnf.getNonbondedMethod()
        acoef = cnf.getTabulatedFunction(0).getFunctionParameters()
        bcoef = cnf.getTabulatedFunction(1).getFunctionParameters()
        numLjTypes = acoef[0]
        epsilon = (np.array(bcoef[2]))**2 / np.array(acoef[2]) / 4
        sigma   = (np.array(acoef[2]) / np.array(bcoef[2]))**(1/6)

        # New Custom Nonbonded Force
        nnf = CustomNonbondedForce("min(rep, LJ) + coul; \
            rep  = A * (cos(pi/2 * r/sig))^2; \
            LJ   = 4 * eps * ((sig/r)^12-(sig/r)^6); \
            coul = C * q1 * q2 * (cos(pi/2 * r / rcut))^2; \
            eps  = epsilon(type1, type2); sig=sigma(type1, type2);")

        nnf.setName('REM')
        nnf.setCutoffDistance(self.rcut * nanometer)
        nnf.setNonbondedMethod(cnf.getNonbondedMethod())
        nnf.setUseLongRangeCorrection(False)

        nnf.addTabulatedFunction('epsilon', Discrete2DFunction(numLjTypes, numLjTypes, epsilon))
        nnf.addTabulatedFunction('sigma',   Discrete2DFunction(numLjTypes, numLjTypes, sigma))
        nnf.addPerParticleParameter('type')
        nnf.addPerParticleParameter('q')

        nnf.addGlobalParameter('pi',   3.141592)
        nnf.addGlobalParameter('A',    self.A * kilojoule/mole)
        nnf.addGlobalParameter('C',    self.C * kilojoule/mole)
        nnf.addGlobalParameter('rcut', self.rcut * nanometer)

        # Add particles
        for i in range(self.system.getNumParticles()):
            # NonbondedForce: ParticleParameters -> partial charge, fake sigma and epsilon
            q = onf.getParticleParameters(i)[0]

            # CustomNonbondedForce: ParticleParameters -> types
            t = cnf.getParticleParameters(i)[0]

            # Add a particle
            nnf.addParticle([t, q])

        # Add exclusion
        bonds = []
        for bond in list(self.final.topology.bonds()):
            bonds.append((bond[0].index, bond[1].index))
        nnf.createExclusionsFromBonds(bonds, excl)

        return nnf


    def updateCustomNonbondedForceOld(self, excl=2):
        forces = { force.__class__.__name__ : force for force in self.system.getForces() }

        ### Make a custom nonbonded force
        customforce = CustomNonbondedForce("min(A * (cos(pi/2 * r/sigma))^2, \
                4*epsilon*((sigma/r)^12-(sigma/r)^6)) + C * q1 * q2 * (cos(pi/2 * r / rcut))^2; \
                sigma=0.5*(sigma1+sigma2); epsilon=sqrt(epsilon1*epsilon2)")

        customforce.setName('REM')
        customforce.setCutoffDistance(self.rcut * nanometer)
        customforce.setNonbondedMethod(forces['CustomNonbondedForce'].getNonbondedMethod())
        customforce.setUseLongRangeCorrection(False)

        customforce.addPerParticleParameter('sigma')
        customforce.addPerParticleParameter('epsilon')
        customforce.addPerParticleParameter('q')

        customforce.addGlobalParameter('pi',   3.141592)
        customforce.addGlobalParameter('A',    self.A * kilojoule/mole)
        customforce.addGlobalParameter('C',    self.C * kilojoule/mole)
        customforce.addGlobalParameter('rcut', self.rcut * nanometer)

        # Add sigma, epsilon, and q
        # sigma, epsilon information is saved in CustomNonbonded
        # A = 4 * eps * sigma**12
        # B = 4 * eps * sigma**6
        atabf = forces['CustomNonbondedForce'].getTabulatedFunction(0).getFunctionParameters()
        acoef = np.array(atabf[2]).reshape(atabf[0], atabf[1])

        btabf = forces['CustomNonbondedForce'].getTabulatedFunction(1).getFunctionParameters()
        bcoef = np.array(btabf[2]).reshape(btabf[0], btabf[1])

        for i in range(self.system.getNumParticles()):
            # NonbondedForce: ParticleParameters -> partial charge, fake sigma and epsilon
            q = forces['NonbondedForce'].getParticleParameters(i)[0]

            # CustomNonbondedForce: ParticleParameters -> types
            type_index = round(forces['CustomNonbondedForce'].getParticleParameters(i)[0])
            A = acoef[type_index, type_index]
            B = bcoef[type_index, type_index]

            sigma   = np.power(A/B, 1/6)
            epsilon = B / 4 / sigma**6
            customforce.addParticle([sigma, epsilon, q])

        # Add exclusion
        bonds = []
        for bond in list(self.final.topology.bonds()):
            bonds.append((bond[0].index, bond[1].index))
        customforce.createExclusionsFromBonds(bonds, excl)

        return customforce


    def runREM(self):
        integrator = LangevinIntegrator(310*kelvin, 1/picosecond, 0.002*picoseconds)
        simulation = Simulation(self.final.topology, self.system, integrator)
        simulation.context.setPositions(self.positions)
        platform = simulation.context.getPlatform().getName()
        print("Platform: ", platform)
        print("E0: %.3e kJ/mol" %simulation.context.getState(getEnergy=True).getPotentialEnergy()._value)

        # for i in range(1, 1000 + 1):
        #     simulation.minimizeEnergy(maxIterations=1)
        #     simulation.context.setPositions(simulation.context.getState(getPositions=True).getPositions())
        #     print("E%d: %.3e kJ/mol" %(i, simulation.context.getState(getEnergy=True).getPotentialEnergy()._value))

        simulation.minimizeEnergy()
        simulation.step(1)
        print("E1: %.3e kJ/mol" %simulation.context.getState(getEnergy=True).getPotentialEnergy()._value)
        print("-------------------------------")

        self.positions      = simulation.context.getState(getPositions=True).getPositions()
        self.numpypositions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)._value * 10



    def runREMSD(self, h0, F_max_tol=1e3, nsteps=1000, verbose=True):
        '''
        https://manual.gromacs.org/current/reference-manual/algorithms/energy-minimization.html
        openMM and gromacs use kJ/mol and nm;
        h0: initial maximum displacement
        '''

        integrator = LangevinIntegrator(310*kelvin, 1/picosecond, 0.002*picoseconds)
        simulation = Simulation(self.final.topology, self.system, integrator)
        simulation.context.setPositions(self.positions)
        #print(simulation.context.getPlatform().getName())

        reports = []
        for i in range(1, nsteps+1):
            state0 = simulation.context.getState(getEnergy=True, getForces=True, getPositions=True)
            E0     = state0.getPotentialEnergy()._value
            F0     = state0.getForces(asNumpy=True)._value
            p0     = state0.getPositions(asNumpy=True)._value
            F0max  = np.linalg.norm(F0, axis=1).max()
            
            if F0max < F_max_tol:
                print(f'Fmax: {F0max:.1e} < F_max_tol: {F_max_tol:.1e}')
                break
            
            p1     = p0 + F0/F0max * h0
            simulation.context.setPositions(p1)
            state1 = simulation.context.getState(getEnergy=True, getForces=True, getPositions=True)
            E1     = state1.getPotentialEnergy()._value
            
            if E1 < E0:
                report = 'accepted'
                h0 = 1.2 * h0
                
            else:
                report = 'rejected'
                h0 = 0.2 * h0
                simulation.context.setPositions(p0)
            
            reports.append(report)
            if verbose:
                print(f'step{i:5d}: {report:s}. E0: {E0:.1e} kJ/mol. E1: {E1:.1e} kJ/mol. Fmax: {F0max:.1e}')
                
        self.positions      = simulation.context.getState(getPositions=True).getPositions()
        self.numpypositions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)._value * 10





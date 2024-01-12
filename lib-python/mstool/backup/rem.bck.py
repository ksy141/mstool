from   openmm        import *
from   openmm.app    import *
from   openmm.unit   import *
from   .universe     import Universe
from   .readmappings import ReadMappings
from   .readxml      import ReadXML
from   .protein_sel  import three2one
from   .rock         import Rock
import numpy as np

class Atom:
    def __init__(self, index):
        self.index = index

class REM:
    # pdb1 = PDBFile('villin.pdb')
    # pdb2 = PDBFile('ala_ala_ala.pdb')
    # modeller = Modeller(pdb1.topology, pdb1.positions)
    # modeller.add(pdb2.topology, pdb2.positions)
    # mergedTopology = modeller.topology
    # mergedPositions = modeller.positions

    def __init__(self, structure=None, out=None, refposre=None,
        rock=None, rockrcut=1.2, rockKbond=5000.0,
        rcut=1.2, pbc=True, 
        A=100, C=50, nsteps=100,
        mapping = [], mapping_add = [], 
        ff      = [], ff_add = [],
        Kchiral=300, Kpeptide=300, Kcistrans=300,
        fcx = 1000.0, fcy = 1000.0, fcz = 1000.0,
        bfactor_posre = 0.5):

        if pbc:
            self.nonbondedMethod = CutoffPeriodic
        else:
            self.nonbondedMethod = CutoffNonPeriodic

        self.u         = Universe(structure)
        self.resnames  = set(self.u.atoms.resname)
        self.mapping   = ReadMappings(mapping=mapping, mapping_add=mapping_add)
        self.Kchiral   = Kchiral
        self.Kpeptide  = Kpeptide
        self.Kcistrans = Kcistrans
        self.A         = A
        self.C         = C
        self.rcut      = rcut
        self.nsteps    = nsteps
        self.fcx       = fcx
        self.fcy       = fcy
        self.fcz       = fcz
        self.bfactor_posre = bfactor_posre
        self.refposre  = refposre


        ### Read a file
        ext = structure.split('.')[-1]
        if ext == 'pdb' or ext == 'PDB':
            self.pdb = PDBFile(structure)
        elif ext == 'dms' or ext == 'DMS':
            self.pdb = DesmondDMSFile(structure)
        else:
            assert 0 == 1, 'Please provide a pdb or dms file'


        ### Add ROCK.xml if rock is given
        if not isinstance(ff_add, list): ff_add = [ff_add]
        if rock: ff_add = ff_add + ['ROCK.xml']


        ### Read XML
        self.xml   = ReadXML(ff=ff, ff_add=ff_add)
        forcefield = ForceField(*self.xml.ff)


        ### Add bonds for non-protein residues in self.u
        self.addBonds()


        ### Merge if rock is given
        if rock:
            rr = Rock(structure=rock, out='ROCK', rcut=rockrcut, Kbond=rockKbond)
            rrdms = DesmondDMSFile(rr.dms)
            modeller = Modeller(rrdms.topology, rrdms.positions)
            modeller.add(self.pdb.topology, self.pdb.positions)
            self.pdb = modeller


        ### Create a system
        self.system = forcefield.createSystem(self.pdb.topology, nonbondedMethod=self.nonbondedMethod, nonbondedCutoff=self.rcut*nanometers)

        # Add posre
        if refposre:
            self.addRefPosre()
        else:
            # rock molecules have bfactor of 1.0 -> will have restraints
            self.addPosre()

        ### OLD PROCEDURE (works better...)
        #self.old_rem()

        ### NEW PROCEDURE
        #self.new_rem()

        ###
        self.updateCustomNonbondedForce2()

        ### Run EM with Additional Torsions
        self.positions = self.pdb.positions
        print("Running REM with additional torsions")
        self.addIsomerTorsions()
        print("-------------------------------")
        self.runREM()

        ### Run EM without Additional Torsions
        print("Running REM without additional torsions")
        self.removeIsomerTorsions()
        self.runREM()

        ### Save
        self.u.atoms[['x','y','z']] = self.numpypositions
        if out is not None: self.u.write(out)
        self.forces = { force.__class__.__name__ : force for force in self.system.getForces() }



    def addBonds(self):
        # add bonds except for protein residues
        protein_resnames = three2one.keys()
        for resname in self.resnames:
            # openMM PDBFile takes care of protein bonds
            if resname in protein_resnames: continue

            if resname not in self.xml.RESI.keys():
                print(f'Warning: openMM xml does not have {resname}')
                continue
            if resname not in self.mapping.RESI.keys():
                print(f'Warning: mapping does not have {resname}')
                continue

            bonds = self.xml.RESI[resname]['bonds']

            for bond in bonds:
                bA  = self.u.atoms.resname == resname
                bA0 = self.u.atoms.name == bond[0]
                bA1 = self.u.atoms.name == bond[1]

                atomA = self.u.atoms[bA & bA0].index
                atomB = self.u.atoms[bA & bA1].index
                #print(resname, bond[0], bond[1], len(atomA), len(atomB))

                assert len(atomA) == len(atomB), "the length of atoms for bonds is different"

                for a, b in zip(atomA, atomB):
                    self.pdb.topology.addBond(Atom(a), Atom(b))


    def old_rem(self):
        forces = { force.__class__.__name__ : force for force in self.system.getForces() }

        ### Make a custom nonbonded force
        customforce = CustomNonbondedForce("min(A * (cos(pi/2 * r/sigma))^2, \
                4*epsilon*((sigma/r)^12-(sigma/r)^6)) + C * q1 * q2 * (cos(pi/2 * r / rcut))^2; \
                sigma=0.5*(sigma1+sigma2); epsilon=sqrt(epsilon1*epsilon2)")

        customforce.addPerParticleParameter('sigma')
        customforce.addPerParticleParameter('epsilon')
        customforce.addPerParticleParameter('q')

        customforce.addGlobalParameter('pi', 3.141592)
        customforce.addGlobalParameter('A',    self.A * kilojoule/mole)
        customforce.addGlobalParameter('C',    self.C * kilojoule/mole)
        customforce.addGlobalParameter('rcut', self.rcut * nanometer)
        customforce.setCutoffDistance(self.rcut * nanometer)
        customforce.setNonbondedMethod((customforce.CutoffPeriodic))

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
        for bond in list(self.pdb.topology.bonds()):
            bonds.append((bond[0].index, bond[1].index))
        customforce.createExclusionsFromBonds(bonds, 2)

        ### Remove Custom-made Nonbonded
        forces = { force.__class__.__name__ : force for force in self.system.getForces() }
        self.system.removeForce(list(forces.keys()).index('CustomNonbondedForce'))

        ### Remove Custom-made Nonbonded1-4
        forces = { force.__class__.__name__ : force for force in self.system.getForces() }
        self.system.removeForce(list(forces.keys()).index('CustomBondForce'))

        ### Remove Nonbonded interactions
        forces = { force.__class__.__name__ : force for force in self.system.getForces() }
        self.system.removeForce(list(forces.keys()).index('NonbondedForce'))

        ### Add Custom-made Nonbonded
        self.system.addForce(customforce)


    def new_rem(self):
        self.updateCustomNonbondedForce2()
        #self.updateCustomNonbondedForce()
        #self.updateCustomBondForce()

        # Update a default improper dihedral (not important)
        #self.updateImproperDihedral()

        # Remove Nonbonded Force
        #self.removeNonbondedForce()


    def updateCustomNonbondedForce(self):
        forces = { force.__class__.__name__ : force for force in self.system.getForces() }

        # Original Nonbonded Force (only carries charges for charmm36)
        onf = forces['NonbondedForce']

        # Custom Nonbonded Force
        cnf = forces['CustomNonbondedForce']
        cnf.addPerParticleParameter('q')
        cnf.addGlobalParameter('pi',   3.141592)
        cnf.addGlobalParameter('A',    self.A * kilojoule/mole)
        cnf.addGlobalParameter('C',    self.C * kilojoule/mole)
        cnf.addGlobalParameter('rcut', self.rcut * nanometer)

        # transfer charge from onf to cnf
        # and then set charge = 0 for onf
        for i in range(onf.getNumParticles()):
            t = cnf.getParticleParameters(i)[0]
            q = onf.getParticleParameters(i)[0]._value
            cnf.setParticleParameters(i, [t, q])
            onf.setParticleParameters(i, 0.0, 0.0, 0.0)

        cnf.setEnergyFunction("min(rep, LJ) + coul; \
                rep = A * (cos(pi/2 * r/sig))^2; \
                LJ = 4*eps*((sig/r)^12-(sig/r)^6); \
                coul = C * q1 * q2 * (cos(pi/2 * r / rcut))^2; \
                sig = (aa/bb)^(1/6); \
                eps = bb^2 / 4 / aa; \
                aa = acoef(type1, type2); \
                bb = bcoef(type1, type2);")



    def updateCustomNonbondedForce2(self):
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
        for bond in list(self.pdb.topology.bonds()):
            bonds.append((bond[0].index, bond[1].index))
        nnf.createExclusionsFromBonds(bonds, 2)

        # Remove CustomNonbondedForce
        for i, force in enumerate(self.system.getForces()):
            if force.getName() == 'LennardJones':
                print("Removing CustomNonbondedForce")
                self.system.removeForce(i)

        # Remove CustomBondForce (1-4 interactions)
        for i, force in enumerate(self.system.getForces()):
            if force.getName() == 'LennardJones14':
                print("Removing CustomBondForce")
                self.system.removeForce(i)

        # Remove NonbondedForce
        for i, force in enumerate(self.system.getForces()):
            if force.getName() == 'NonbondedForce':
                print("Removing NonbondedForce")
                self.system.removeForce(i)

        # Add REM force
        self.system.addForce(nnf)



    # def updateCustomNonbondedForceNew(self):
    #     forces = { force.__class__.__name__ : force for force in self.system.getForces() }

    #     # Original Nonbonded Force (only carries charges for charmm36)
    #     onf = forces['NonbondedForce']

    #     # Custom Nonbonded Force
    #     cnf = CustomNonbondedForce("min(A * (cos(pi/2 * r/sigma))^2, \
    #             4*epsilon*((sigma/r)^12-(sigma/r)^6)) + C * q1 * q2 * (cos(pi/2 * r / rcut))^2; \
    #             sigma=0.5*(sigma1+sigma2); epsilon=sqrt(epsilon1*epsilon2)")

    #     cnf.setName('CustomNonbondedForceNew')
    #     cnf.setCutoffDistance(self.rcut * nanometer)
    #     cnf.setNonbondedMethod(cnf.CutoffPeriodic)
    #     cnf.addPerParticleParameter('sigma')
    #     cnf.addPerParticleParameter('epsilon')
    #     cnf.addPerParticleParameter('q')
    #     cnf.addGlobalParameter('pi',   3.141592)
    #     cnf.addGlobalParameter('A',    self.A * kilojoule/mole)
    #     cnf.addGlobalParameter('C',    self.C * kilojoule/mole)
    #     cnf.addGlobalParameter('rcut', self.rcut * nanometer)

    #     # Add sigma, epsilon, and q
    #     # sigma, epsilon information is saved in CustomNonbonded
    #     # A = 4 * eps * sigma**12
    #     # B = 4 * eps * sigma**6
    #     atabf = forces['CustomNonbondedForce'].getTabulatedFunction(0).getFunctionParameters()
    #     acoef = np.array(atabf[2]).reshape(atabf[0], atabf[1])

    #     btabf = forces['CustomNonbondedForce'].getTabulatedFunction(1).getFunctionParameters()
    #     bcoef = np.array(btabf[2]).reshape(btabf[0], btabf[1])

    #     for i in range(self.system.getNumParticles()):
    #         # NonbondedForce: ParticleParameters -> partial charge, fake sigma and epsilon
    #         q = onf.getParticleParameters(i)[0]

    #         # CustomNonbondedForce: ParticleParameters -> types
    #         type_index = round(forces['CustomNonbondedForce'].getParticleParameters(i)[0])
    #         A = acoef[type_index, type_index]
    #         B = bcoef[type_index, type_index]

    #         sigma   = np.power(A/B, 1/6)
    #         epsilon = B / 4 / sigma**6
    #         cnf.addParticle([sigma, epsilon, q])

    #     # Add exclusion
    #     bonds = []
    #     for bond in list(self.pdb.topology.bonds()):
    #         bonds.append((bond[0].index, bond[1].index))
    #     cnf.createExclusionsFromBonds(bonds, 3)

    #     self.system.addForce(cnf)


    def updateCustomBondForce(self):
        # Custom Bond Force (1-4 interactions)
        # (sigma, epsilon) -> 4*epsilon*((sigma/r)^12-(sigma/r)^6)
        # Charged interactions = 0.0
        forces = { force.__class__.__name__ : force for force in self.system.getForces() }

        cbf = forces['CustomBondForce']
        cbf.setEnergyFunction("min(rep, LJ); \
            rep = A * (cos(pi/2 * r/sigma))^2; \
            LJ  = 4 * epsilon * ((sigma/r)^12-(sigma/r)^6);")
        cbf.addGlobalParameter('pi', 3.141592)
        cbf.addGlobalParameter('A',  self.A * kilojoule/mole)


    def updateCustomNonbondedForceOld(self):
        forces = { force.__class__.__name__ : force for force in self.system.getForces() }

        customforce = CustomNonbondedForce("min(A * (cos(pi/2 * r/sigma))^2, \
                4*epsilon*((sigma/r)^12-(sigma/r)^6)) + C * q1 * q2 * (cos(pi/2 * r / rcut))^2; \
                sigma=0.5*(sigma1+sigma2); epsilon=sqrt(epsilon1*epsilon2)")

        customforce.setName('CustomNonbondedForceOld')
        customforce.addPerParticleParameter('sigma')
        customforce.addPerParticleParameter('epsilon')
        customforce.addPerParticleParameter('q')

        customforce.addGlobalParameter('pi',   3.141592)
        customforce.addGlobalParameter('A',    self.A * kilojoule/mole)
        customforce.addGlobalParameter('C',    self.C * kilojoule/mole)
        customforce.addGlobalParameter('rcut', self.rcut * nanometer)
        customforce.setCutoffDistance(self.rcut * nanometer)
        customforce.setNonbondedMethod((customforce.CutoffPeriodic))

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
        for bond in list(self.pdb.topology.bonds()):
            bonds.append((bond[0].index, bond[1].index))
        customforce.createExclusionsFromBonds(bonds, 2)

        # Remove CustomNonbondedForce and CustomBondForce
        for i, force in enumerate(self.system.getForces()):
            if force.getName() == 'CustomNonbondedForce':
                self.system.removeForce(i)

        # CustomBondForce
        for i, force in enumerate(self.system.getForces()):
            if force.getName() == 'LennardJones14':
                self.system.removeForce(i)

        # NonbondedForce
        self.removeNonbondedForce()

        # Add
        self.system.addForce(customforce)
       

    def updateImproperDihedral(self):
        # Update CustomTorsionForce (improper dihedrals)
        # k * (theta - theta0)^2
        forces = { force.__class__.__name__ : force for force in self.system.getForces() }

        for i, force in enumerate(self.system.getForces()):
            if force.getName() == 'CustomTorsionForce':
                force.setEnergyFunction('k * (acos(cos(theta-theta0)))^2')


    def runREM(self):
        integrator = LangevinIntegrator(310*kelvin, 1/picosecond, 0.002*picoseconds)
        simulation = Simulation(self.pdb.topology, self.system, integrator)
        simulation.context.setPositions(self.positions)

        print("E0: %.3e kJ/mol" %simulation.context.getState(getEnergy=True).getPotentialEnergy()._value)
        simulation.minimizeEnergy()
        simulation.step(self.nsteps)
        print("E1: %.3e kJ/mol" %simulation.context.getState(getEnergy=True).getPotentialEnergy()._value)
        print("-------------------------------")

        self.positions      = simulation.context.getState(getPositions=True).getPositions()
        self.numpypositions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True) * 10


    def addIsomerTorsions(self):
        self.system.addForce(self.addChiralTorsions())
        self.system.addForce(self.addPeptideTorsions())
        self.system.addForce(self.addCisTransTorsions())


    def addPosre(self):
        cef = CustomExternalForce("hkx*(x-x0)^2+hky*(y-y0)^2+hkz*(z-z0)^2")
        cef.addPerParticleParameter("x0")
        cef.addPerParticleParameter("y0")
        cef.addPerParticleParameter("z0")
        cef.addPerParticleParameter("hkx")
        cef.addPerParticleParameter("hky")
        cef.addPerParticleParameter("hkz")

        bA = self.u.atoms.bfactor > self.bfactor_posre
        df = self.u.atoms[bA]

        for index, row in df.iterrows():
            x0d = (row.x * angstrom).value_in_unit(nanometer)
            y0d = (row.y * angstrom).value_in_unit(nanometer)
            z0d = (row.z * angstrom).value_in_unit(nanometer)
            hfcxd = self.fcx*kilojoule_per_mole/nanometer**2
            hfcyd = self.fcy*kilojoule_per_mole/nanometer**2
            hfczd = self.fcz*kilojoule_per_mole/nanometer**2
            cef.addParticle(index,[ x0d, y0d, z0d, hfcxd,  hfcyd,  hfczd])

        print('Adding Posre for {:d} atoms whose bfactor > {:.2f})'.format(len(df), self.bfactor_posre))
        self.system.addForce(cef)


    def addRefPosre(self):
        ref = Universe(self.refposre)

        cef = CustomExternalForce("hkx*(x-x0)^2+hky*(y-y0)^2+hkz*(z-z0)^2")
        cef.addPerParticleParameter("x0")
        cef.addPerParticleParameter("y0")
        cef.addPerParticleParameter("z0")
        cef.addPerParticleParameter("hkx")
        cef.addPerParticleParameter("hky")
        cef.addPerParticleParameter("hkz")

        N = 0
        for index, atom in self.u.atoms.iterrows():
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
                hfcxd = self.fcx*kilojoule_per_mole/nanometer**2
                hfcyd = self.fcy*kilojoule_per_mole/nanometer**2
                hfczd = self.fcz*kilojoule_per_mole/nanometer**2
                cef.addParticle(index,[ x0d, y0d, z0d, hfcxd,  hfcyd,  hfczd])
                N += 1

            else:
                assert 0 == 1, 'more than one atom with the same name, resid, chain?'

        print('Adding RefPosre for {:d} atoms that exist in {:s}'.format(N, self.refposre))
        self.system.addForce(cef)


    def removeIsomerTorsions(self):
        # Need to remove it one by one because of indexing
        for i, force in enumerate(self.system.getForces()):
            if force.getName() == 'ChiralTorsion':
                self.system.removeForce(i)

        for i, force in enumerate(self.system.getForces()):
            if force.getName() == 'PeptideTorsion':
                self.system.removeForce(i)

        for i, force in enumerate(self.system.getForces()):
            if force.getName() == 'CisTransTorsion':
                self.system.removeForce(i)


    def removeNonbondedForce(self):
        # Somehow this is still required even though I set all charges = 0.0
        # I do not know why
        for i, force in enumerate(self.system.getForces()):
            if force.getName() == 'NonbondedForce':
                self.system.removeForce(i)


    def addChiralTorsions(self):
        ctf = CustomTorsionForce('k * (acos(cos(theta-theta0)))^2')
        #ctf = CustomTorsionForce('k * (cos(theta) - cos(theta0))^2')
        #ctf = CustomTorsionForce('k * (theta-theta0)^2')
        #ctf = CustomTorsionForce('-k * (1 + cos(theta-theta0))')
        ctf.setName('ChiralTorsion')
        ctf.addPerTorsionParameter("k")
        ctf.addPerTorsionParameter("theta0")

        N = 0
        for resname in self.resnames:
            chirals = self.mapping.RESI[resname]['chiral']
            for chiral in chirals:
                bA  = self.u.atoms.resname == resname
                bA0 = self.u.atoms.name    == chiral[0]
                bA1 = self.u.atoms.name    == chiral[1]
                bA2 = self.u.atoms.name    == chiral[2]
                bA3 = self.u.atoms.name    == chiral[3]              
                bA4 = self.u.atoms.name    == chiral[4]

                target = self.u.atoms[bA & bA0].index
                center = self.u.atoms[bA & bA1].index
                atomC  = self.u.atoms[bA & bA2].index
                atomD  = self.u.atoms[bA & bA3].index
                atomE  = self.u.atoms[bA & bA4].index

                assert len(target) == len(center) == \
                    len(atomC) == len(atomD) == len(atomE), \
                    "the length of atoms for chiral torsions is different"

                for a, b, c, d, e in zip(target, center, atomC, atomD, atomE):
                    ctf.addTorsion(b, c, d, e, [self.Kchiral, -35.26 * 3.141592/180])
                    ctf.addTorsion(b, d, e, c, [self.Kchiral, -35.26 * 3.141592/180])
                    ctf.addTorsion(b, e, c, d, [self.Kchiral, -35.26 * 3.141592/180])

                    ctf.addTorsion(a, c, d, e, [self.Kchiral, -70.53 * 3.141592/180])
                    ctf.addTorsion(a, d, e, c, [self.Kchiral, -70.53 * 3.141592/180])
                    ctf.addTorsion(a, e, c, d, [self.Kchiral, -70.53 * 3.141592/180])
                    N += 1

        print('Adding ChiralTorsion for {:d} chirals'.format(N))
        return ctf

    def addCisTransTorsions(self):
        ctf = CustomTorsionForce('k * (acos(cos(theta-theta0)))^2')
        #ctf = CustomTorsionForce('k * (cos(theta) - cos(theta0))^2')
        #ctf = CustomTorsionForce('k * (theta-theta0)^2')
        ctf.setName('CisTransTorsion')
        ctf.addPerTorsionParameter("k")
        ctf.addPerTorsionParameter("theta0")

        N = 0
        for resname in self.resnames:
            for isomer in ['cis', 'trans']:
                atomset = self.mapping.RESI[resname][isomer]
                for atoms in atomset:
                    bA  = self.u.atoms.resname == resname
                    bA0 = self.u.atoms.name    == atoms[0]
                    bA1 = self.u.atoms.name    == atoms[1]
                    bA2 = self.u.atoms.name    == atoms[2]
                    bA3 = self.u.atoms.name    == atoms[3]

                    atomA = self.u.atoms[bA & bA0].index
                    atomB = self.u.atoms[bA & bA1].index
                    atomC = self.u.atoms[bA & bA2].index
                    atomD = self.u.atoms[bA & bA3].index

                    assert len(atomA) == len(atomB) == \
                        len(atomC) == len(atomD), \
                        "the length of atoms for cistrans torsions is different"

                    for a, b, c, d in zip(atomA, atomB, atomC, atomD):
                        N += 1
                        if isomer == 'cis':
                            ctf.addTorsion(a, b, c, d, [self.Kcistrans, 0.000])
                            ctf.addTorsion(a, c, b, d, [self.Kcistrans, 0.000])
                        elif isomer == 'trans':
                            ctf.addTorsion(a, b, c, d, [self.Kcistrans, 3.141592])
                            ctf.addTorsion(a, c, b, d, [self.Kcistrans, 3.141592])

        print('Adding CisTransTorsion for {:d} isomers'.format(N))
        return ctf


    def addPeptideTorsions(self):
        ctf = CustomTorsionForce('k * (acos(cos(theta-theta0)))^2')
        #ctf = CustomTorsionForce('k * (cos(theta) - cos(theta0))^2')
        #ctf = CustomTorsionForce('k * (theta-theta0)^2')
        ctf.setName('PeptideTorsion')
        ctf.addPerTorsionParameter("k")
        ctf.addPerTorsionParameter("theta0")

        protein_resnames = three2one.keys()
        bA  = self.u.atoms.resname.isin(protein_resnames)
        bAN = self.u.atoms.name == 'N'
        bAH = (self.u.atoms.name.isin(['HN', 'H'])) | ((self.u.atoms.resname == 'PRO') & (self.u.atoms.name == 'CD'))
        bAC = self.u.atoms.name == 'C'
        bAO = self.u.atoms.name == 'O'

        N = 0
        Natoms = self.u.atoms[bA & bAN]
        Hatoms = self.u.atoms[bA & bAH]
        Catoms = self.u.atoms[bA & bAC]
        Oatoms = self.u.atoms[bA & bAO]

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
            ctf.addTorsion(OatomIndex, CatomIndex, NatomIndex, HatomIndex, [self.Kpeptide, 3.141592])
            N += 1

        print('Adding PeptideTorsion for {:d} isomers'.format(N))
        return ctf


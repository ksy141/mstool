from   openmm        import *
from   openmm.app    import *
from   openmm.unit   import *

from   .universe     import Universe
from   .readmappings import ReadMappings
from   .readxml      import ReadXML
from   .checktetrahedron import CheckTetrahedron

from   ..utils.protein_sel import three2one
from   ..utils.rock        import Rock
from   ..utils.openmmutils import *

import numpy  as np
import pandas as pd
import sys
sys.setrecursionlimit(1000000)


class REM:
    # pdb1 = PDBFile('villin.pdb')
    # pdb2 = PDBFile('ala_ala_ala.pdb')
    # modeller = Modeller(pdb1.topology, pdb1.positions)
    # modeller.add(pdb2.topology, pdb2.positions)
    # mergedTopology = modeller.topology
    # mergedPositions = modeller.positions

    def __init__(self, structure=None, out=None, protein=None, refposre=None, outrem=None,
        rock=None, rockrcut=1.2, rockKbond=5000.0, rockresname='ROCK', rockout='ROCK_rem.pdb',
        nonrockout='NONROCK.dms',
        rocktype='HL',
        rcut=1.2, pbc=True, 
        A=100, C=50,
        mapping = [], mapping_add = [], 
        ff      = [], ff_add = [],
        Kchiral=300, Kpeptide=300, Kcistrans=300, Kdihedral=300,
        fcx = 1000.0, fcy = 1000.0, fcz = 1000.0,
        bfactor_posre = 0.5, add_bonds=True, sort=False, version='v4',
        cospower=2, turn_off_torsion_warning=False,
        nsteps=10000,
        turn_off_EMNVT=False):
        
        # v3 should not be used
        # protein: version = 'v4' seems the best
        # lipid:   version = 'v1' seems the best


        ### NonbondedMethod
        if pbc:
            self.nonbondedMethod = CutoffPeriodic
        else:
            self.nonbondedMethod = CutoffNonPeriodic


        ### Parameters
        self.A     = A
        self.C     = C
        self.rcut  = rcut
        self.cospower = cospower
        self.nsteps   = nsteps
        u          = Universe(structure)


        if not isinstance(ff_add, list):
            ff_add = [ff_add]


        ### Make a rock file
        if rock:
            rr         = Rock(structure=rock, out='ROCK', 
                              rcut=rockrcut, 
                              Kbond=rockKbond, 
                              resname=rockresname,
                              type=rocktype)

            rrdms      = DesmondDMSFile(rr.dms)
            ff_add    += ['ROCK.xml']


        ### Read XML
        mapping    = ReadMappings(mapping=mapping, mapping_add=mapping_add)
        xml        = ReadXML(ff=ff, ff_add=ff_add)
        forcefield = ForceField(*xml.ff)
        self.forcefield = forcefield


        ### Read a structure file
        if structure:
            ext = structure.split('.')[-1]
            if ext == 'pdb' or ext == 'PDB':
                pdb = PDBFile(structure)
                # already build standard bonds
                ### Add bonds for non-protein residues in u
                if add_bonds:
                    bonds, pdb = addBonds(u, xml, pdb)

            if ext == 'dms' or ext == 'DMS':
                #from .dmsfile import DMSFile
                #pdb = DMSFile(structure)
                pdb   = DesmondDMSFile(structure)
                bonds = getBonds(structure, ff=ff, ff_add=ff_add)
                pdbatoms = [atom for atom in pdb.topology.atoms()]
                for bond in bonds:
                    pdb.topology.addBond(pdbatoms[bond[0]], pdbatoms[bond[1]])

        else:
            raise IOError('Please provide a pdb or dms file')


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
        
        if structure:
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


        ### Bonds (which whill not include ROCK)
        self.bonds = []
        for bond in self.final.topology.bonds():
            i0 = bond[0].index
            i1 = bond[1].index
            if i0 > i1: i0, i1 = i1, i0
            if [i0, i1] not in self.bonds:
                self.bonds.append([i0, i1])
        u.bonds = self.bonds


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
        if version == 'v1':
            print("-------------------------------")
            print('using REM version 1')
            self.system.addForce(self.updateCustomNonbondedForceOld())
            self.removeForces(['LennardJones', 'LennardJones14', 'NonbondedForce'])
        elif version == 'v2':
            print('using REM version 2')
            self.system.addForce(self.updateCustomNonbondedForce(excl=2))
            self.removeForces(['LennardJones', 'LennardJones14', 'NonbondedForce'])
        elif version == 'v3':
            print('using REM version 3')
            self.updateCustomBondForce()
            self.removeForces(['LennardJones', 'NonbondedForce']) #LennardJones14
        else:
            print('using REM version 4')
            self.system.addForce(self.updateCustomNonbondedForce(excl=3))
            self.updateCustomBondForce()
            self.removeForces(['LennardJones', 'NonbondedForce']) #LennardJones14
 
        print("Adding Isomer Torsions - started")
        self.system.addForce(addPeptideTorsions(  u, Kpeptide))
        self.system.addForce(addCisTransTorsions( u, Kcistrans, mapping, exclude=[rockresname], turn_off_torsion_warning=turn_off_torsion_warning))
        self.system.addForce(addChiralTorsions(   u, Kchiral,   mapping, exclude=[rockresname], turn_off_torsion_warning=turn_off_torsion_warning))
        self.system.addForce(addDihedralTorsions( u, Kdihedral, mapping, exclude=[rockresname], turn_off_torsion_warning=turn_off_torsion_warning))
        print("Adding Isomer Torsions - finished")


        ### Run REM with Additional Torsions
        self.positions = self.final.positions
        print("Running REM with isomeric torsions")
        self.runREM()

        ### Run REM without Additional Torsions
        print("Running REM without isomeric torsions")
        self.removeForces(['PeptideTorsion', 'CisTransTorsion', 'ChiralTorsion'])
        self.runREM()

        ### Save outREM
        if outrem:
            u.atoms[['x','y','z']] = self.numpypositions
            u.write(outrem)
            #CheckTetrahedron(outrem, ff=ff, ff_add=ff_add)

        ### Run EM + NVT
        if turn_off_EMNVT:
            print("EM+NVT is turned off because turn_off_EMNVT=True")
        elif not rock:
            print(f"Running EM+NVT with unmodified force field\nand without isomeric torsions for {self.nsteps} steps")
            self.runEMNVT()
        elif rock:
            print("EM+NVT is turned off with ROCK (AA)")

        ### Save
        u.atoms[['x','y','z']] = self.numpypositions

        if rock:
            rockstruct = Universe(data = u.atoms[u.atoms.resname == 'ROCK'])
            rockstruct.dimensions = u.dimensions
            rockstruct.cell       = u.cell 
            rockstruct.bodns      = addBonds(rockstruct, xml)
            rockstruct.write(rockout)

            nonrockstruct = Universe(data = u.atoms[u.atoms.resname != 'ROCK'])
            nonrockstruct.dimensions = u.dimensions
            nonrockstruct.cell       = u.cell
            nonrockstruct.bonds      = addBonds(nonrockstruct, xml)
            nonrockstruct.write(nonrockout)

            new = Universe(data = pd.concat([Universe(rock).atoms, nonrockstruct.atoms], ignore_index=True))
            new.dimensions = u.dimensions
            new.cell       = u.cell
            if out: new.write(out)
            ff_add.remove('ROCK.xml')

        else:
            #if sort: u.sort()
            #if out.split('.')[-1] == 'dms':
                #u.update_anum_mass()
                #u.bonds = addBonds(u, xml)
                #u.bonds = self.bonds

            if out:
                u.write(out, wrap=True)
        
        #CheckTetrahedron(out, ff=ff, ff_add=ff_add)
        self.universe = u
        self.forces   = { force.__class__.__name__ : force for force in self.system.getForces() }
        #print(self.forces.keys())



    def removeForces(self, removes):
        for remove in removes:
            for i, force in enumerate(self.system.getForces()):
                if force.getName() == remove:
                    self.system.removeForce(i)


    def updateCustomNonbondedForce(self, excl=3):
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
        nnf = CustomNonbondedForce(f"min(rep, LJ) + coul; \
            rep  = A * (cos(pi/2 * r/sig))^{self.cospower}; \
            LJ   = 4 * eps * ((sig/r)^12-(sig/r)^6); \
            coul = C * q1 * q2 * (cos(pi/2 * r / rcut))^{self.cospower}; \
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


    def updateCustomBondForce(self):
        # Custom Bond Force (1-4 interactions)
        # (sigma, epsilon) -> 4*epsilon*((sigma/r)^12-(sigma/r)^6)
        # Charged interactions = 0.0

        forces = { force.__class__.__name__ : force for force in self.system.getForces() }

        if 'CustomBondForce' in forces.keys():
            cbf = forces['CustomBondForce']
            cbf.setEnergyFunction(f"min(rep, LJ); \
                rep = A * (cos(pi/2 * r/sigma))^{self.cospower}; \
                LJ  = 4 * epsilon * ((sigma/r)^12-(sigma/r)^6);")
            cbf.addGlobalParameter('pi', 3.141592)
            cbf.addGlobalParameter('A',  self.A * kilojoule/mole)



    def updateCustomNonbondedForceOld(self, excl=2):
        forces = { force.__class__.__name__ : force for force in self.system.getForces() }

        ### Make a custom nonbonded force
        customforce = CustomNonbondedForce(f"min(A * (cos(pi/2 * r/sigma))^{self.cospower}, \
                4*epsilon*((sigma/r)^12-(sigma/r)^6)) + C * q1 * q2 * (cos(pi/2 * r / rcut))^{self.cospower}; \
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
        print("-------------------------------")
        print("Platform: ", platform)
        print("E0: %.3e kJ/mol" %simulation.context.getState(getEnergy=True).getPotentialEnergy()._value)
        simulation.minimizeEnergy()
        print("E1: %.3e kJ/mol" %simulation.context.getState(getEnergy=True).getPotentialEnergy()._value)
        print("-------------------------------")

        self.positions      = simulation.context.getState(getPositions=True).getPositions()
        self.numpypositions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)._value * 10


    def runEMNVT(self):
        integrator = LangevinIntegrator(310*kelvin, 1/picosecond, 0.002*picoseconds)
        system     = self.forcefield.createSystem(self.final.topology, 
                                                  nonbondedMethod=self.nonbondedMethod, 
                                                  nonbondedCutoff=self.rcut*nanometers)
        simulation = Simulation(self.final.topology, system, integrator)
        simulation.context.setPositions(self.positions)
        #print(simulation.context.getState().getPeriodicBoxVectors(asNumpy=True)._value * 10)
        platform = simulation.context.getPlatform().getName()
        print("-------------------------------")
        print("Platform: ", platform)
        print("E0: %.3e kJ/mol" %simulation.context.getState(getEnergy=True).getPotentialEnergy()._value)
        simulation.minimizeEnergy()
        simulation.step(self.nsteps)
        print("E1: %.3e kJ/mol" %simulation.context.getState(getEnergy=True).getPotentialEnergy()._value)
        print("-------------------------------")

        self.positions      = simulation.context.getState(getPositions=True).getPositions()
        self.numpypositions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)._value * 10



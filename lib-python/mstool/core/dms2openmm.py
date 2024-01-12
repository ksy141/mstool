import sqlite3
import pandas as pd
from   .universe import Universe
from   openmm.app import *
from   openmm import *
from   openmm.unit import *

class DMS2openmm:
    def __init__(self, dms_in, nonbondedMethod='CutoffPeriodic', nonbondedCutoff=1.1, soft=False, A=200.0, C=50.0):

        ### Determine nonbondedMethod
        if nonbondedMethod   == 'CutoffPeriodic':
            nbM = CutoffPeriodic
        elif nonbondedMethod == 'NoCutoff':
            nbM = NoCutoff
        elif nonbondedMethod == 'CutoffNonPeriodic':
            nbM = CutoffNonPeriodic
        elif nonbondedMethod == 'PME':
            nbM = PME
        else:
            raise AssertionError('unrecognized nonbonded method')

        ### Nonbonded parameters
        self.soft = soft
        self.A    = A
        self.C    = C


        ### read DMS 
        self.conn = sqlite3.connect(dms_in, isolation_level=None, detect_types=sqlite3.PARSE_COLNAMES)


        ### DMS to openMM
        self.dms       = DesmondDMSFile(dms_in, verbose=False)
        self.topology  = self.dms.topology
        self.system    = self.dms.createSystem(
            nonbondedMethod = nbM,
            nonbondedCutoff = nonbondedCutoff * nanometer)


        ### Make NBFIX force
        self.Nonbonded(nonbondedMethod, nonbondedCutoff)

        ### Fetch cosine angle
        self.Angles()

        ### Correct Proper Dihedrals
        self.Dihedrals()

        ### Correct Improper Dihedrals
        self.Impropers()

        ### turn off DispersionCorrection (default = True)
        self.forces = { force.__class__.__name__ : force for force in self.system.getForces() }
        #self.forces['NonbondedForce'].setUseDispersionCorrection(False)



    def Nonbonded(self, nonbondedMethod, nonbondedCutoff):
        # read DMS - nonbonded_combined_param
        df         = pd.read_sql_query("SELECT * FROM nonbonded_combined_param", self.conn)
        numLjTypes = len(df.param1.unique())
        epsilon    = [0] * (numLjTypes * numLjTypes)
        sigma      = epsilon[:]

        for m in range(numLjTypes):
            for n in range(numLjTypes):
                bA1 = (df.param1 == m) & (df.param2 == n)
                bA2 = (df.param1 == n) & (df.param2 == m)
                s   = df[bA1 | bA2].sigma   * 0.1   #A to nm
                e   = df[bA1 | bA2].epsilon * 4.184 #kcal to kJ

                sigma[m   + numLjTypes * n] = s
                epsilon[m + numLjTypes * n] = e


        # Make NBFIX
        if not self.soft:
            force = CustomNonbondedForce('4*eps*((sig/r)^12-(sig/r)^6) + 138.911*q1*q2/r; eps=epsilon(type1, type2); sig=sigma(type1, type2);')
            force.addTabulatedFunction('epsilon', Discrete2DFunction(numLjTypes, numLjTypes, epsilon))
            force.addTabulatedFunction('sigma',   Discrete2DFunction(numLjTypes, numLjTypes, sigma))
            force.addPerParticleParameter('type')
            force.addPerParticleParameter('q')


        if self.soft:
            force = CustomNonbondedForce("min(rep, LJ) + coul; \
                rep = A * (cos(pi/2 * r/sig))^2; \
                LJ = 4*eps*((sig/r)^12-(sig/r)^6); \
                coul = C * q1 * q2 * (cos(pi/2 * r / rcut))^2; \
                eps=epsilon(type1, type2); sig=sigma(type1, type2);")

            force.addTabulatedFunction('epsilon', Discrete2DFunction(numLjTypes, numLjTypes, epsilon))
            force.addTabulatedFunction('sigma',   Discrete2DFunction(numLjTypes, numLjTypes, sigma))
            force.addPerParticleParameter('type')
            force.addPerParticleParameter('q')

            force.addGlobalParameter('pi', 3.141592)
            force.addGlobalParameter('A',  self.A * kilojoule/mole)
            force.addGlobalParameter('C',  self.C * kilojoule/mole)
            force.addGlobalParameter('rcut', nonbondedCutoff * nanometer)

        force.setName('Nonbonded')
        force.setUseLongRangeCorrection(False)


        if nonbondedMethod == 'CutoffPeriodic':
            force.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
            force.setCutoffDistance(nonbondedCutoff * nanometer)
        elif nonbondedMethod == 'NoCutOff':
            force.setNonbondedMethod(CustomNonbondedForce.NoCutoff)
        elif nonbondedMethod == 'CutoffNonPeriodic':
            force.setNonbondedMethod(CustomNonbondedForce.CutoffNonPeriodic)
            force.setCutoffDistance(nonbondedCutoff * nanometer)
        else:
            raise AssertionError('unrecognized nonbonded method')

        # read DMS - particle
        df = pd.read_sql_query("SELECT * FROM PARTICLE", self.conn)

        # add nbtype
        for index, atom in df.iterrows():
            force.addParticle((atom.nbtype,atom.charge))

        # add exclusion
        forces = { force.__class__.__name__ : force for force in self.system.getForces() }
        for i in range(forces['NonbondedForce'].getNumExceptions()):
            exception = forces['NonbondedForce'].getExceptionParameters(i)
            force.addExclusion(exception[0], exception[1]) 

        # remove force
        self.system.removeForce(list(forces.keys()).index('NonbondedForce'))

        # add force
        self.system.addForce(force)



    def Angles(self):

        force = CustomAngleForce("0.5 * k * (cos(theta) - cos_theta0)^2;")
        force.addPerAngleParameter("cos_theta0");
        force.addPerAngleParameter("k");

        # fetch
        q = """SELECT p0, p1, p2, cos_theta0, fc 
        FROM angle_harmcos_term INNER JOIN angle_harmcos_param 
        ON angle_harmcos_term.param=angle_harmcos_param.id"""
        angs = self.conn.execute(q).fetchall()

        for ang in angs:
            p0 = ang[0]
            p1 = ang[1]
            p2 = ang[2]
            cos_theta0 = ang[3]
            fc = ang[4] * 2 * kilocalorie_per_mole
            force.addAngle(p0, p1, p2, (cos_theta0, fc))

        self.system.addForce(force)


    def Dihedrals(self):
        # remove force
        # openmm has an incorrect converter for PeriodicTorsionsToSystem
        # incorrect:      sum_{n=0-6} (kn * (1 + cos(n * theta - theta0)))
        # correct:   k0 + sum_{n=1-6} (kn * cos(n * theta - theta0))
        # plus, openmm complains that it cannot run simulations because of n=0
        forces = { force.__class__.__name__ : force for force in self.system.getForces() }
        self.system.removeForce(list(forces.keys()).index('PeriodicTorsionForce'))

        # Make a CustomTorsionForce instance
        energy  = 'fc0'
        for n in range(1, 7):
            energy += f' + fc{n} * cos({n} * theta - phi0)'
        force = CustomTorsionForce(energy)

        # Set a name because openmm uses CustomTorsionForce to make an improper dihedral
        force.setName('ProperDihedrals')

        # Add parameters
        force.addPerTorsionParameter('phi0')
        for n in range(0, 7):
            force.addPerTorsionParameter(f'fc{n}')

        # Fetch
        q = """SELECT p0, p1, p2, p3, phi0, fc0, fc1, fc2, fc3, fc4, fc5, fc6
        FROM dihedral_trig_term INNER JOIN dihedral_trig_param
        ON dihedral_trig_term.param=dihedral_trig_param.id"""
        dihes = self.conn.execute(q).fetchall()

        for dihe in dihes:
            p0   = dihe[0]
            p1   = dihe[1]
            p2   = dihe[2]
            p3   = dihe[3]
            phi0 = dihe[4]  * degree
            fc0  = dihe[5]  * kilocalorie_per_mole
            fc1  = dihe[6]  * kilocalorie_per_mole
            fc2  = dihe[7]  * kilocalorie_per_mole
            fc3  = dihe[8]  * kilocalorie_per_mole
            fc4  = dihe[9]  * kilocalorie_per_mole
            fc5  = dihe[10] * kilocalorie_per_mole
            fc6  = dihe[11] * kilocalorie_per_mole
            force.addTorsion(p0, p1, p2, p3, (phi0, fc0, fc1, fc2, fc3, fc4, fc5, fc6))

        self.system.addForce(force)


    def Impropers(self):
        # periodicity at phi0 +- 180 is not taken into account in openMM.
        for i, force in enumerate(self.system.getForces()):
            if force.getName() == 'CustomTorsionForce':
                force.setEnergyFunction('k * (acos(cos(theta-theta0)))^2')
                force.setName('ImproperDihedrals')


    def make(self):
        return self.system, self.dms

        



"""
dmsfile.py: Load Desmond dms files and run EM/NVT/NPT

Portions copyright (c) 2013 Stanford University and the Authors
Authors: Robert McGibbon
Contributors: Emilio Gallicchio, Baofeng Zhang, Tony Zhao

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.

Further modified by Siyoung Kim
- Fix PeriodicTorsions
- Change ImproperTorsions from k(theta-theta0)^2 -> k(acos(0.99*cos(theta-theta0)))^2 
  to avoid discontinuity at +-180
- Nonbonded using TabulatedFunction (Supports "nonbonded_combined_param")
- Reduced Nonbonded Energy Minimization (REM) is supported
- Energy Minimization + NVT + NPT can be run
"""
from __future__ import absolute_import
from __future__ import print_function

import os
import math
import numpy as np

import openmm as mm
from openmm     import MonteCarloBarostat, MonteCarloMembraneBarostat
from openmm.app import forcefield as ff
from openmm.app import Element, Topology, PDBFile, Simulation
from openmm.app import DCDReporter, StateDataReporter
from openmm.app.element import hydrogen
from openmm.unit import (nanometer, angstrom, dalton, radian,
                        kilocalorie_per_mole, kilojoule_per_mole,
                        degree, elementary_charge, femtosecond, 
                        picosecond, picoseconds, kelvin, bar)

from openmm.app.amberprmtopfile import HCT
from openmm.app.internal.customgbforces import GBSAHCTForce
import sqlite3
import shutil

# {'NoCutoff': 0, 'CutoffNonPeriodic': 1, 'CutoffPeriodic': 2}
methodMap = {'NoCutoff':          mm.NonbondedForce.NoCutoff,
             'CutoffNonPeriodic': mm.NonbondedForce.CutoffNonPeriodic,
             'CutoffPeriodic':    mm.NonbondedForce.CutoffPeriodic,
             'Ewald':             mm.NonbondedForce.Ewald,
             'PME':               mm.NonbondedForce.PME,
             'LJPME':             mm.NonbondedForce.LJPME}


class DMSFile(object):
    """DesmondDMSFile parses a Desmond DMS (desmond molecular system) and
    constructs a topology and (optionally) an OpenMM System from it
    """

    def __init__(self, file, verbose=False):
        """Load a DMS file

        Parameters
        ----------
        file : string (for one file)
        or
        file : list of strings (multiple files, each containing a molecule)
        the name(s) of the file to load
        """

        # sqlite3 is included in the standard lib, but at python
        # compile time, you can disable support (I think), so it's
        # not *guarenteed* to be available. Doing the import here
        # means we only raise an ImportError if people try to use
        # this class, so the module can be safely imported

        self._verbose = verbose
        self._file = []
        self._open = []
        self._tables = []
        self._conn = []
        self._provenance = []
        self._natoms = [] #number of atoms in each file

        #constructs list of files even if only one file
        if isinstance(file, list):
            self._file = file
        else:
            self._file.append(file)

        #open sqlite connectors for each of the files
        nfiles = len(self._file)
        for (fcounter,f) in zip(range(0,nfiles),self._file):
            self._open.append(False)
            if not  os.path.exists(str(f)):
                raise IOError("No such file or directory: '%s'" % str(f))
            conn = sqlite3.connect(f)
            self._conn.append(conn)
            self._open[fcounter] = True
            tables = self._readSchemas(conn)
            if len(tables) == 0:
                raise IOError('DMS file %s was not loaded sucessfully. No tables found' % str(f))
            if 'nbtype' not in tables['particle']:
                raise ValueError('No nonbonded parameters associated with '
                                 'DMS file %s. You can add a forcefield with the '
                                 'viparr command line tool distributed with desmond' % str(f))
            self._tables.append(tables)
            self._natoms.append(0)

        # Build the topology
        self.topology, self.positions, self.velocities = self._createTopology()
        self._topologyAtoms = list(self.topology.atoms())
        self._atomBonds = [{} for x in range(len(self._topologyAtoms))]
        self._angleConstraints = [{} for x in range(len(self._topologyAtoms))]

        #offset to each molecule/file in the topology
        self._offset = self._prefixsum(self._natoms)

        #retrieve provenance
        self._createProvenance()

    def getPositions(self):
        """Get the positions of each atom in the system
        """
        return self.positions

    def getVelocities(self):
        """Get the positions of each atom in the system
        """
        return self.velocities

    def getTopology(self):
        """Get the topology of the system
        """
        return self.topology

    def getProvenance(self):
        """Get the provenance string of this system
        """
        if len(self._conn) == 1:
            return self._provenance[0]
        else:
            return self._provenance

    def _createTopology(self):
        """Build the topology of the system
        """
        top = Topology()
        positions = []
        velocities = []
        boxVectors = []
        self.cell  = []

        #assume cell dimensions are set in the first file
        #the other molecules inherit the same cell
        conn = self._conn[0]
        for x, y, z in conn.execute('SELECT x, y, z FROM global_cell'):
            boxVectors.append(mm.Vec3(x, y, z))
            self.cell.append([x * 0.1, y * 0.1, z * 0.1])

        unitCellDimensions = [boxVectors[0][0], boxVectors[1][1], boxVectors[2][2]]
        top.setUnitCellDimensions(unitCellDimensions*angstrom)

        #process each file
        nfiles = len(self._conn)
        for (fcounter, conn, tables) in zip(range(0,nfiles),self._conn,self._tables):

            atoms = {}
            lastChain = None
            lastResId = None
            c = top.addChain()
            q = """SELECT id, name, anum, resname, resid, chain, x, y, z, vx, vy, vz, mass
                FROM particle ORDER BY id"""
            for (atomId, atomName, atomNumber, resName, resId, chain, x, y, z, vx, vy, vz, mass) in conn.execute(q):
                newChain = False
                if chain != lastChain:
                    lastChain = chain
                    c = top.addChain()
                    newChain = True
                if resId != lastResId or newChain:
                    lastResId = resId
                    if resName in PDBFile._residueNameReplacements:
                        resName = PDBFile._residueNameReplacements[resName]
                    r = top.addResidue(resName, c)
                    if resName in PDBFile._atomNameReplacements:
                        atomReplacements = PDBFile._atomNameReplacements[resName]
                    else:
                        atomReplacements = {}

                if atomNumber == 0 and atomName.startswith('Vrt'):
                    elem = None
                else:
                    elem = Element.getByAtomicNumber(atomNumber)

                if atomName in atomReplacements:
                    atomName = atomReplacements[atomName]

                atoms[atomId] = top.addAtom(atomName, elem, r)
                positions.append(mm.Vec3(x, y, z))

                velocities.append(mm.Vec3(vx, vy, vz))

            self._natoms[fcounter] = len(atoms)

            for p0, p1 in conn.execute('SELECT p0, p1 FROM bond'):
                top.addBond(atoms[p0], atoms[p1])

        positions = positions*angstrom
        velocities = velocities*angstrom/femtosecond
        return top, positions, velocities

    def setPositions(self, positions):
        """Update atomic positions in attached DMS files
        """
        q = """UPDATE particle SET x = ?1, y = ?2, z = ?3 WHERE id == ?4"""
        #file counter
        iat = 0 #global atom counter
        for (fcounter,conn,tables,offset) in self._localVars():
            natoms = self._natoms[fcounter]
            for iat_in_file in range(0,natoms):
                vec = positions[iat]
                try:
                    (x, y , z) = vec.value_in_unit(angstrom)
                except:
                    (x, y , z) = vec
                conn.execute(q, (x,y,z,iat_in_file))
                iat += 1
            conn.commit()

        return iat

    def setVelocities(self, velocities):
        """Update atomic velocities in attached DMS files
        """
        q = """UPDATE particle SET vx = ?1, vy = ?2, vz = ?3 WHERE id == ?4"""
         #file counter
        iat = 0 #global atom counter
        for (fcounter,conn,tables,offset) in self._localVars():
            natoms = self._natoms[fcounter]
            for iat_in_file in range(0,natoms):
                vec = velocities[iat]
                (vx, vy , vz) = vec.value_in_unit(angstrom/femtosecond)
                conn.execute(q, (vx,vy,vz,iat_in_file))
                iat += 1
            conn.commit()

        return iat

    def _get_gb_params(self):
        """
        get charge, radius, screened_radius from hct table in the .dms files
        --radiusN=radius-0.009  # Offset radius, the radius needs to be converted to nanometer unit
        --screenN=screened_radius*radiusN #Scaled offset radius
        """
        length_conv = u.angstrom.conversion_factor_to(u.nanometer)
        q="""SELECT charge,radius,screened_radius FROM hct ORDER BY id"""

        gb_p = []

        try:
            for (fcounter,conn,tables,offset) in self._localVars():
                for charge,radius,screened_radius in conn.execute(q):
                    chargeN = charge
                    radiusN = radius
                    screenN = screened_radius
                    radiusN *=length_conv
                    radiusN -=0.009
                    screenN *=radiusN
                    gb_p.append((chargeN,radiusN,screenN))
            return gb_p
        except:
            return None


    def _get_agbnp2_params(self):
        """
        get charge, radius, etc. from AGBNP2 table in the .dms file and computes AGBNP2 parameters for each atom
        """
        length_conv = u.angstrom.conversion_factor_to(u.nanometer)
        en_conv = u.kilocalorie_per_mole.conversion_factor_to(u.kilojoule_per_mole)
        gamma_conv = en_conv/(length_conv*length_conv)
        alpha_conv = en_conv*length_conv*length_conv*length_conv

        # gather AGBNP2 parameters from agbnp2 table
        #   atomic number and charge from particle table matching atom id
        #   sigma and epsilon from nonbonded_param table matching hbtype from particle
        q="""SELECT anum,charge,radius,igamma,ialpha,idelta,sgamma,salpha,sdelta,hbtype,hbw,sigma,epsilon from particle INNER JOIN agbnp2 ON particle.id==agbnp2.id INNER JOIN nonbonded_param ON particle.nbtype==nonbonded_param.id ORDER BY particle.id"""

        gb_p = []

        try:
            for (fcounter,conn,tables,offset) in self._localVars():
                for anum,charge,radius,igamma,ialpha,idelta,sgamma,salpha,sdelta,hbtype,hbw,sigma,epsilon in conn.execute(q):
                    if anum == 1:
                        ishydrogenN = 1
                    else:
                        ishydrogenN = 0
                    radiusN = length_conv*radius
                    chargeN = charge
                    gammaN = gamma_conv*(igamma+0.*sgamma)#AGBNP must have only one gamma
                    alphaN = alpha_conv*(ialpha+salpha)
                    # delta parameter is ignored
                    hbwN = en_conv * hbw
                    gb_p.append([radiusN,chargeN,gammaN,alphaN,hbtype,hbwN,ishydrogenN])
            return gb_p
        except:
            return None

    def _add_agbnp2_ct(self,gb):
        """
        adds connection table information to AGBNP3 force
        """
        q = """SELECT p0,p1 FROM bond"""
        for (fcounter,conn,tables,offset) in self._localVars():
            for p0, p1 in conn.execute(q):
                gb.addParticleConnection(p0+offset,p1+offset)
        
    def createSystem(self, nonbondedMethod='CutoffPeriodic', nonbondedCutoff=1.2,
                     ewaldErrorTolerance=0.0005, removeCMMotion=True, hydrogenMass=None,
                     OPLS=False, implicitSolvent=None, AGBNPVersion=1, REM=False, A=100, C=50, martini=False,
                     improper_prefactor=0.99, tapering='shift', addForces=[]):
        """Construct an OpenMM System representing the topology described by this
        DMS file. tapering='shift' is must for Martini simulations because
        openMM uses MCBarostat, and for LJ-dominant systems like Martini 
        (rather than charged interactions dominant), MCBarostat will not work well.

        Parameters
        ----------
        nonbondedMethod : object=NoCutoff
            The method to use for nonbonded interactions.  Allowed values are
            NoCutoff, CutoffNonPeriodic, CutoffPeriodic, Ewald, PME, or LJPME.
        nonbondedCutoff : distance=1.2*nanometer
            The cutoff distance to use for nonbonded interactions
        ewaldErrorTolerance : float=0.0005
            The error tolerance to use if nonbondedMethod is Ewald, PME, or LJPME.
        removeCMMotion : boolean=True
            If true, a CMMotionRemover will be added to the System
        hydrogenMass : mass=None
            The mass to use for hydrogen atoms bound to heavy atoms.  Any mass
            added to a hydrogen is subtracted from the heavy atom to keep their
            total mass the same.
        OPLS : boolean=False
            If True, forces OPLS combining rules
        implicitSolvent: string=None
            If not None, creates implicit solvent force of the given name
            Allowed values are: HCT and 'AGBNP'
            (the corresponding tables must be present in the DMS file)
        AGBNPVersion: int=1
            AGBNP implicit solvent version
        REM: bool=False
            Turn on Reduced Nonbonded Energy Minimization (REM)
        A: float=100
            Parameter for repulsion in REM (kJ/mol)
        C: float=50
            Parameter for charged interaction in REM (kJ/mol)
        martini: bool=False
            Martini system (create CustomNonbondedForce, containing LJ+charged with a hard cutoff)
        improper_prefactor: float=0.99
            Prefactor for improper dihedrals (if 1.0, no problem in Mac but problem in Linux...)
        """
        self._checkForUnsupportedTerms()
        sys = mm.System()

        # Build the box dimensions
        boxSize = self.topology.getUnitCellDimensions()
        print('DMS boxSize:', boxSize)
        if boxSize is not None:
            if boxSize[0]._value * boxSize[1]._value * boxSize[2]._value == 0:
                print('DMS boxSize is set for a unit cell')
                sys.setDefaultPeriodicBoxVectors((1.0, 0, 0), (0, 1.0, 0), (0, 0, 1.0))
                self.cell = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
            else:
                sys.setDefaultPeriodicBoxVectors((boxSize[0], 0, 0), (0, boxSize[1], 0), (0, 0, boxSize[2]))
        elif nonbondedMethod in ('CutoffPeriodic', 'Ewald', 'PME', 'LJPME'):
            raise ValueError('Illegal nonbonded method for a non-periodic system')

        # Create all of the particles
        for (fcounter,conn,tables,offset) in self._localVars():
            for mass in conn.execute('SELECT mass FROM particle ORDER BY id'):
                sys.addParticle(mass[0]*dalton)

        # Add all of the forces
        self._addBondsToSystem(sys)
        self._addAnglesToSystem(sys)
        self._addCosineAnglesToSystem(sys)
        self._addConstraintsToSystem(sys)
        #self._addPeriodicTorsionsToSystem(sys, OPLS)
        self._addPeriodicTorsionsToSystemModified(sys)
        #self._addImproperHarmonicTorsionsToSystem(sys)
        self._addImproperHarmonicTorsionsToSystemModified(sys, improper_prefactor)
        self._addCMAPToSystem(sys)
        self._addVirtualSitesToSystem(sys)
        self._addPositionalHarmonicRestraints(sys)
        for addForce in addForces:
            sys.addForce(addForce)

        if REM and martini:
            # REM martini simulation
            self._addNonbondedForceToSystemTableMartiniREM(sys, nonbondedCutoff, nonbondedMethod, A, C, tapering)

        elif not REM and martini:
            # Normal martini simulations
            self._addNonbondedForceToSystemTableMartini(sys, nonbondedCutoff, nonbondedMethod, tapering)

        elif REM and not martini:
            # REM all-atom simulations
            self._addNonbondedForceToSystemTableREM(sys, nonbondedCutoff, nonbondedMethod, A, C)

        else:
            # Normal all-atom simulations
            self._addNonbondedForceToSystemTable(sys, nonbondedCutoff, nonbondedMethod, ewaldErrorTolerance)


        #add implicit solvent model.
        if implicitSolvent is not None:

            if not (implicitSolvent in (HCT, 'AGBNP', 'GVolSA', 'AGBNP3')):
                raise ValueError('Illegal implicit solvent method')
            
            if self._verbose:
                print('Adding implicit solvent ...')

            #with implicit solvent turn off native reaction field
            #However note that this does not affect the shifted Coulomb potential of the Nonbonded force
            #(it affects the only the energy, not the forces and equation of motion)
            nb.setReactionFieldDielectric(1.0)

            if implicitSolvent is HCT:
                gb_parms = self._get_gb_params()

                if gb_parms:
                    if self._verbose:
                        print('Adding HCT GB force ...')
                    gb = GBSAHCTForce(SA='ACE')
                    for i in range(len(gb_parms)):
                        gb.addParticle(list(gb_parms[i]))
                    gb.finalize()
                    sys.addForce(gb)
                else:
                    raise IOError("No HCT parameters found in DMS file")
                    
            if implicitSolvent == 'AGBNP3':
                #load AGBNP3 plugin if available
                try:
                    from AGBNP3plugin import AGBNP3Force
                except ImportError:
                    raise NotImplementedError('AGBNP3 is not supported in this version')
                #sets up AGBNP3
                gb_parms = self._get_agbnp2_params()
                if gb_parms:
                    if self._verbose:
                        print('Adding AGBNP3 force ...')
                    gb = AGBNP3Force()
                    # add particles
                    for i in range(len(gb_parms)):
                        p = gb_parms[i]
                        gb.addParticle(p[0],p[1],p[2],p[3],p[4],p[5],p[6])
                    # connection table (from bonds)
                    self._add_agbnp2_ct(gb)
                    sys.addForce(gb)
                else:
                    raise IOError("No AGBNP parameters found in DMS file")

            if implicitSolvent == 'GVolSA':
                #implemented as AGBNP version 0
                implicitSolvent = 'AGBNP'
                AGBNPVersion = 0
                if self._verbose:
                    print('Using GVolSA')

            if implicitSolvent == 'AGBNP':
                #load AGBNP plugin if available
                try:
                    from AGBNPplugin import AGBNPForce
                except ImportError:
                    raise NotImplementedError('AGBNP is not supported in this version')
                #sets up AGBNP
                gb_parms = self._get_agbnp2_params()
                if gb_parms:
                    gb = AGBNPForce()
                    gb.setNonbondedMethod(methodMap[nonbondedMethod])
                    gb.setCutoffDistance(nonbondedCutoff)
                    gb.setVersion(AGBNPVersion)
                    if self._verbose:
                        print('Using AGBNP force version %d ...' % AGBNPVersion)
                    # add particles
                    for i in range(len(gb_parms)):
                        [radiusN,chargeN,gammaN,alphaN,hbtype,hbwN,ishydrogenN] = gb_parms[i]
                        h_flag = ishydrogenN > 0
                        gb.addParticle(radiusN, gammaN, alphaN, chargeN, h_flag)
                    sys.addForce(gb)
                    self.gb_parms = gb_parms
                    self.agbnp = gb
                else:
                    raise IOError("No AGBNP parameters found in DMS file")

                    
        # Adjust masses.
        if hydrogenMass is not None:
            for atom1, atom2 in self.topology.bonds():

                if atom1.element == hydrogen:
                    (atom1, atom2) = (atom2, atom1)
                if atom2.element == hydrogen and atom1.element not in (hydrogen, None):
                    transferMass = hydrogenMass-sys.getParticleMass(atom2.index)
                    sys.setParticleMass(atom2.index, hydrogenMass)
                    sys.setParticleMass(atom1.index, sys.getParticleMass(atom1.index)-transferMass)

        # Add a CMMotionRemover.
        if removeCMMotion:
            sys.addForce(mm.CMMotionRemover())
        
        # Added
        self.system = sys
        return sys

    def _addBondsToSystem(self, sys):
        """Create the harmonic bonds
        """
        bonds = mm.HarmonicBondForce()
        sys.addForce(bonds)

        q = """SELECT p0, p1, r0, fc, constrained
        FROM stretch_harm_term INNER JOIN stretch_harm_param
        ON stretch_harm_term.param=stretch_harm_param.id"""
        for (fcounter,conn,tables,offset) in self._localVars():
            for p0, p1, r0, fc, constrained in conn.execute(q):
                p0 += offset
                p1 += offset
                if constrained:
                    sys.addConstraint(p0, p1, r0*angstrom)
                else:
                    # Desmond writes the harmonic bond force without 1/2
                    # so we need to to double the force constant
                    bonds.addBond(p0, p1, r0*angstrom, 2*fc*kilocalorie_per_mole/angstrom**2)

                # Record information that will be needed for constraining angles.
                self._atomBonds[p0][p1] = r0*angstrom
                self._atomBonds[p1][p0] = r0*angstrom

        return bonds

    def _addAnglesToSystem(self, sys):
        """Create the harmonic angles
        """
        angles = mm.HarmonicAngleForce()
        sys.addForce(angles)
        degToRad = math.pi/180

        q = """SELECT p0, p1, p2, theta0, fc, constrained
        FROM angle_harm_term INNER JOIN angle_harm_param
        ON angle_harm_term.param=angle_harm_param.id"""
        for (fcounter,conn,tables,offset) in self._localVars():
            for p0, p1, p2, theta0, fc, constrained in conn.execute(q):
                p0 += offset
                p1 += offset
                p2 += offset
                if constrained:
                    l1 = self._atomBonds[p1][p0]
                    l2 = self._atomBonds[p1][p2]
                    length = (l1*l1 + l2*l2 - 2*l1*l2*math.cos(theta0*degToRad)).sqrt()
                    sys.addConstraint(p0, p2, length)
                    self._angleConstraints[p1][p0] = p2
                    self._angleConstraints[p1][p2] = p0
                else:
                    # Desmond writes the harmonic angle force without 1/2
                    # so we need to to double the force constant
                    angles.addAngle(p0, p1, p2, theta0*degToRad, 2*fc*kilocalorie_per_mole/radian**2)

        return angles

    def _addCosineAnglesToSystem(self, sys):
        """Create the cosine angles for Martini simulations
        """
        if 'angle_harmcos_param' in self._tables[0].keys():
            angles = mm.CustomAngleForce("k * (cos(theta) - cos_theta0)^2;")
            sys.addForce(angles)
            angles.addPerAngleParameter("cos_theta0");
            angles.addPerAngleParameter("k");
            
            # need to add if table exists
            q = """SELECT p0, p1, p2, cos_theta0, fc 
            FROM angle_harmcos_term INNER JOIN angle_harmcos_param 
            ON angle_harmcos_term.param=angle_harmcos_param.id"""
            for (fcounter,conn,tables,offset) in self._localVars():
                for p0, p1, p2, cos_theta0, fc in conn.execute(q):
                    p0 += offset
                    p1 += offset
                    p2 += offset
                    angles.addAngle(p0, p1, p2, (cos_theta0, fc * kilocalorie_per_mole))

            return angles
 

    def _addConstraintsToSystem(self, sys):
        """Add constraints to system. Normally these should already be
        added by the bonds table, but we want to make sure that there's
        no extra information in the constraints table that we're not
        including in the system"""

        for (fcounter,conn,tables,offset) in self._localVars():
            for term_table in [n for n in list(tables.keys()) if n.startswith('constraint_a') and n.endswith('term')]:
                param_table = term_table.replace('term', 'param')
                q = """SELECT p0, p1, r1
                FROM %(term)s INNER JOIN %(param)s
                ON %(term)s.param=%(param)s.id""" % \
                    {'term': term_table, 'param': param_table}
                for p0, p1, r1 in conn.execute(q):
                    p0 += offset
                    p1 += offset
                    if not p1 in self._atomBonds[p0]:
                        sys.addConstraint(p0, p1, r1*angstrom)
                        self._atomBonds[p0][p1] = r1*angstrom
                        self._atomBonds[p1][p0] = r1*angstrom

            if 'constraint_hoh_term' in tables:
                degToRad = math.pi/180
                q = """SELECT p0, p1, p2, r1, r2, theta
                FROM constraint_hoh_term INNER JOIN constraint_hoh_param
                ON constraint_hoh_term.param=constraint_hoh_param.id"""
                for p0, p1, p2, r1, r2, theta in conn.execute(q):
                    p0 += offset
                    p1 += offset
                    p2 += offset
                    # Here, p0 is the heavy atom and p1 and p2 are the H1 and H2
                    # wihth O-H1 and O-H2 distances r1 and r2
                    if not (self._angleConstraints[p0].get(p1, None) == p2):
                        length = (r1*r1 + r2*r2 - 2*r1*r2*math.cos(theta*degToRad)).sqrt()
                        sys.addConstraint(p1, p2, length)

    def _addPeriodicTorsionsToSystem(self, sys, OPLS):
        """Create the torsion terms
        """
        if OPLS:
            periodic = mm.CustomTorsionForce('f * cos(n * theta - phi0)')
            periodic.addPerTorsionParameter('n')
            periodic.addPerTorsionParameter('phi0')
            periodic.addPerTorsionParameter('f')
        else:
            periodic = mm.PeriodicTorsionForce()
        sys.addForce(periodic)

        q = """SELECT p0, p1, p2, p3, phi0, fc0, fc1, fc2, fc3, fc4, fc5, fc6
        FROM dihedral_trig_term INNER JOIN dihedral_trig_param
        ON dihedral_trig_term.param=dihedral_trig_param.id"""

        for (fcounter,conn,tables,offset) in self._localVars():
            for p0, p1, p2, p3, phi0, fc0, fc1, fc2, fc3, fc4, fc5, fc6 in conn.execute(q):
                p0 += offset
                p1 += offset
                p2 += offset
                p3 += offset
                for order, fc in enumerate([fc0, fc1, fc2, fc3, fc4, fc5, fc6]):
                    if fc == 0:
                        continue
                    if OPLS:
                        periodic.addTorsion(p0, p1, p2, p3, [order, phi0*degree, fc*kilocalorie_per_mole])
                    else:
                        periodic.addTorsion(p0, p1, p2, p3, order, phi0*degree, fc*kilocalorie_per_mole)

    def _addPeriodicTorsionsToSystemModified(self, sys):
        """Create the torsion terms
        I think openmm has an incorrect converter for PeriodicTorsionsToSystem
        openMM:  sum_{n=0-6} (kn * (1 + cos(n * theta - theta0)))
        correct: k0 + sum_{n=1-6} (kn * cos(n * theta - theta0))
        """
        
        energy = 'fc0'
        for n in range(1, 7): energy += f' + fc{n} * cos({n} * theta - phi0)'
        periodic = mm.CustomTorsionForce(energy)
        periodic.addPerTorsionParameter('phi0')
        for n in range(0, 7): periodic.addPerTorsionParameter(f'fc{n}')
        sys.addForce(periodic)

        q = """SELECT p0, p1, p2, p3, phi0, fc0, fc1, fc2, fc3, fc4, fc5, fc6
        FROM dihedral_trig_term INNER JOIN dihedral_trig_param
        ON dihedral_trig_term.param=dihedral_trig_param.id"""

        for (fcounter,conn,tables,offset) in self._localVars():
            for p0, p1, p2, p3, phi0, fc0, fc1, fc2, fc3, fc4, fc5, fc6 in conn.execute(q):
                p0 += offset
                p1 += offset
                p2 += offset
                p3 += offset
                phi0 *= degree
                fc0 *= kilocalorie_per_mole
                fc1 *= kilocalorie_per_mole
                fc2 *= kilocalorie_per_mole
                fc3 *= kilocalorie_per_mole
                fc4 *= kilocalorie_per_mole
                fc5 *= kilocalorie_per_mole
                fc6 *= kilocalorie_per_mole

                periodic.addTorsion(p0, p1, p2, p3, (phi0, fc0, fc1, fc2, fc3, fc4, fc5, fc6))


    def _addImproperHarmonicTorsionsToSystem(self, sys):
        """Create the improper harmonic torsion terms
        """
        harmonicTorsion = mm.CustomTorsionForce('k*(theta-theta0)^2')
        harmonicTorsion.addPerTorsionParameter('theta0')
        harmonicTorsion.addPerTorsionParameter('k')

        go = []

        for (fcounter,conn,tables,offset) in self._localVars():
            if not self._hasTable('improper_harm_term', tables):
                go.append(False)
            else:
                go.append(True)

        if any(go):
            sys.addForce(harmonicTorsion)
        else:
            return

        q = """SELECT p0, p1, p2, p3, phi0, fc
        FROM improper_harm_term INNER JOIN improper_harm_param
        ON improper_harm_term.param=improper_harm_param.id"""

        for (fcounter,conn,tables,offset) in self._localVars():
            if not go[fcounter]:
                continue
            for p0, p1, p2, p3, phi0, fc in conn.execute(q):
                p0 += offset
                p1 += offset
                p2 += offset
                p3 += offset
                harmonicTorsion.addTorsion(p0, p1, p2, p3, [phi0*degree, fc*kilocalorie_per_mole])

    def _addImproperHarmonicTorsionsToSystemModified(self, sys, improper_prefactor=0.99):
        """Create the improper harmonic torsion terms
        periodicity at phi0 +- 180 is not taken into account with k * (theta-theta0)^2
        """
        harmonicTorsion = mm.CustomTorsionForce(f'k * (acos({improper_prefactor} * cos(theta-theta0)))^2')
        harmonicTorsion.addPerTorsionParameter('theta0')
        harmonicTorsion.addPerTorsionParameter('k')

        go = []

        for (fcounter,conn,tables,offset) in self._localVars():
            if not self._hasTable('improper_harm_term', tables):
                go.append(False)
            else:
                go.append(True)

        if any(go):
            sys.addForce(harmonicTorsion)
        else:
            return

        q = """SELECT p0, p1, p2, p3, phi0, fc
        FROM improper_harm_term INNER JOIN improper_harm_param
        ON improper_harm_term.param=improper_harm_param.id"""

        for (fcounter,conn,tables,offset) in self._localVars():
            if not go[fcounter]:
                continue
            for p0, p1, p2, p3, phi0, fc in conn.execute(q):
                p0 += offset
                p1 += offset
                p2 += offset
                p3 += offset
                harmonicTorsion.addTorsion(p0, p1, p2, p3, [phi0*degree, fc*kilocalorie_per_mole])

    def _addCMAPToSystem(self, sys):
        """Create the CMAP terms
        """
        go = []

        for (fcounter,conn,tables,offset) in self._localVars():
            if not self._hasTable('torsiontorsion_cmap_term', tables):
                go.append(False)
            else:
                go.append(True)

        if any(go):
            # Create CMAP torsion terms
            cmap = mm.CMAPTorsionForce()
            sys.addForce(cmap)
        else:
            return

        cmap_indices = {}

        for (fcounter,conn,tables,offset) in self._localVars():
            if not go[fcounter]:
                continue
            for name in [k for k in list(tables.keys()) if k.startswith('cmap')]:
                size2 = conn.execute('SELECT COUNT(*) FROM %s' % name).fetchone()[0]
                fsize = math.sqrt(size2)
                if fsize != int(fsize):
                    raise ValueError('Non-square CMAPs are not supported')
                size = int(fsize)

                map = [0 for i in range(size2)]
                for phi, psi, energy in conn.execute("SELECT phi, psi, energy FROM %s" % name):
                    i = int((phi % 360) / (360.0 / size))
                    j = int((psi % 360) / (360.0 / size))
                    map[i+size*j] = energy
                index = cmap.addMap(size, map*kilocalorie_per_mole)
                cmap_indices[name] = index

            q = """SELECT p0, p1, p2, p3, p4, p5, p6, p7, cmapid
            FROM torsiontorsion_cmap_term INNER JOIN torsiontorsion_cmap_param
            ON torsiontorsion_cmap_term.param=torsiontorsion_cmap_param.id"""
            for p0, p1, p2, p3, p4, p5, p6, p7, cmapid in conn.execute(q):
                p0 += offset
                p1 += offset
                p2 += offset
                p3 += offset
                p4 += offset
                p5 += offset
                p6 += offset
                p7 += offset
                cmap.addTorsion(cmap_indices[cmapid], p0, p1, p2, p3, p4, p5, p6, p7)


    def _addNonbondedForceToSystemTable(self, sys, 
                                        nonbondedCutoff=1.2, 
                                        nonbondedMethod='CutoffPeriodic',
                                        ewaldErrorTolerance=0.0005):
        """Create the nonbonded force using TabulatedFunction"""
        
        energy = ("4*eps*((sig/r)^12-(sig/r)^6); "
                  "eps = epsilon(type1, type2); "
                  "sig = sigma(type1, type2); ")
        name   = 'LennardJones'
        cnb = mm.CustomNonbondedForce(energy)
        cnb.setNonbondedMethod(methodMap[nonbondedMethod])
        cnb.setCutoffDistance(nonbondedCutoff)
        cnb.setUseSwitchingFunction(False)
        cnb.setUseLongRangeCorrection(False)
 
        ### Get sigma and epsilon for each nbtype
        q = """SELECT id, sigma, epsilon FROM nonbonded_param ORDER BY id"""
        for (fcounter,conn,tables,offset) in self._localVars():
            numLjTypes = len(conn.execute(q).fetchall())
            sigma2D   = np.zeros((numLjTypes, numLjTypes))
            epsilon2D = np.zeros((numLjTypes, numLjTypes))
            
            for id, sigma, epsilon in conn.execute(q):
                sigma2D[id, id]   = sigma
                epsilon2D[id, id] = epsilon
        
        ### Complete 2D matrix
        for i in range(numLjTypes):
            for j in range(numLjTypes):
                if i == j: continue
                sigma2D[i, j]   = 0.5 *   ( sigma2D[i, i]   + sigma2D[j, j] )
                epsilon2D[i, j] = np.sqrt ( epsilon2D[i, i] * epsilon2D[j, j] )

        ### NBFIX
        try:
            q = """SELECT param1, param2, sigma, epsilon FROM nonbonded_combined_param"""
            for (fcounter,conn,tables,offset) in self._localVars():
                for param1, param2, sigma, epsilon in conn.execute(q):
                    sigma2D[param1, param2]   = sigma
                    sigma2D[param2, param1]   = sigma
                    epsilon2D[param1, param2] = epsilon
                    epsilon2D[param2, param1] = epsilon
        except:
            pass

        sigma2D   *= 0.1
        epsilon2D *= 4.184

        cnb.addTabulatedFunction('sigma',   mm.Discrete2DFunction(numLjTypes, numLjTypes, sigma2D.flatten()))
        cnb.addTabulatedFunction('epsilon', mm.Discrete2DFunction(numLjTypes, numLjTypes, epsilon2D.flatten()))
        cnb.addPerParticleParameter('type')

        ### Nonbonded is needed because of charged interactions
        nb = mm.NonbondedForce()
        nb.setNonbondedMethod(methodMap[nonbondedMethod])
        nb.setCutoffDistance(nonbondedCutoff)
        nb.setEwaldErrorTolerance(ewaldErrorTolerance)
        nb.setUseSwitchingFunction(False)
        nb.setUseDispersionCorrection(False)

        ### Add charge to nb and add nbtype to cnb for each particle
        q = """SELECT charge, nbtype FROM particle ORDER BY id"""
        for (fcounter,conn,tables,offset) in self._localVars():
            for charge, nbtype in conn.execute(q):
                nb.addParticle(charge, 1.0*angstrom, 0.0*kilocalorie_per_mole)
                cnb.addParticle((nbtype,))

        ### addException to nb and addExclusion to cnb
        for (fcounter,conn,tables,offset) in self._localVars():
            for p0, p1 in conn.execute('SELECT p0, p1 FROM exclusion'):
                p0 += offset
                p1 += offset
                nb.addException(p0, p1, 0.0, 1.0, 0.0)
                cnb.addExclusion(p0, p1)

        ### 1-4 interactions. No need to make CustomBondForce if you use NonbondedForce.addException
        q = """SELECT p0, p1, aij, bij, qij
        FROM pair_12_6_es_term INNER JOIN pair_12_6_es_param
        ON pair_12_6_es_term.param=pair_12_6_es_param.id"""
        for (fcounter,conn,tables,offset) in self._localVars():
            for p0, p1, a_ij, b_ij, q_ij in conn.execute(q):
                p0 += offset
                p1 += offset
                a_ij = (a_ij*kilocalorie_per_mole*(angstrom**12)).in_units_of(kilojoule_per_mole*(nanometer**12))
                b_ij = (b_ij*kilocalorie_per_mole*(angstrom**6)).in_units_of(kilojoule_per_mole*(nanometer**6))
                q_ij = q_ij*elementary_charge**2
                if (b_ij._value == 0.0) or (a_ij._value == 0.0):
                    new_epsilon = 0
                    new_sigma = 1
                else:
                    new_epsilon =  b_ij**2/(4*a_ij)
                    new_sigma = (a_ij / b_ij)**(1.0/6.0)
                nb.addException(p0, p1, q_ij, new_sigma, new_epsilon, True)

            n_total = conn.execute("""SELECT COUNT(*) FROM pair_12_6_es_term""").fetchone()
            n_in_exclusions = conn.execute("""SELECT COUNT(*)
            FROM exclusion INNER JOIN pair_12_6_es_term
            ON (    ( exclusion.p0==pair_12_6_es_term.p0 AND exclusion.p1==pair_12_6_es_term.p1)
                 OR ( exclusion.p0==pair_12_6_es_term.p1 AND exclusion.p1==pair_12_6_es_term.p0) 
               )""").fetchone()
            if not n_total == n_in_exclusions:
                raise NotImplementedError('All pair_12_6_es_terms must have a corresponding exclusion')
        
        sys.addForce(nb)
        sys.addForce(cnb)


    def _addNonbondedForceToSystemTableREM(self, sys, 
                                           nonbondedCutoff=1.2, 
                                           nonbondedMethod='CutoffPeriodic',
                                           A=100, C=50):
        """Create the REM nonbonded force using TabulatedFunction"""

        energy = ("min(rep, LJ) + coul; "
                  "rep  = A * (cos(pi/2 * r/sig))^2; "
                  "LJ   = 4 * eps * ((sig/r)^12-(sig/r)^6); "
                  "coul = C * q1 * q2 * (cos(pi/2 * r / rcut))^2; "
                  "eps  = epsilon(type1, type2); "
                  "sig=sigma(type1, type2);")
        
        cnb = mm.CustomNonbondedForce(energy)
        cnb.setNonbondedMethod(methodMap[nonbondedMethod])
        cnb.setCutoffDistance(nonbondedCutoff)
        cnb.setUseSwitchingFunction(False)
        cnb.setUseLongRangeCorrection(False)
        cnb.addGlobalParameter('pi',   3.141592)
        cnb.addGlobalParameter('A',    A)
        cnb.addGlobalParameter('C',    C)
        cnb.addGlobalParameter('rcut', nonbondedCutoff)

 
        ### Get sigma and epsilon for each nbtype
        q = """SELECT id, sigma, epsilon FROM nonbonded_param ORDER BY id"""
        for (fcounter,conn,tables,offset) in self._localVars():
            numLjTypes = len(conn.execute(q).fetchall())
            sigma2D   = np.zeros((numLjTypes, numLjTypes))
            epsilon2D = np.zeros((numLjTypes, numLjTypes))
            
            for id, sigma, epsilon in conn.execute(q):
                sigma2D[id, id]   = sigma
                epsilon2D[id, id] = epsilon
        
        ### Complete 2D matrix
        for i in range(numLjTypes):
            for j in range(numLjTypes):
                if i == j: continue
                sigma2D[i, j]   = 0.5 *  ( sigma2D[i, i]   + sigma2D[j, j] )
                epsilon2D[i, j] = np.sqrt( epsilon2D[i, i] * epsilon2D[j, j] )

        ### NBFIX
        try:
            q = """SELECT param1, param2, sigma, epsilon FROM nonbonded_combined_param"""
            for (fcounter,conn,tables,offset) in self._localVars():
                for param1, param2, sigma, epsilon in conn.execute(q):
                    sigma2D[param1, param2]   = sigma
                    sigma2D[param2, param1]   = sigma
                    epsilon2D[param1, param2] = epsilon
                    epsilon2D[param2, param1] = epsilon

        except:
            pass

        sigma2D   *= 0.1
        epsilon2D *= 4.184

        cnb.addTabulatedFunction('sigma',   mm.Discrete2DFunction(numLjTypes, numLjTypes, sigma2D.flatten()))
        cnb.addTabulatedFunction('epsilon', mm.Discrete2DFunction(numLjTypes, numLjTypes, epsilon2D.flatten()))
        cnb.addPerParticleParameter('q')
        cnb.addPerParticleParameter('type')

        ### Add charge to nb and add nbtype to cnb for each particle
        q = """SELECT charge, nbtype FROM particle ORDER BY id"""
        for (fcounter,conn,tables,offset) in self._localVars():
            for charge, nbtype in conn.execute(q):
                cnb.addParticle((charge, nbtype))

        ### addExclusion to cnb
        for (fcounter,conn,tables,offset) in self._localVars():
            for p0, p1 in conn.execute('SELECT p0, p1 FROM exclusion'):
                p0 += offset
                p1 += offset
                cnb.addExclusion(p0, p1)

        # New Custom Bonded Force (1-4 interactions)
        cb = mm.CustomBondForce("min(rep, LJ) + coul; \
            rep  = A * (cos(pi/2 * r/sigma))^2; \
            LJ   = 4 * epsilon * ((sigma/r)^12-(sigma/r)^6); \
            coul = C * q1q2 * (cos(pi/2 * r / rcut))^2;")
        cb.addPerBondParameter("sigma")
        cb.addPerBondParameter("epsilon")
        cb.addPerBondParameter('q1q2')
        cb.addGlobalParameter('pi',   3.141592)
        cb.addGlobalParameter('A',    A)
        cb.addGlobalParameter('C',    C)
        cb.addGlobalParameter('rcut', nonbondedCutoff)

        q = """SELECT p0, p1, aij, bij, qij
        FROM pair_12_6_es_term INNER JOIN pair_12_6_es_param
        ON pair_12_6_es_term.param=pair_12_6_es_param.id"""
        for (fcounter,conn,tables,offset) in self._localVars():
            for p0, p1, a_ij, b_ij, q_ij in conn.execute(q):
                p0 += offset
                p1 += offset
                a_ij = (a_ij*kilocalorie_per_mole*(angstrom**12)).in_units_of(kilojoule_per_mole*(nanometer**12))
                b_ij = (b_ij*kilocalorie_per_mole*(angstrom**6)).in_units_of(kilojoule_per_mole*(nanometer**6))
                q_ij = q_ij*elementary_charge**2
                if (b_ij._value == 0.0) or (a_ij._value == 0.0):
                    new_epsilon = 0
                    new_sigma = 1
                else:
                    new_epsilon =  b_ij**2/(4*a_ij)
                    new_sigma = (a_ij / b_ij)**(1.0/6.0)
                cb.addBond(p0, p1, [new_sigma, new_epsilon, q_ij])

            n_total = conn.execute("""SELECT COUNT(*) FROM pair_12_6_es_term""").fetchone()
            n_in_exclusions = conn.execute("""SELECT COUNT(*)
            FROM exclusion INNER JOIN pair_12_6_es_term
            ON (    ( exclusion.p0==pair_12_6_es_term.p0 AND exclusion.p1==pair_12_6_es_term.p1)
                 OR ( exclusion.p0==pair_12_6_es_term.p1 AND exclusion.p1==pair_12_6_es_term.p0) 
               )""").fetchone()
            if not n_total == n_in_exclusions:
                raise NotImplementedError('All pair_12_6_es_terms must have a corresponding exclusion')

        sys.addForce(cb)
        sys.addForce(cnb)


    def _addNonbondedForceToSystemTableMartini(self, sys, 
                                               nonbondedCutoff=1.1, 
                                               nonbondedMethod='CutoffPeriodic',
                                               tapering='shift'):
        """Create the nonbonded force using TabulatedFunction for Martini"""
        
        if tapering:
            if tapering == 'shift':
                energy = ("4*eps*((sig/r)^12-(sig/r)^6) + 138.911*q1*q2/r "
                         f"- 4*eps*((sig/{nonbondedCutoff})^12-(sig/{nonbondedCutoff})^6) "
                         f"- 138.911*q1*q2/{nonbondedCutoff}; "
                          "eps = epsilon(type1, type2); "
                          "sig = sigma(type1, type2); ")
            else:
                raise ValueError('tapering should be either "shift" or False')

        else:
            energy = ("4*eps*((sig/r)^12-(sig/r)^6) + 138.911*q1*q2/r; "
                      "eps = epsilon(type1, type2); "
                      "sig = sigma(type1, type2); ")
        
        cnb = mm.CustomNonbondedForce(energy)
        cnb.setNonbondedMethod(methodMap[nonbondedMethod])
        cnb.setCutoffDistance(nonbondedCutoff)
        cnb.setUseSwitchingFunction(False)
        cnb.setUseLongRangeCorrection(False)
 
        ### Get sigma and epsilon for each nbtype
        q = """SELECT id, sigma, epsilon FROM nonbonded_param ORDER BY id"""
        for (fcounter,conn,tables,offset) in self._localVars():
            numLjTypes = len(conn.execute(q).fetchall())
            sigma2D   = np.zeros((numLjTypes, numLjTypes))
            epsilon2D = np.zeros((numLjTypes, numLjTypes))
            
            for id, sigma, epsilon in conn.execute(q):
                sigma2D[id, id]   = sigma
                epsilon2D[id, id] = epsilon
        
        ### Complete 2D matrix
        for i in range(numLjTypes):
            for j in range(numLjTypes):
                if i == j: continue
                sigma2D[i, j]   = 0.5 *   ( sigma2D[i, i]   + sigma2D[j, j] )
                epsilon2D[i, j] = np.sqrt ( epsilon2D[i, i] * epsilon2D[j, j] )

        ### NBFIX
        q = """SELECT param1, param2, sigma, epsilon FROM nonbonded_combined_param"""
        for (fcounter,conn,tables,offset) in self._localVars():
            for param1, param2, sigma, epsilon in conn.execute(q):
                sigma2D[param1, param2]   = sigma
                sigma2D[param2, param1]   = sigma
                epsilon2D[param1, param2] = epsilon
                epsilon2D[param2, param1] = epsilon

        sigma2D   *= 0.1
        epsilon2D *= 4.184

        cnb.addTabulatedFunction('sigma',   mm.Discrete2DFunction(numLjTypes, numLjTypes, sigma2D.flatten()))
        cnb.addTabulatedFunction('epsilon', mm.Discrete2DFunction(numLjTypes, numLjTypes, epsilon2D.flatten()))
        cnb.addPerParticleParameter('q')
        cnb.addPerParticleParameter('type')

        ### Add charge and nbtype to cnb for each particle
        q = """SELECT charge, nbtype FROM particle ORDER BY id"""
        for (fcounter,conn,tables,offset) in self._localVars():
            for charge, nbtype in conn.execute(q):
                cnb.addParticle((charge, nbtype))

        ### addExclusion to cnb
        for (fcounter,conn,tables,offset) in self._localVars():
            for p0, p1 in conn.execute('SELECT p0, p1 FROM exclusion'):
                p0 += offset
                p1 += offset
                cnb.addExclusion(p0, p1)

        sys.addForce(cnb)


    def _addNonbondedForceToSystemTableMartiniREM(self, sys, 
                                                  nonbondedCutoff=1.1, 
                                                  nonbondedMethod='CutoffPeriodic',
                                                  A=100, C=50, tapering='shift'):
        """Create the REM nonbonded force using TabulatedFunction for Martini"""

        if tapering:
            if tapering == 'shift':
                # Both rep and coul have pot = 0 at cutoff
                # adjust only LJ
                energy = ("min(rep, LJ) + coul - LJ0; "
                          "rep  = A * (cos(pi/2 * r/sig))^2; "
                          "LJ   = 4*eps*((sig/r)^12-(sig/r)^6); "
                          "coul = C * q1 * q2 * (cos(pi/2 * r / rcut))^2; "
                         f"LJ0  = 4*eps*((sig/{nonbondedCutoff})^12-(sig/{nonbondedCutoff})^6);"
                          "eps=epsilon(type1, type2); sig=sigma(type1, type2); ")


            else:
                raise ValueError('tapering should be either "shift" or False')

        else:
            energy = ("min(rep, LJ) + coul; "
                      "rep = A * (cos(pi/2 * r/sig))^2; "
                      "LJ = 4*eps*((sig/r)^12-(sig/r)^6); "
                      "coul = C * q1 * q2 * (cos(pi/2 * r / rcut))^2; "
                      "eps=epsilon(type1, type2); sig=sigma(type1, type2); ")
        
        cnb = mm.CustomNonbondedForce(energy)
        cnb.setNonbondedMethod(methodMap[nonbondedMethod])
        cnb.setCutoffDistance(nonbondedCutoff)
        cnb.setUseSwitchingFunction(False)
        cnb.setUseLongRangeCorrection(False)
        cnb.addGlobalParameter('pi',   3.141592)
        cnb.addGlobalParameter('A',    A)
        cnb.addGlobalParameter('C',    C)
        cnb.addGlobalParameter('rcut', nonbondedCutoff)

        ### Get sigma and epsilon for each nbtype
        q = """SELECT id, sigma, epsilon FROM nonbonded_param ORDER BY id"""
        for (fcounter,conn,tables,offset) in self._localVars():
            numLjTypes = len(conn.execute(q).fetchall())
            sigma2D   = np.zeros((numLjTypes, numLjTypes))
            epsilon2D = np.zeros((numLjTypes, numLjTypes))
            
            for id, sigma, epsilon in conn.execute(q):
                sigma2D[id, id]   = sigma
                epsilon2D[id, id] = epsilon
        
        ### Complete 2D matrix
        for i in range(numLjTypes):
            for j in range(numLjTypes):
                if i == j: continue
                sigma2D[i, j]   = 0.5 *   ( sigma2D[i, i]   + sigma2D[j, j] )
                epsilon2D[i, j] = np.sqrt ( epsilon2D[i, i] * epsilon2D[j, j] )

        ### NBFIX
        q = """SELECT param1, param2, sigma, epsilon FROM nonbonded_combined_param"""
        for (fcounter,conn,tables,offset) in self._localVars():
            for param1, param2, sigma, epsilon in conn.execute(q):
                sigma2D[param1, param2]   = sigma
                sigma2D[param2, param1]   = sigma
                epsilon2D[param1, param2] = epsilon
                epsilon2D[param2, param1] = epsilon

        sigma2D   *= 0.1
        epsilon2D *= 4.184

        cnb.addTabulatedFunction('sigma',   mm.Discrete2DFunction(numLjTypes, numLjTypes, sigma2D.flatten()))
        cnb.addTabulatedFunction('epsilon', mm.Discrete2DFunction(numLjTypes, numLjTypes, epsilon2D.flatten()))
        cnb.addPerParticleParameter('q')
        cnb.addPerParticleParameter('type')

        ### Add charge and nbtype to cnb for each particle
        q = """SELECT charge, nbtype FROM particle ORDER BY id"""
        for (fcounter,conn,tables,offset) in self._localVars():
            for charge, nbtype in conn.execute(q):
                cnb.addParticle((charge, nbtype))

        ### addExclusion to cnb
        for (fcounter,conn,tables,offset) in self._localVars():
            for p0, p1 in conn.execute('SELECT p0, p1 FROM exclusion'):
                p0 += offset
                p1 += offset
                cnb.addExclusion(p0, p1)

        sys.addForce(cnb)


    def _addNonbondedForceToSystem(self, sys, OPLS):
        """Create the nonbonded force
        """
        cnb = None
        nb = mm.NonbondedForce()
        sys.addForce(nb)

        if OPLS:
            cnb = mm.CustomNonbondedForce("4.0*epsilon12*((sigma12/r)^12 - (sigma12/r)^6); sigma12=sqrt(sigma1*sigma2); epsilon12=sqrt(epsilon1*epsilon2)")
            cnb.addPerParticleParameter("sigma")
            cnb.addPerParticleParameter("epsilon")
            sys.addForce(cnb)

        if OPLS:
            q = """SELECT sigma, epsilon
            FROM particle INNER JOIN nonbonded_param
            ON particle.nbtype=nonbonded_param.id ORDER BY particle.id"""
            for (fcounter,conn,tables,offset) in self._localVars():
                for sigma, epsilon in conn.execute(q):
                    cnb.addParticle([sigma*angstrom, epsilon*kilocalorie_per_mole])

        q = """SELECT charge, sigma, epsilon
        FROM particle INNER JOIN nonbonded_param
        ON particle.nbtype=nonbonded_param.id ORDER BY particle.id"""
        for (fcounter,conn,tables,offset) in self._localVars():
            for charge, sigma, epsilon in conn.execute(q):
                if OPLS:
                    epsilon = 0
                nb.addParticle(charge, sigma*angstrom, epsilon*kilocalorie_per_mole)

        for (fcounter,conn,tables,offset) in self._localVars():
            for p0, p1 in conn.execute('SELECT p0, p1 FROM exclusion'):
                p0 += offset
                p1 += offset
                nb.addException(p0, p1, 0.0, 1.0, 0.0)
                if OPLS:
                    cnb.addExclusion(p0, p1)

        q = """SELECT p0, p1, aij, bij, qij
        FROM pair_12_6_es_term INNER JOIN pair_12_6_es_param
        ON pair_12_6_es_term.param=pair_12_6_es_param.id"""
        for (fcounter,conn,tables,offset) in self._localVars():
            for p0, p1, a_ij, b_ij, q_ij in conn.execute(q):
                p0 += offset
                p1 += offset
                a_ij = (a_ij*kilocalorie_per_mole*(angstrom**12)).in_units_of(kilojoule_per_mole*(nanometer**12))
                b_ij = (b_ij*kilocalorie_per_mole*(angstrom**6)).in_units_of(kilojoule_per_mole*(nanometer**6))
                q_ij = q_ij*elementary_charge**2
                if (b_ij._value == 0.0) or (a_ij._value == 0.0):
                    new_epsilon = 0
                    new_sigma = 1
                else:
                    new_epsilon =  b_ij**2/(4*a_ij)
                    new_sigma = (a_ij / b_ij)**(1.0/6.0)
                nb.addException(p0, p1, q_ij, new_sigma, new_epsilon, True)

            n_total = conn.execute("""SELECT COUNT(*) FROM pair_12_6_es_term""").fetchone()
            n_in_exclusions = conn.execute("""SELECT COUNT(*)
            FROM exclusion INNER JOIN pair_12_6_es_term
            ON (    ( exclusion.p0==pair_12_6_es_term.p0 AND exclusion.p1==pair_12_6_es_term.p1)
                 OR ( exclusion.p0==pair_12_6_es_term.p1 AND exclusion.p1==pair_12_6_es_term.p0) 
               )""").fetchone()
            if not n_total == n_in_exclusions:
                raise NotImplementedError('All pair_12_6_es_terms must have a corresponding exclusion')

        return nb, cnb

    def _addNonbondedForceToSystemREM(self, sys, A, C, rcut):
        """Create the Reduced Nonbonded Energy Minimization (REM) nonbonded force
        NonbondedForce is not created. CustomNonbondedForce will take care of 
        both electrostatic and vdw interactions.
        """
        
        energy = ("min(rep, LJ) + coul; "
                  "rep = A * (cos(pi/2 * r/sigma))^2; "
                  "LJ = 4 * epsilon * ((sigma/r)^12-(sigma/r)^6); "
                  "coul = C * q1 * q2 * (cos(pi/2 * r / rcut))^2; "
                  "sigma = 0.5 * (sigma1 + sigma2); "
                  "epsilon = sqrt(epsilon1*epsilon2); ")

        # New Custom Nonbonded Force
        cnb = mm.CustomNonbondedForce(energy)

        cnb.addPerParticleParameter("sigma")
        cnb.addPerParticleParameter("epsilon")
        cnb.addPerParticleParameter('q')
        cnb.addGlobalParameter('pi',   3.141592)
        cnb.addGlobalParameter('A',    A)
        cnb.addGlobalParameter('C',    C)
        cnb.addGlobalParameter('rcut', rcut)
        cnb.setName('REM')
        cnb.setUseLongRangeCorrection(False)
        #cnb.setNonbondedMethod(mm.CutoffPeriodic)
        cnb.setCutoffDistance(rcut)
        cnb.setUseSwitchingFunction(False)
        sys.addForce(cnb)

        q = """SELECT charge, sigma, epsilon
        FROM particle INNER JOIN nonbonded_param
        ON particle.nbtype=nonbonded_param.id ORDER BY particle.id"""
        for (fcounter,conn,tables,offset) in self._localVars():
            for charge, sigma, epsilon in conn.execute(q):
                cnb.addParticle([sigma*0.1, epsilon*4.184, charge])

        for (fcounter,conn,tables,offset) in self._localVars():
            for p0, p1 in conn.execute('SELECT p0, p1 FROM exclusion'):
                p0 += offset
                p1 += offset
                cnb.addExclusion(p0, p1)


        # New Custom Bonded Force (1-4 interactions)
        cb = mm.CustomBondForce("min(rep, LJ) + coul; \
            rep  = A * (cos(pi/2 * r/sigma))^2; \
            LJ   = 4 * epsilon * ((sigma/r)^12-(sigma/r)^6); \
            coul = C * q1q2 * (cos(pi/2 * r / rcut))^2;")
        cb.addPerBondParameter("sigma")
        cb.addPerBondParameter("epsilon")
        cb.addPerBondParameter('q1q2')
        cb.addGlobalParameter('pi',   3.141592)
        cb.addGlobalParameter('A',    A)
        cb.addGlobalParameter('C',    C)
        cb.addGlobalParameter('rcut', rcut * nanometer)
        cb.setName('REM14')
        sys.addForce(cb)

        q = """SELECT p0, p1, aij, bij, qij
        FROM pair_12_6_es_term INNER JOIN pair_12_6_es_param
        ON pair_12_6_es_term.param=pair_12_6_es_param.id"""
        for (fcounter,conn,tables,offset) in self._localVars():
            for p0, p1, a_ij, b_ij, q_ij in conn.execute(q):
                p0 += offset
                p1 += offset
                a_ij = (a_ij*kilocalorie_per_mole*(angstrom**12)).in_units_of(kilojoule_per_mole*(nanometer**12))
                b_ij = (b_ij*kilocalorie_per_mole*(angstrom**6)).in_units_of(kilojoule_per_mole*(nanometer**6))
                q_ij = q_ij*elementary_charge**2
                if (b_ij._value == 0.0) or (a_ij._value == 0.0):
                    new_epsilon = 0
                    new_sigma = 1
                else:
                    new_epsilon =  b_ij**2/(4*a_ij)
                    new_sigma = (a_ij / b_ij)**(1.0/6.0)
                cb.addBond(p0, p1, [new_sigma, new_epsilon, q_ij])

            n_total = conn.execute("""SELECT COUNT(*) FROM pair_12_6_es_term""").fetchone()
            n_in_exclusions = conn.execute("""SELECT COUNT(*)
            FROM exclusion INNER JOIN pair_12_6_es_term
            ON (    ( exclusion.p0==pair_12_6_es_term.p0 AND exclusion.p1==pair_12_6_es_term.p1)
                 OR ( exclusion.p0==pair_12_6_es_term.p1 AND exclusion.p1==pair_12_6_es_term.p0) 
               )""").fetchone()
            if not n_total == n_in_exclusions:
                raise NotImplementedError('All pair_12_6_es_terms must have a corresponding exclusion')

        return cnb

    def _addVirtualSitesToSystem(self, sys):
        """Create any virtual sites in the system
        """
        go = []

        for (fcounter,conn,tables,offset) in self._localVars():
            if not any(t.startswith('virtual_') for t in list(tables.keys())):
                go.append(False)
            else:
                go.append(True)

        if not any(go):
            return

        for (fcounter,conn,tables,offset) in self._localVars():
            if not go[fcounter]:
                continue
            if 'virtual_lc2_term' in tables:
                q = """SELECT p0, p1, p2, c1
                FROM virtual_lc2_term INNER JOIN virtual_lc2_param
                ON virtual_lc2_term.param=virtual_lc2_param.id"""
                for p0, p1, p2, c1 in conn.execute(q):
                    p0 += offset
                    p1 += offset
                    p2 += offset
                    vsite = mm.TwoParticleAverageSite(p1, p2, (1-c1), c1)
                    sys.setVirtualSite(p0, vsite)

        for (fcounter,conn,tables,offset) in self._localVars():
            if not go[fcounter]:
                continue
            if 'virtual_lc3_term' in tables:
                q = """SELECT p0, p1, p2, p3, c1, c2
                FROM virtual_lc3_term INNER JOIN virtual_lc3_param
                ON virtual_lc3_term.param=virtual_lc3_param.id"""
                for p0, p1, p2, p3, c1, c2 in conn.execute(q):
                    p0 += offset
                    p1 += offset
                    p2 += offset
                    p3 += offset
                    vsite = mm.ThreeParticleAverageSite(p1, p2, p3, (1-c1-c2), c1, c2)
                    sys.setVirtualSite(p0, vsite)

        for (fcounter,conn,tables,offset) in self._localVars():
            if not go[fcounter]:
                continue
            if 'virtual_out3_term' in tables:
                q = """SELECT p0, p1, p2, p3, c1, c2, c3
                FROM virtual_out3_term INNER JOIN virtual_out3_param
                ON virtual_out3_term.param=virtual_out3_param.id"""
                for p0, p1, p2, p3, c1, c2, c3 in conn.execute(q):
                    p0 += offset
                    p1 += offset
                    p2 += offset
                    p3 += offset
                    vsite = mm.OutOfPlaneSite(p1, p2, p3, c1, c2, c3)
                    sys.setVirtualSite(p0, vsite)

        for (fcounter,conn,tables,offset) in self._localVars():
            if not go[fcounter]:
                continue
            if 'virtual_fdat3_term' in tables:
                raise NotImplementedError('OpenMM does not currently support '
                                          'fdat3-style virtual sites')

    def _addPositionalHarmonicRestraints(self, sys):

        go = []

        for (fcounter,conn,tables,offset) in self._localVars():
            if not self._hasTable('posre_harm_term',tables):
                go.append(False)
            else:
                go.append(True)
            if go[fcounter] and (not self._hasTable('posre_harm_param',tables)):
                raise IOError('DMS file lacks posre_harm_param table even though posre_harm_term table is present.')

        if not any(go):
            return

        if self._verbose:
            print("Using positional harmonic restraints.")

        force = mm.CustomExternalForce("hkx*(x-x0)^2+hky*(y-y0)^2+hkz*(z-z0)^2")
        force.addPerParticleParameter("x0")
        force.addPerParticleParameter("y0")
        force.addPerParticleParameter("z0")
        force.addPerParticleParameter("hkx")
        force.addPerParticleParameter("hky")
        force.addPerParticleParameter("hkz")
        sys.addForce(force)

        q = """SELECT p0, x0, y0, z0, fcx, fcy, fcz FROM posre_harm_term INNER JOIN posre_harm_param ON posre_harm_term.param=posre_harm_param.id"""

        for (fcounter,conn,tables,offset) in self._localVars():
            if not go[fcounter]:
                continue
            for p0, x0, y0, z0, fcx, fcy, fcz in conn.execute(q):
                p0 += offset
                x0d = (x0*angstrom).value_in_unit(nanometer)
                y0d = (y0*angstrom).value_in_unit(nanometer)
                z0d = (z0*angstrom).value_in_unit(nanometer)
                hfcxd = (0.5*fcx*kilocalorie_per_mole/angstrom**2).value_in_unit(kilojoule_per_mole/(nanometer**2))
                hfcyd = (0.5*fcy*kilocalorie_per_mole/angstrom**2).value_in_unit(kilojoule_per_mole/(nanometer**2))
                hfczd = (0.5*fcz*kilocalorie_per_mole/angstrom**2).value_in_unit(kilojoule_per_mole/(nanometer**2))
                force.addParticle(p0,[ x0d, y0d, z0d, hfcxd,  hfcyd,  hfczd])

                
    def _hasTable(self, table_name, tables):
        """check existence of a table
        """
        return table_name in tables

    
    def _readSchemas(self, conn):
        """Read and return the schemas of each of the tables in the dms file connection 'conn'"""
        tables = {}
        for table in conn.execute("SELECT name FROM sqlite_master WHERE type='table'"):
            names = []
            for e in conn.execute('PRAGMA table_info(%s)' % table):
                names.append(str(e[1]))
            tables[str(table[0])] = names
        return tables


    def _checkForUnsupportedTerms(self):
        """Check the file for forcefield terms that are not currenty supported,
        raising a NotImplementedError
        """
        flat_bottom_potential_terms = ['stretch_fbhw_term', 'angle_fbhw_term',
                                       'improper_fbhw_term', 'posre_fbhw_term']

        for (fcounter,conn,tables,offset) in self._localVars():
            if any((t in tables) for t in flat_bottom_potential_terms):
                raise NotImplementedError('Flat bottom potential terms '
                                          'are not implemeneted')
            nbinfo = dict(zip(tables['nonbonded_info'],
                              conn.execute('SELECT * FROM nonbonded_info').fetchone()))
            self._opls_combining_rules = False
            if nbinfo['vdw_funct'] != u'vdw_12_6':
                raise NotImplementedError('Only Lennard-Jones van der Waals '
                                          'interactions are currently supported')
            if (nbinfo['vdw_rule'] != u'arithmetic') and (nbinfo['vdw_rule'] != u'geometric') and (nbinfo['vdw_rule'] != u'arithmetic/geometric'):
                raise NotImplementedError('Unknown vdw combination rule')
            if nbinfo['vdw_rule'] == u'geometric':
                #need to generalize this
                self._opls_combining_rules = True
            if nbinfo['vdw_rule'] == u'arithmetic':
                self._opls_combining_rules = False
            """
            if 'nonbonded_combined_param' in tables:
                raise NotImplementedError('nonbonded_combined_param interactions '
                                          'are not currently supported')
            """
            if 'alchemical_particle' in tables:
                raise NotImplementedError('Alchemical particles are not supported')
            if 'alchemical_stretch_harm' in tables:
                raise NotImplementedError('Alchemical bonds are not supported')
            if 'polar_term' in tables:
                if conn.execute("SELECT COUNT(*) FROM polar_term").fetchone()[0] != 0:
                    raise NotImplementedError('Drude particles are not currently supported')

    def _createProvenance(self):
        for (fcounter, conn, tables, offset) in self._localVars():
            try:
                # build the provenance string
                provenance = []
                q = """SELECT id, user, timestamp, version, workdir, cmdline, executable
                FROM provenance"""
                for id, user, timestamp, version, workdir, cmdline, executable in conn.execute(q):
                    for row in conn.execute('SELECT * FROM provenance'):
                        rowdict = dict(list(zip(tables['provenance'], row)))
                        provenance.append('%(id)d) %(timestamp)s: %(user)s\n  version: %(version)s\n  '
                                          'cmdline: %(cmdline)s\n  executable: %(executable)s\n' % rowdict)
                    self._provenance.append(''.join(provenance))
            except:
                self._provenance.append('')
                if self._verbose:
                    print('Warning: unable to retrieve provenance information from DMS file %s' % str(f))

    def _prefixsum(self, values):
        """ exclusive prefix sum of 'values' """
        total =0
        sum = [total]
        for v in values:
            total += v
            sum.append(total)
        return sum

    def _localVars(self):
        counter = range(0,len(self._conn))
        return zip(counter, self._conn, self._tables, self._offset)

    def close(self):
        """Close the SQL connections
        """
        for (fcounter,conn,tables,offset) in self._localVars():
            if self._open[fcounter]:
                conn.close()

    def __del__(self):
        self.close()

    def getEnergy(self):
        integrator = mm.LangevinMiddleIntegrator(310*kelvin, 1/picosecond, 0.002*picoseconds)
        simulation = Simulation(self.topology, self.system, integrator)
        simulation.context.setPositions(self.getPositions())
        print('{:7s}: {:10.3f} kJ/mol'.format('openMM', getEnergy(simulation)))
 

    def runEM(self, out=None):
        """Perform Energy Minimization

        Parameters
        ----------
        out : str
            Filename for the output DMS
        """

        integrator = mm.LangevinMiddleIntegrator(310*kelvin, 1/picosecond, 0.002*picoseconds)
        simulation = Simulation(self.topology, self.system, integrator)
        simulation.context.setPositions(self.getPositions())
        platform = simulation.context.getPlatform().getName()
        print('-------------------------------')
        print("Platform: ", platform)
        print('E0: %.5e kJ/mol' %getEnergy(simulation))
        simulation.minimizeEnergy()
        print('E1: %.5e kJ/mol' %getEnergy(simulation))
        print('-------------------------------')

        ### SAVE THE LAST FRAME
        self.positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)._value
        
        if out:
            shutil.copy(self._file[0], out)
            conn    = sqlite3.connect(out)
            n_atoms = conn.execute('SELECT COUNT(DISTINCT id) FROM particle;').fetchone()[0]
            assert n_atoms == len(self.positions), 'n_atoms is different'

            for index in range(n_atoms):
                conn.execute('UPDATE particle SET x = ?, y = ?, z = ? WHERE id = ?', 
                    (*self.positions[index][0:3] * 10, index))

            conn.commit()
            conn.close()


    def runEMNPT(self, out=None, emout=None,
        nonbondedCutoff=1.1, 
        nsteps=10000, dcdfreq=1000, csvfreq=1000, dt=0.02, P=1.0, T=310, tension=0.0,
        semiisotropic=False, barfreq=100, addForces=[], frictionCoeff=5.0, EM=True):
        '''Perform EM, followed by NPT. I used this mostly to run Martini simulations.
        Therefore, the default parameters (dt=0.02) are tuned for martini simulations.
        However, this does not mean that you cannot use this function in all-atom simulations.
        Just change the parameters appropriately for the all-atom simulations.
        For all-atom, dt = 0.002 and frictionCoeff = 1
    
        If you want to use EM + NVT (not NPT) simulations, set barfreq=0 
    
        Parameters
        ----------
        out : str
            Filename for the dms (e.g., npt.dms)
        nonbondedCutoff : float=1.1
            Cutoff distance for vdw in nm.
        dt : float=0.02
            Integration time step. Note that you have to change this to 0.002 if you use all-atom resolution.
        nsteps : int=10000
            Number of NPT steps
        dcdfreq : int=1000
            dcd frequency
        csvfreq : int=1000
            Simulaitonal data frequencuy
        '''
        
        ### ADD BAROSTAT
        if barfreq > 0:
            if semiisotropic:
                self.system.addForce(MembranePressure(P=P, T=T, r=tension, barfreq=barfreq))
            else:
                self.system.addForce(Pressure(P=P, T=T, barfreq=barfreq))
    
        ### ADD additional forces
        for addf in addForces:
            self.system.addForce(copy.deepcopy(addf))
        
        ### PREPARE SIMS
        integrator = mm.LangevinMiddleIntegrator(T*kelvin, frictionCoeff/picosecond, dt*picoseconds)
        simulation = Simulation(self.topology, self.system, integrator)
        simulation.context.setPositions(self.getPositions())
        simulation.context.setPeriodicBoxVectors(self.cell[0], self.cell[1], self.cell[2])
        platform = simulation.context.getPlatform().getName()
        
        # boxSize = self.topology.getUnitCellDimensions()
        # unitCellDimensions = [boxVectors[0][0], boxVectors[1][1], boxVectors[2][2]]
        # top.setUnitCellDimensions(unitCellDimensions*angstrom)

        ### RUN EM
        print('-------------------------------')
        print("Platform: ", platform)
        print('dt: %.1f fs' %(dt * 1e3))
        print('E0: %.5e kJ/mol' %getEnergy(simulation))
        if EM: simulation.minimizeEnergy()
        print('E1: %.5e kJ/mol' %getEnergy(simulation))
        print('-------------------------------')

        ### SAVE EM
        if emout:
            self.positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)._value
            shutil.copy(self._file[0], emout)
            conn    = sqlite3.connect(emout)
            cursor  = conn.cursor()
            n_atoms = cursor.execute('SELECT COUNT(DISTINCT id) FROM particle;').fetchone()[0]
    
            for index in range(n_atoms):
                cursor.execute('UPDATE particle SET x = ?, y = ?, z = ? WHERE id = ?', 
                    (*self.positions[index][0:3] * 10, index))
 
            conn.commit()
            conn.close()
    
        ### RUN NVT/NPT
        if out:
            prefix = '.'.join(out.split('.')[:-1])
            simulation.reporters.append(DCDReporter(      prefix + '.dcd', dcdfreq))
            simulation.reporters.append(StateDataReporter(prefix + '.csv', csvfreq, step=True, potentialEnergy=True, temperature=True))
        
        simulation.context.setVelocitiesToTemperature(T)
        simulation.step(nsteps)
        self.positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)._value
        self.cell = getCell(simulation)
    
        ### SAVE THE LAST FRAME
        if out is not None and out.split('.')[-1] == 'dms':
            shutil.copy(self._file[0], out)
    
            conn    = sqlite3.connect(out)
            cursor  = conn.cursor()
            n_atoms = cursor.execute('SELECT COUNT(DISTINCT id) FROM particle;').fetchone()[0]
    
            assert n_atoms == len(self.positions), 'n_atoms is different'
    
            for index in range(n_atoms):
                cursor.execute('UPDATE particle SET x = ?, y = ?, z = ? WHERE id = ?', 
                    (*self.positions[index][0:3] * 10, index))
            
            # global_cell id starts from 1 while particle id starts from 0 (why?)
            ids = sorted([tmp[0] for tmp in cursor.execute('SELECT id FROM global_cell;').fetchall()])
            assert len(set(ids)) == 3, f'global_cell has {len(set(ids))} entries (it should have been 3)'
            cursor.execute('UPDATE global_cell SET x=?, y=?, z=? WHERE id=?', (*self.cell[0] * 10, ids[0]))
            cursor.execute('UPDATE global_cell SET x=?, y=?, z=? WHERE id=?', (*self.cell[1] * 10, ids[1]))
            cursor.execute('UPDATE global_cell SET x=?, y=?, z=? WHERE id=?', (*self.cell[2] * 10, ids[2]))
    
            conn.commit()
            conn.close()


def getEnergy(simulation):
    return simulation.context.getState(getEnergy=True).getPotentialEnergy()._value

def getCell(simulation):
    return simulation.context.getState().getPeriodicBoxVectors(asNumpy=True)._value

def getPositions(simulation):
    return simulation.context.getState(getPositions=True).getPositions(asNumpy=True)._value

def Pressure(P=1.0, T=310.0, barfreq=100):
    return MonteCarloBarostat(P*bar, T*kelvin)

def MembranePressure(P=1.0, T=310.0, r=0.0, barfreq=100):
    return MonteCarloMembraneBarostat(P*bar, r*bar*nanometer, T*kelvin,
             MonteCarloMembraneBarostat.XYIsotropic, MonteCarloMembraneBarostat.ZFree)

Backmapping
===========

**mstool** is a multiscale simulation tool that backmaps a coarse-grained structure into an all-atom structure (`Backmap`) and reviews the resulting structure (`CheckStructure`). It requires minimal user input (mapping and isomeric information of molecules) and is more powerful than the previous backmapping tools. For more details about the methodology, please check out the `paper <https://pubs.acs.org/doi/abs/10.1021/acs.jpcb.3c05593>`_. The beta version is available from `github.com/ksy141/mstool <https://github.com/ksy141/mstool>`_. Please report bugs or ask questions on the `Github Issues <https://github.com/ksy141/mstool/issues>`_.

Input
-----

Structure
^^^^^^^^^
mstool can read and write **PDB** or **DMS** files. **GRO** files are not supported because they do not contain chain information. Note that each residue in a structure file should have a unique set of (resname, resid, chain) to use mstool. **In other words, no more than one residue should have the same resname, resid, and chain name.** Large systems can bump into this issue because PDB only supports a residue number from 0 to 9999 because of a limit of PDB formatting on the length of resid. For instance, water molecules can easily hit this limit if they all have the same chain name. To avoid this problem, reassign the chain names of water for every 10,000 water molecules or use DMS files that do not limit the length of chain, resid, or resname.

Mapping
^^^^^^^
Files containing **mapping** and **isomeric** information, referred to as simply *mapping files* throughout the documentation, are central and required input files in mstool. They describe which atoms belong to which coarse-grained beads and isomeric information (**cis/trans/chiral/dihedral**) if molecules are isomers. If no user mapping files are provided, mstool will read the predefined mapping files for Martini force fields, ``$mstoolpath/mapping/martini.protein.c36m.dat`` and ``$mstoolpath/mapping/martini.lipid.c36.dat``.

For each residue type, write a residue name, following a keyword, RESI. Put the name of a coarse-grained bead in a square bracket and an atom name that belongs to the coarse-grained bead in the following line. For instance, consider a methane molecule, CH4. Suppose a coarse-grained force field has one bead per each methane molecule whose coarse-grained bead name is CG1. A methane mapping file (residue name METH) will look something like:

.. code-block:: text

   RESI METH
   [ CG1 ]
   C1 H1 H2 H3 H4

If a molecule has a chiral center, and you want backmapped molecules that are enantiomerically pure rather than racemic, specify five atoms using their atomic names to describe each chiral center.

.. code-block:: text

   [ chiral ]
   A B C D E

In the above example, Atom B is the central atom of a tetrahedron. Curl your right fingers in the order of C -> D -> E. Your right thumb should point in the direction of Atom A, perpendicular to the C-D-E triangle. 

Users can also set a dihedral angle. A dihedral of A, B, C, and D atoms will be 90 degrees in a backmapped structure with the below input.

.. code-block:: text

   [ dihedral ]
   A B C D 90.0

Geometric isomerism can be defined using *dihedral*. For instance, the below two represent a *cis* arrangement of four atoms and are equivalent. 

.. code-block:: text

   [ cis ]
   A B C D

.. code-block:: text

   [ dihedral ]
   A B C D 0.0

Similarly, the below two define a *trans* arrangement of four atoms and are equivalent.

.. code-block:: text

   [ trans ]
   A B C D

.. code-block:: text

   [ dihedral ]
   A B C D 180.0


Atomic arrangements of the supported isomeric properties in a mapping file are shown in the below :ref:`figure <isomer>`.

.. _isomer:

.. figure:: _static/isomer.png
   :scale: 70%
    
   Atomic arrangements of chiral, dihedral, cis, and trans.


Force field
^^^^^^^^^^^
All-atom, openMM-formatted (**XML**) force fields should be provided to use mstool. If no user force fields are provided, mstool will read the default charmm36 force fields, ``$mstoolpath/FF/charmm36/charmm36.xml`` and ``$mstoolpath/FF/charmm36/water.xml``.

.. note::

   Residue names, coarse-grained bead names, and atomic names should be consistent between these input files, structures, mapping files, and force fields.


Workflows
---------

Backmapping consists of three steps in mstool: Ungrouping, Energy Minimization, and Reviewing a backmapped structure. 

Ungrouping is placing atoms near their corresponding coarse-grained beads, implemented in ``Ungroup``, based on an input mapping file. The atom placement is random except for protein, which we will discuss later.

The all-atom structure after ``Ungroup`` has very high potential energy and cannot be relaxed with standard energy minimization. The structure should undergo special energy minimization, called Reduced Nonbonded Energy Minimization (REM), implemented in ``REM``. During REM, nonbonded interactions, both Lennard-Jones and Coulomb interactions, are replaced with a soft, cosine function for numerical stability. If molecules are isomers and their isomeric properties are specified in the mapping files (e.g., *cis/trans/chiral/dihedral*), dihedral potentials will be applied to the relevant atoms to keep the desired isomeric properties. For instance, the below will apply a set of dihedral potentials to keep a protein residue, ALA, to be L-alanine (or S-alanine) during REM. 

.. code-block:: text

   RESI ALA
   [ BB ]
   N HN CA HA C O
   CB HB1 HB2 HB3

   [ chiral ]
   HA CA N CB C

Finally, the relaxed structure will be reviewed by ``CheckStructure``. It confirms whether the isomeric properties of backmapped molecules are consistent with those written in the input mapping files. Besides the isomeric information in the input mapping files, the tool will automatically detect protein peptide bonds based on resname, resid, chain, and name and report any *cis* peptide bonds. 

The workflow that does the three tasks is implemented in ``Backmap``.

Examples
--------

Example 1. Methane
^^^^^^^^^^^^^^^^^^
The first example is to backmap a coarse-grained system containing 125 molecules, each of which is represnted by a single coarse-grained bead. In this example, let's consider this molecule as methane, CH4. Generic residue and coarse-grained bead names were assigned: ONE and CG1, respectively. To copy the coarse-grained structure file into your current directory, execute the following Python script (or copy it by yourself. The file is ``$mstoolpath/examples/Backmapping/Example1_methane/cg.pdb``):

.. code-block:: python

   import mstool
   import shutil
   shutil.copyfile(mstool.ONECGBEADSTRUCTURE, './cg.pdb')

A mapping file is a central input in mstool. Each coarse-grained bead whose residue name is ONE and whose name is CG1 should represent five atoms whose atomic names are set to C1, H1, H2, H3, and H4. Make the following mapping file and save it as ``mapping.dat``. Note that the residue name is ONE, and the bead name is CG1 to be consistent with the structure.

.. code-block:: text

   RESI ONE
   [ CG1 ]
   C1 H1 H2 H3 H4

Finally, the methane topology is not defined in the standard CHARMM36 force field; therefore, should be added for backmapping. Make the following force field file and save it as ``ff.xml``. Note that the residue name is ONE, and the bead name is CG1 to be consistent with the structure and mapping files. Also, the atomic names in the force field file should be consistent with those defined in the mapping file, which are C1, H1, H2, H3, and H4.

.. code-block:: xml

  <ForceField>
    <Residues>
      <Residue name="ONE">
        <Atom charge="-0.36" name="C1" type="CT3"/>
        <Atom charge="+0.09" name="H1" type="HA3"/>
        <Atom charge="+0.09" name="H2" type="HA3"/>
        <Atom charge="+0.09" name="H3" type="HA3"/>
        <Atom charge="+0.09" name="H4" type="HA3"/>
        <Bond atomName1="C1" atomName2="H1"/>
        <Bond atomName1="C1" atomName2="H2"/>
        <Bond atomName1="C1" atomName2="H3"/>
        <Bond atomName1="C1" atomName2="H4"/>
      </Residue>
    </Residues>
  </ForceField>
           
All the ingredients are ready! Consistency in atomic names, bead names, and residue names between structure, mapping, and force field files are important in mstool. Please double-check whether they have consistent names like the :ref:`the below figure <consistency>`:

.. _consistency:
.. figure:: _static/consistency.png
    
    Consistency in names between inputs. 

Let's backmap the coarse-grained structure.

.. code-block:: python

   import mstool
   mstool.Ungroup('cg.pdb', 'aa.pdb', mapping='mapping.dat')
   mstool.REM('aa.pdb', 'aa_final.pdb', mapping='mapping.dat', ff_add='ff.xml')
   mstool.CheckStructure('aa_final.pdb', mapping='mapping.dat')

``Ungroup`` makes an intermediate all-atom structure. ``REM`` relaxes the structure using Reduced Nonbonded Energy Minimization. ``CheckStructure`` reviews the isomeric properties of the relaxed structure, although in this case, there is nothing to review because methane is not an isomer, and no isomeric information is defined in ``mapping.dat``. The initial coarse-grained, intermediate all-atom, and final all-atom structures are visualized in :ref:`Methane backmapping <methane>`.

.. _methane:
.. figure:: _static/methane.png

   Methane backmapping. (left) Initial coarse-grained structure. (center) Intermediate all-atom structure. (right) Final all-atom structure.

The one-line backmapping procedure is available. The below executes all of the three steps inside the workflow:

.. code-block:: python
   
   import mstool
   mstool.Backmap('cg.pdb', mapping='mapping.dat', ff_add='ff.xml')

The final all-atom structure is ``workdir/step4_final.pdb``.

.. note::

  The methane force field was created for backmapping. Since this force field was not validated against experimental data, it should not be used for production molecular dynamics simulations.

.. note:: 

    ``ff='ff.xml'`` and ``ff_add='ff.xml'`` arguments in ``REM`` are not equivalent. The former only reads ``ff.xml``. The latter reads the default CHARMM36 force field files, ``$mstoolpath/FF/charmm36/charmm36.xml`` and ``$mstoolpath/FF/charmm36/water.xml``, and then reads the additionally provided ``ff.xml``. In this example, the methane force field should be provided as an additional force field because it only defines the topology of methane but not the necessary parameters for simulations (e.g., epsilon and sigma of CT3 and HA3 and bond parameters of CT3-HA3). In most cases, you do not need ``ff=ff.xml`` but ``ff_add=ff.xml`` when you have a new molecule not defined in the standard CHARMM36 force field.

Example 2. Ethane
^^^^^^^^^^^^^^^^^

In the second example, we will backmap the same coarse-grained system of the first example to ethane, C2H6. As explained above, this toy coarse-grained system has 125 molecules, each of which is represented by a single coarse-grained bead. The residue name and bead name are ONE and CG1, respectively. The ethane force field is available in the standard CHARMM36 force field. The residue name in the standard force field is ETHA. Check this in ``$mstoolpath/FF/charmm36/charmm36.xml``. To make the residue name in our structure consistent with the force field, let's change the residue name from ONE to ETHA:

.. code-block:: python

   import mstool
   u = mstool.Universe(mstool.ONECGBEADSTRUCTURE)
   u.atoms.resname = 'ETHA'
   u.write('cg.pdb')

The next task is to make a mapping file for ethane. A coarse-grained bead named CG1 should be backmapped to eight atoms whose names are C1 H11 H12 H13 C2 H21 H22 H23 in the standard CHARMM36 force field. Make the following mapping file and save it as ``mapping.dat``.

.. code-block:: text

   RESI ETHA
   [ CG1 ]
   C1 H11 H12 H13
   C2 H21 H22 H23

An additional force field is not required in this case because ethane is already defined in the CHARMM36 force field. Let's backmap the structure:

.. code-block:: python

   import mstool
   mstool.Ungroup('cg.pdb', 'aa.pdb', mapping='mapping.dat')
   mstool.REM('aa.pdb', 'aa_final.pdb', mapping='mapping.dat')
   mstool.CheckStructure('aa_final.pdb', mapping='mapping.dat')

The one-line backmapping procedure also works. The final structure is ``workdir/step4_final.pdb``.

.. code-block:: python

    import mstool
    mstool.Backmap('cg.pdb', mapping='mapping.dat')

The initial coarse-grained, intermediate all-atom, and final all-atom structures are visualized in :ref:`Ethane backmapping <ethane>`.

.. _ethane:

.. figure:: _static/ethane.png

   Ethane backmapping. (left) Initial coarse-grained structure. (center) Intermediate all-atom structure. (right) Final all-atom structure.

Example 3. trans-2-butene
^^^^^^^^^^^^^^^^^^^^^^^^^

In this example, we will backmap the same coarse-grained system of the first example to trans-2-butene, C4H8. As explained above, this toy coarse-grained system has 125 molecules, each of which is represented by a single coarse-grained bead. The residue name and bead name are ONE and CG1, respectively. The 2-butene force field is available in the standard CHARMM36 force field. The residue name in the standard force field is BTE2. Check this in ``$mstoolpath/FF/charmm36/charmm36.xml``. To make the residue name in our structure consistent with the force field, let's change the residue name from ONE to BTE2:

.. code-block:: python

   import mstool
   u = mstool.Universe(mstool.ONECGBEADSTRUCTURE)
   u.atoms.resname = 'BTE2'
   u.write('cg.pdb')

Let's make a mapping file for trans-2-butene and save it as ``mapping.xml``. Note that trans-2-butene has a double bond between the two central carbon atoms; therefore, it is a geometric isomer. A desired isomeric property should be written in the mapping file, in this case, *trans*. To be consistent with the coarse-grained structure file, the residue name, BTE2, and the coarse-grained bead name, CG1, should be used in the mapping file. The atomic names should be consistent with the force field. 

.. code-block:: text

   RESI BTE2
   [ CG1 ]
   C1 H11 H12 H13
   C2 H21
   C3 H31
   C4 H41 H42 H43

   [ trans ]
   C1 C2 C3 C4

Let's backmap the structure:

.. code-block:: python
   
   import mstool
   mstool.Backmap('cg.pdb', mapping='mapping.dat')

The final structure is ``workdir/step4_final.pdb``. At the end of the workflow, the tool reviews (``CheckStructure``) whether there are any cis-2-butene. Your backmapped structure should be good if you see the following message:

.. code-block:: text

  ####################################################################
  workdir/step3_em.dms was reviewed
  ####################################################################
  
  The following isomers were reviewed:
  trans: resname BTE2 - C1 C2 C3 C4
  ####################################################################
  
  ####################################################################
  No molecules had flipped isomers
  ####################################################################
  
  ####################################################################
  In summary, the number of residues with the flipped isomers:
  peptide   :          0
  cistrans  :          0
  chiral    :          0
  dihedral  :          0
  ####################################################################
  
  Adding bonds for non-protein residues - started
  Adding bonds for non-protein residues - finished
  ####################################################################
  Tetrahedron checking - started
  Tetrahedron checking - finished
  ####################################################################

.. note::

   Geometric isomerism is not specified in force fields. In other words, cis-2-butene and trans-2-butene have exactly the same force field. A dihedral potential is internally applied during REM to ensure a backmapped molecule has the desired isomeric property as written in the mapping file.


Example 4. Martini POPC
^^^^^^^^^^^^^^^^^^^^^^^

In this example, we will backmap a Martini POPC bilayer. Copy the coarse-grained structure into your current directory by executing the following Python script (or copy it by yourself. The file is ``$mstoolpath/examples/Backmapping/Example4_POPC/cg.pdb``):

.. code-block:: python

   import mstool
   import shutil
   shutil.copyfile(mstool.POPCSTRUCTURE, './cg.pdb')

There is no need to define a mapping file for Martini POPC because it is already available in ``$mstoolpath/mapping/martini.lipid.c36.dat``. Note that :ref:`POPC <popc>` has one chiral center and one cis bond, which is already defined in the default mapping file. Also, an additional force field is not required because POPC is defined in ``$mstoolpath/FF/charmm36/charmm36.xml``. Backmapping is as simple as the following:

.. code-block:: python

   import mstool
   mstool.Backmap('cg.pdb')

The backmapped structure is ``workdir/step4_final.pdb``. The isomeric properties of backmapped molecules are reviewed at the end of the workflow.

.. _popc:
.. figure:: _static/popc.png

   POPC backmapping. (A) Molecular structure of POPC. (B) Initial coarse-grained and backmapped structures.


One thing to note is that mstool assumes the coarse-grained water resname is W. If the water resname is not W in your coarse-grained structure, it should be provided to ``Backmap`` or ``Ungroup``. Let's assume your water resname is WAT in your coarse-grained structure.
Also, mstool uses the 1-to-4 mapping for water by default to be consistent with the Martini force field. That is, each coarse-grained water bead represents four all-atom water molecules. If you want to change this to 1-to-n mapping (e.g., using other coarse-grained force fields or wanting more or less water in a backmapped structure):

.. code-block:: python

   mstool.Backmap(..., water_resname='WAT', water_number=n)
   mstool.Ungroup(..., water_resname='WAT', water_number=n)


Example 5. Martini bilayer
^^^^^^^^^^^^^^^^^^^^^^^^^^

In the previous example, we backmapped a Martini POPC bilayer. In this example, let's backmap a multi-component, spherical bilayer at Martini resolution, which is shown in `Figure 5B of the mstool publication <https://pubs.acs.org/doi/abs/10.1021/acs.jpcb.3c05593>`_. Copy the coarse-grained strucutre (``$mstoolpath/examples/Backmapping/Example5_Sphere/cg.pdb``) into the current directory:

.. code-block:: python

   import mstool
   import shutil
   shutil.copyfile(mstool.MULTISTRUCTURE, './cg.pdb')

Like the previous POPC bilayer, all the lipids included in the system have mapping information in the default mapping file, ``$mstoolpath/mapping/martini.lipid.c36.dat``. :ref:`Lipids <lipids>` with predefined mapping files are shown below:

.. _lipids:
.. figure:: _static/lipids.png

   Lipids with predefined mapping files.

Because all of our lipids are already defined in the default mapping file, making a new mapping file is unnecessary. Also, the CHARMM36 force field already has parameters and topologies for these lipids. Let's backmap by simply executing the following Python script (because this is a large system, it will take ~20 mins or more):

.. code-block:: python

   import mstool
   mstool.Backmap('cg.pdb')

This example has many isomeric properties to be reviewed because each lipid has at least one chiral center. If you see the report that no molecules had flipped isomers, your structure should be good to start a production run.

.. code-block:: text

   ####################################################################
   workdir/step3_em.dms was reviewed
   ####################################################################
   
   The following isomers were reviewed:
   chiral: resname POPG - HS C2 O21 C1 C3
   chiral: resname POPG - O13 P O11 O14 O12
   chiral: resname POPG - H12A C12 OC2 C13 C11
   chiral: resname CHL1 - H3 C3 O3 C2 C4
   chiral: resname CHL1 - C19 C10 C1 C5 C9
   chiral: resname CHL1 - H9 C9 C8 C10 C11
   chiral: resname CHL1 - H8 C8 C9 C7 C14
   chiral: resname CHL1 - H14 C14 C8 C13 C15
   chiral: resname CHL1 - C18 C13 C12 C14 C17
   chiral: resname CHL1 - H17 C17 C13 C20 C16
   chiral: resname CHL1 - H20 C20 C21 C17 C22
   chiral: resname DOPA - HS C2 O21 C1 C3
   chiral: resname DOPA - O13 P O11 O14 O12
   chiral: resname POPC - HS C2 O21 C1 C3
   chiral: resname POPC - O13 P O11 O14 O12
   chiral: resname DOPG - HS C2 O21 C1 C3
   chiral: resname DOPG - O13 P O11 O14 O12
   chiral: resname DOPG - H12A C12 OC2 C13 C11
   chiral: resname POPS - HS C2 O21 C1 C3
   chiral: resname POPS - O13 P O11 O14 O12
   chiral: resname POPS - H12A C12 N C11 C13
   chiral: resname DPPC - HS C2 O21 C1 C3
   chiral: resname DPPC - O13 P O11 O14 O12
   chiral: resname DOPE - HS C2 O21 C1 C3
   chiral: resname DOPE - O13 P O11 O14 O12
   chiral: resname POPA - HS C2 O21 C1 C3
   chiral: resname POPA - O13 P O11 O14 O12
   chiral: resname DOPC - HS C2 O21 C1 C3
   chiral: resname DOPC - O13 P O11 O14 O12
   chiral: resname DOPS - HS C2 O21 C1 C3
   chiral: resname DOPS - O13 P O11 O14 O12
   chiral: resname DOPS - H12A C12 N C11 C13
   chiral: resname POPE - HS C2 O21 C1 C3
   chiral: resname POPE - O13 P O11 O14 O12
   cis: resname POPG - C28 C29 C210 C211
   cis: resname DOPA - C28 C29 C210 C211
   cis: resname DOPA - C38 C39 C310 C311
   cis: resname POPC - C28 C29 C210 C211
   cis: resname DOPG - C28 C29 C210 C211
   cis: resname DOPG - C38 C39 C310 C311
   cis: resname POPS - C28 C29 C210 C211
   cis: resname POPS - H91 C29 C210 H101
   cis: resname DOPE - C28 C29 C210 C211
   cis: resname DOPE - C38 C39 C310 C311
   cis: resname POPA - C28 C29 C210 C211
   cis: resname DOPC - C28 C29 C210 C211
   cis: resname DOPC - C38 C39 C310 C311
   cis: resname DOPS - C28 C29 C210 C211
   cis: resname DOPS - C38 C39 C310 C311
   cis: resname POPE - C28 C29 C210 C211
   ####################################################################
   
   ####################################################################
   No molecules had flipped isomers
   ####################################################################
   
   ####################################################################
   In summary, the number of residues with the flipped isomers:
   peptide   :          0
   cistrans  :          0
   chiral    :          0
   dihedral  :          0
   ####################################################################

   Adding bonds for non-protein residues - started
   Adding bonds for non-protein residues - finished
   ####################################################################
   Tetrahedron checking - started
   <Atom 76603 (C25) of chain 1 residue 614 (DOPE)> 1.5542432896731566
   <Atom 170544 (C214) of chain 1 residue 1372 (DOPS)> 1.6635004939312479
   Tetrahedron checking - finished
   ####################################################################

At the end of the check, the tool also reviews whether atoms have a good tetrahedron geometry. I got two warnings from my backmapped structure. However, this is harmless as it is not an isomeric property and will be quickly fixed within ~0.1 ns of MD simulations. If you want your backmapped structure more equilibrated, which will lower the chance of getting these tetrahedron warnings, increase the number of NVT steps. The default is 10000, which is 2 ps.

.. code-block:: python
   
   mstool.Backmap(..., nsteps=)
   mstool.REM(..., nsteps=)

Initial coarse-grained, intermediate all-atom, and final all-atom structures are shown in the following :ref:`figure <bilayer>`.

.. _bilayer:

.. figure:: _static/bilayer.png

   Multi-component, Martini, bilayer backmapping. (left) Initial coarse-grained structure. (center) Intermediate all-atom structure. (right) Final all-atom structure.


Example 6. Triolein
^^^^^^^^^^^^^^^^^^^

The previous spherical bilayer contains the Martini lipids with the default mapping files. What should you do if you have a new Martini lipid not supported by default in mstool? In this example, a bilayer membrane contains POPC and a neutral lipid, triolein (resname TRIO). The mapping file and force field for POPC already exist; However, TRIO is a new molecule not defined in the default mapping files and the standard CHARMM36 force field. Therefore, new files should be made for TRIO. Let's copy the coarse-grained structure, mapping file of TRIO, and force field of TRIO into the current directory (or copy them by yourself. The path is ``$mstoolpath/examples/Backmapping/Example6_TRIO``

.. code-block:: python
   
   import mstool
   import shutil
   shutil.copyfile(mstool.TRIOSTRUCTURE, './cg.pdb')
   shutil.copyfile(mstool.TRIOMAPPING, './mapping.dat')
   shutil.copyfile(mstool.TRIOFF, './ff.xml')

Review whether the mapping file and force field of TRIO look reasonable to you. The partial charges of TRIO were obtained from `Biophys. Rep., 2021, 1, 2, 100034. <https://www.cell.com/biophysreports/fulltext/S2667-0747(21)00034-3>`_

Let's backmap the structure. Provide the TRIO mapping file as an additional mapping file (``mapping_add='mapping.dat'``). If you provide this as a standalone mapping file (``mapping='mapping.dat'``), the default files, which contain the mapping information of POPC, will not be read. Similarly, provide the TRIO force field file as an additional force field file (``ff_add='ff.xml'``) rather than as a standalone force field file (``ff='ff.xml'``). The final structure is ``workdir/step4_final.pdb``.

.. code-block:: python

   import mstool
   mstool.Backmap('cg.pdb', mapping_add='mapping.dat', ff_add='ff.xml')

.. note::

   openMM does not allow two residues that have the same resid and chain. Review whether your structure has two or more residues with the same resid and chain. The coarse-grained structure in this example has two ions: SOD and CLA. For each NaCl pair, SOD and CLA have the same resid in this structure. However, their chains differ, so openMM does not complain about this.


Example 7. Membrane Protein
^^^^^^^^^^^^^^^^^^^^^^^^^^^

We will backmap a membrane protein, ompF porin, in this example. Lipids, in general, are flexible molecules and can be easily and quickly equilibrated (< 50 ns) even if their starting structures are not at equilibrium. Therefore, lipids will be backmapped in the same way we have done.

Protein is different because it is a very long molecule, unlike lipids. Its equilibrium timescale is beyond the all-atom timescale. In other words, if you mis-backmap your protein, your protein will likely have the wrong structure throughout your all-atom trajectory. Think of protein as solid and lipid as liquid. 

What should we do? No backmapping is better than using a true structure, which is an experimentally resolved structure. Assuming that your protein structure does not change too much in your coarse-grained trajectory, you can simply copy your **all-atom** protein structure and then align it against the **coarse-grained** protein structure. Let's copy the coarse-grained and all-atom protein structure into the current directory.

.. code-block:: python

   import mstool
   import shutil
   shutil.copyfile(mstool.MPCG, './cg.pdb')
   shutil.copyfile(mstool.MPAA, './protein_AA.pdb')

The all-atom protein structure is prealigned against the coarse-grained protein structure. Review whether the conformation and location of these two structures are reasonably the same. We have to separate the coarse-grained protein structure from the coarse-grained nonprotein structure.

.. code-block:: python

   import mstool
   u = mstool.Universe('cg.pdb')

   # select nonprotein
   non_protein_bA = ~u.atoms.resname.isin(mstool.three2one.keys())

   # make a nonprotein universe and save it
   non_protein = mstool.Universe(data=u.atoms[non_protein_bA])
   non_protein.dimensions = u.dimensions
   non_protein.write('cg_nonprotein.pdb')

Let's backmap a structure. Nonprotein molecules are ungrouped from their coarse-grained beads. The all-atom protein structure given as an argument will be used instead of ungrouping coarse-grained beads. The all-atom system then will undergo REM for relaxation. The final structure is ``workdir/step4_final.pdb``.

.. code-block:: python
   
   import mstool
   mstool.Backmap(AA='protein_AA.pdb', structure='cg_nonprotein.pdb')




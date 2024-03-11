MembraneBuilder
===============

The tool incorporates a membrane-building function at an all-atom resolution, using a provided lipid composition similar to `CHARMM-GUI <https://charmm-gui.org>`_. Initially, it constructs a membrane at Martini coarse-grained resolution, conducts coarse-grained simulations for equilibrium, and subsequently converts the coarse-grained structure into an all-atom structure. If a protein structure is provided, the tool will then solvate the protein with lipids. However, instead of converting the coarse-grained protein into an all-atom protein during the backmapping process, the coarse-grained protein will be substituted with the initial input all-atom protein. This ensures that the tool respects the input protein structure.

.. note::

   The tool exclusively produces all-atom structures and does not generate input files, such as topology or parameters, for molecular dynamics engines.

.. note::

   The tool will ionize the system assuming a net charge of zero. As a result, the final system will consistently contain an equal number of positive and negative ions. Prior to running simulations, ensure that you adjust the ion count to neutralize the system accordingly.

.. note::

   The center of the resulting membrane is the origin. If you are using gromacs or openMM, it is recommended for visual clarity to shift the membrane by half the periodic boundary box before running simulations, although this will not affect the simulations themselves.

Example 1. POPC bilayer
-----------------------

The first example is to make a POPC bilayer. The below example will create 40 POPC molecules in each leaflet. The final structure is ``workdir/step7_final.pdb``.

.. code-block:: python

   import mstool
   mstool.BilayerBuilder(upper={'POPC': 40}, lower={'POPC': 40})

Example 2. Heterogeneous bilayer
--------------------------------

A heterogeneous bilayer membrane can be constructed by specifying the type and quantity of lipids. The following example will generate a bilayer membrane containing all the default lipid types supported by the tool.

.. code-block:: python

   import mstool
   lipids = {'DPPC': 5, 'DOPC': 5, 'DMPC': 5, 'DSPC': 5, 'POPC': 5,
             'DOPS': 5, 'POPS': 5, 'POPG': 5, 'DOPG': 5, 'CHL1': 5,
             'POPA': 5, 'DOPA': 5, 'POPE': 5', DOPE': 5}
   mstool.BilayerBuilder(upper=lipids, lower=lipids)

Example 3. Triolein
-------------------

If you wish to create a bilayer membrane using a lipid not included by default in mstool, you must provide a Martini force field specific to the lipid, along with a mapping file containing isomeric details. Additionally, if the lipid force field is not available in the default openMM force field file, you will need to supply an all-atom force field for the lipid. It is crucial to ensure consistency in the residue name and atomic names across these input files. 

Below is an example of creating a POPC bilayer with triolein, which is not included by default in the mstool. The all-atom and coarse-grained force fields for triolein, along with its mapping file, have been pre-made. Therefore, you can copy them and review.

.. code-block:: python

   import mstool
   import shutil
   shutil.copyfile(mstool.TRIOMAPPING, 'mapping.dat')
   shutil.copyfile(mstool.TRIOMARTINI, 'martini.itp')
   shutil.copyfile(mstool.TRIOFF,      'ff.xml')
    
   mstool.BilayerBuilder(upper={'POPC':100, 'TRIO':5},
                         lower={'POPC':100, 'TRIO':5},
                         mapping_add='mapping.dat',
                         martini_add='martini.itp',
                         ff_add='ff.xml')


Example 4. Membrane Protein
---------------------------

In this example, we will construct a membrane containing a membrane protein, specifically ompF porin. *Please note that the input membrane protein must be prealigned relative to a membrane with its center at the origin and its normal oriented along the +z axis.* The tool does not backmap the coarse-grained protein structure but instead employs the input protein structure during the backmapping step. Consequently, the protein structure incorporated into the final structure should closely resemble the input structure. Please review `workdir/protein.pdb` and `workdir/step7_final.pdb`.

.. code-block:: python

   import mstool
   import shutil
   shutil.copyfile(mstool.MPAA2, 'protein.pdb')
   mstool.BilayerBuilder(protein='protein.pdb',
                         upper={'POPC':100},
                         lower={'POPC':100})


Example 5. GPCR
---------------

This is an another example featuring a membrane protein, GPCR. The protein structure contains a ligand molecule that is not supported by the default CHARMM force field. Nevertherless, the tool can still make a bilayer membrane containing the input protein structure.


.. code-block:: python
    
   import mstool
   import shutil
   shutil.copyfile(mstool.GPCR, 'gpcr.pdb')
   mstool.BilayerBuilder(protein='gpcr.pdb', 
                         upper={'POPC':150}, 
                         lower={'POPC':150})


Example 6. Sphere Membrane
--------------------------

In this example, we will create a spherical membrane with a radius of 60 Å. The process is akin to constructing a planar bilayer membrane, with the only distinction being the utilization of ``SphereBuilder`` instead of ``BilayerBuilder``. Due to the large size, I suggest omitting solvent during the backmapping step in the tutorial to save time.

.. code-block:: python

   import mstool
   mstool.SphereBuilder(radius=60,
                        upper={'POPC': 1090, 'CHL1': 10},
                        lower={'POPC':  390, 'CHL1': 10},
                        remove_solvent=True)

Example 7. Sphere Protein
-------------------------

In this example, we will incorporate five copies of ompF and five copies of GPCR into a spherical membrane with a radius of 100Å. Since the protein structures are prealigned relative to a membrane centered at the origin, with its normal direction along the +z axis, our initial step involves distributing 10 proteins across the surface of a hypothetical sphere with a radius of 100Å, followed by the addition of lipids.

.. code-block:: python

   import mstool
   mstool.SphereProtein(radius=100,
                        protein={mstool.Universe(mstool.MPAA2): 5,
                                 mstool.Universe(mstool.GPCR):  5},
                        out='protein.dms')
    
   mstool.SphereBuilder(radius=100,
                        protein='protein.dms',
                        upper={'POPC': 1600, 'DOPS': 990, 'CHL1': 10},
                        lower={'POPC': 1000, 'DOPS': 400, 'CHL1': 10})



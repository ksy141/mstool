.. mstool documentation master file, created by
   sphinx-quickstart on Sun Dec 31 18:27:34 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

**mstool** is a multiscale simulation tool for molecular dynamics simulations. It is a Python library that backmaps a coarse-grained structure into an all-atom structure (`Backmap`) and reviews the resulting structure (`CheckStructure`). It also can relax a system via Reduced Nonbonded Energy Minimization (`REM`) with unphysical contacts between atoms, which normal energy minimization fails to do.

All source code is available from `github.com/ksy141/mstool <https://github.com/ksy141/mstool>`_. Please report bugs or ask questions on the `Github Issues <https://github.com/ksy141/mstool/issues>`_.

When using **mstool** in published work, please cite the following paper:

* Siyoung Kim, `Backmapping with Mapping and Isomeric Information <https://pubs.acs.org/doi/abs/10.1021/acs.jpcb.3c05593>`_, *J. Phys. Chem. B.* 2023, 127, 49, 10488.

.. note::
   This project is under active development.

Contents
--------
.. toctree::
    installation
    input
    backmap
    membranebuilder


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

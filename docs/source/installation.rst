Installation
============

Users can download a stable version of **mstool** from :ref:`PyPi <pypi>` or the up-to-date version from :ref:`GitHub <github>`. In either case, **openMM** should be installed with conda first.

.. code-block:: console

   $ conda install -c conda-forge openmm


.. _pypi:

Option 1. PyPi
--------------

Once you install openMM with conda, mstool can be easily installed through PyPi:

.. code-block:: console

   $ pip install mstool

Optionally, if one wishes to register a shortcut for the base path of mstool to a shell config (``$HOME/.bashrc`` or ``$HOME/.zshrc``):

.. code-block:: console

   $ dir=$(pip show mstool | grep Location | awk '{print $2}')
   $ echo "mstoolpath=$dir/mstool" >> yourshellconfig


.. _github:

Option 2. GitHub
----------------

Once you install openMM with conda, download mstool from GitHub and install it:

.. code-block:: console

   $ git clone git@github.com:ksy141/mstool.git
   $ cd mstool
   $ pip install .

Optionally, if one wishes to register a shortcut for the base path of mstool to a shell config (``$HOME/.bashrc`` or ``$HOME/.zshrc``):

.. code-block:: console

   $ dir=$(pip show mstool | grep Location | awk '{print $2}')
   $ echo "mstoolpath=$dir/mstool" >> yourshellconfig


Option 3. Manual
----------------

mstool can be installed manually from the source code. After installing openMM with conda, Install the following dependencies:

* python (>=3.7)
* cython
* pandas
* numpy
* matplotlib

Once you install the dependencies, download the source code from GitHub and compile the cython code. 

.. code-block:: console

   $ cd ~
   $ git clone git@github.com:ksy141/mstool.git
   $ cd mstool/src/mstool/lib
   $ python setup.py build_ext --inplace &> /dev/null

Finally, add mstool to your shell configuration: ``$HOME/.bashrc`` or ``$HOME/.zshrc``

.. code-block:: console

   $ echo "export PYTHONPATH=$HOME/mstool/src:\$PYTHONPATH" >> yourshellconfig
   $ echo "mstoolpath=$HOME/mstool/src/mstool" >> yourshellconfig


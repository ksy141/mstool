Installation
============

GitHub
------
To use the up-to-date version of mstool, please download the package from `github.com/ksy141/mstool <https://github.com/ksy141/mstool>`_. The following Python packages are required:

* python (>=3.7)
* openmm
* sqlite3
* cython
* pandas
* numpy

.. code-block:: console

   $ cd ~
   $ git clone git@github.com:ksy141/mstool.git
   $ cd mstool
   $ mstoolpath=$(realpath .)
   $ bash install.sh

**install.sh** only works when using a **bash** or **zsh** shell, and your shell configuration is located in the home directory (`~/.bashrc` or `~/.zshrc`), which will likely be the case. If this is not the case, the following jobs should be done manually:

* Building Cython code:

.. code-block:: console

   $ cd $mstoolpath/lib-python/mstool/lib
   $ python setup.py build_ext --inplace &> /dev/null


* Adding mstool to ``PYTHONPATH``:

.. code-block:: console

   $ echo "export PYTHONPATH=$mstoolpath/lib-python:\$PYTHONPATH" >> yourshellconfig



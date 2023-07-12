# mstool

Multiscale Simulation Tools for Mapping, Backmapping, and Loop Modeling.

## Dependencies
openmm, sqlite3, cython, pandas, numpy. To install openmm:
```
conda install -c conda-forge openmm
```

## Installation
Simply run the below for installation.
```
$SHELL install.sh
```

The two actions are included in `install.sh`.
1. Add `PYTHONPATH` to the environment (`~/.bashrc` and `~/.zshrc`) and `alias`.
   ```
   export PYTHONPATH=$mstoolpath/lib-python:\$PYTHONPATH
   alias  mstool=$mstoolpath/
   ```
2. Build a distance matrix, which uses cython.
   ```
   cd $mstool/lib-python/mstool/lib
   python setup.py build_ext --inplace
   ```

## Mapping scheme
`mstool` will automatically read the default mapping files:
`$mstool/mapping/martini.protein.dat` and `$mstool/mapping/martini.lipid.c36.dat`

These are tailored toward martini. Please note that only a few Martini molecules are present in these mapping files at this moment. You should make a new file if you want to add molecules or your coarse-grained structures are not at martini resolution. Under each RESI `resname`, define a coarse-grained bead in a square bracket ([]) and its corresponding atomistic atoms. If a molecule has any isomers (cis/trans/chiral), define them in a square bracket. 

## Forcefield
Currently, it only supports charmm36 forcefield. The names of atomistic atoms defined in the mapping schemes should match these defined in the forcefield. All of the standard charmm36 molecules are defined at `$mstool/FF/charmm36/charmm36.xml`, which will be autmoatically read.

## How to use
```
import mstool
mstool.Backmap('input.pdb',  'backmapped.pdb', mapping_add = [])
mstool.REM('backmapped.pdb', 'output.pdb', mapping_add = [])
```

`input.pdb` is a coarse-grained structure. This will randomly place atomistic atoms near their corresponding coarse-grained bead. If you want to add mapping schemes because your molecules are not defined in the default mapping schemes, use `mapping_add = ['your_new_mapping1.dat', 'your_new_mapping2.dat']`.
If you want to overwrite the default mapping schemes, use `mapping = ['your_new_mapping1.dat', 'your_new_mapping2.dat']`.

`REM` will minimize a system. This is an important step as the above step just randomly ungroups atoms. cis/trans/chiral information is read from the default mapping schemes. If you want to add mapping schemes for new molecules, use `mapping_add = ['your_new_mapping1.dat', 'your_new_mapping2.dat']`.

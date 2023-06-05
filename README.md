# mstool

Multiscale Simulation Tools for Mapping, Backmapping, and Loop Modeling.

## Dependencies
openmm, sqlite3, pandas, numpy. To install openmm
```
conda install -c conda-forge openmm
```

## Installation
No installation is required. Download the package and add PYTHONPATH.
```
cd $HOME
git clone https://github.com/ksy141/mstool.git
export PYTHONPATH="$HOME:$PYTHONPATH"
```

## Mapping scheme
`mstool` will automatically read default mapping files:
`$HOME/mstool/mapping/martini.protein.dat` and `$HOME/mstool/mapping/martini.lipid.c36.dat`

These are tailored toward martini. Not all Martini molecules are present. You should make a new file if you want to add molecules or your coarse-grained structures are not at martini resolution. Under each RESI `resname`, define a coarse-grained bead in a square bracket ([]) and its corresponding atomistic atoms. If a molecule has any isomers (cis/trans/chiral), define them in a square bracket. 

## Forcefield
Currently, it only supports charmm36 forcefield. The names of atomistic atoms defined in the mapping schemes should match these defined in the forcefield. All of the standard C36 molecules are defined at `$HOME/mstool/FF/charmm36/charmm36.xml`, which will autmoatically read.

## How to use
```
import mstool

# input.pdb is a coarse-grained structure.
# This will randomly place atomistic atoms near their corresponding coarse-grained bead.
# If you want to add mapping schemes because your molecules are not defined in the default mapping schemes,
# use mapping_add = ['your_new_mapping1.dat', 'your_new_mapping2.dat']
# If you want to overwrite the default mapping schemes, use mapping = ['your_new_mapping1.dat', 'your_new_mapping2.dat'].
mstool.Backmap('input.pdb', 'backmapped.pdb', mapping_add = [])

# This will minimize a system. This is an important step as the above step just randomly ungroups atoms.
# cis/trans/chiral information is read from the default mapping schemes. If you want to add mapping schemes for new molecules,
# use mapping_add = ['your_new_mapping1.dat', 'your_new_mapping2.dat']
mstool.REM('backmapped.pdb', 'output.pdb', mapping_add = [])
```



import mstool
import shutil
import os

u = mstool.Universe('lipid.pdb')

mstool.REM(protein='../protein/workdir/step9_final.pdb',
           structure='lipid.pdb', 
           out='final.pdb', 
           mapping_add='../mapping.dat',
           ff_add='../trio.xml', 
           pbc=False)


mstool.CheckStructure('final.pdb')

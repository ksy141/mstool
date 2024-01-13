# copy this script to somewhere else
# and then execute this script
# to save time, I reduced the number of NVT steps (nsteps=100)
# The default is nsteps=10000, which is 2 ps.

import mstool
import shutil
import os

### Example 1
os.makedirs('Example1_methane')
os.chdir('Example1_methane')
shutil.copyfile(mstool.ONECGBEADSTRUCTURE, './cg.pdb')

with open('mapping.dat', 'w') as W:
    W.write('''RESI ONE
[ CG1 ]
C1 H1 H2 H3 H4
''')

with open('ff.xml', 'w') as W:
    W.write('''<ForceField>
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
''')

mstool.Backmap('cg.pdb', mapping='mapping.dat', ff_add='ff.xml', nsteps=100)
os.chdir('../')

### Example 2
os.makedirs('Example2_ethane')
os.chdir('Example2_ethane')
u = mstool.Universe(mstool.ONECGBEADSTRUCTURE)
u.atoms.resname = 'ETHA'
u.write('cg.pdb')

with open('mapping.dat', 'w') as W:
    W.write('''RESI ETHA
[ CG1 ]
C1 H11 H12 H13
C2 H21 H22 H23
''')

mstool.Backmap('cg.pdb', mapping='mapping.dat', nsteps=100)
os.chdir('../')

### Example 3
os.makedirs('Example3_butene')
os.chdir('Example3_butene')
u = mstool.Universe(mstool.ONECGBEADSTRUCTURE)
u.atoms.resname = 'BTE2'
u.write('cg.pdb')

with open('mapping.dat', 'w') as W:
    W.write('''RESI BTE2
[ CG1 ]
C1 H11 H12 H13
C2 H21
C3 H31
C4 H41 H42 H43

[ trans ]
C1 C2 C3 C4
''')

mstool.Backmap('cg.pdb', mapping='mapping.dat', nsteps=100)
os.chdir('../')

### Example 4
os.makedirs('Example4_POPC')
os.chdir('Example4_POPC')
shutil.copyfile(mstool.POPCSTRUCTURE, './cg.pdb')
mstool.Backmap('cg.pdb', nsteps=100)
os.chdir('../')

### Example 5
os.makedirs('Example5_Sphere')
os.chdir('Example5_Sphere')
shutil.copyfile(mstool.MULTISTRUCTURE, './cg.pdb')
mstool.Backmap('cg.pdb', nsteps=100)
os.chdir('../')

### Example 6
os.makedirs('Example6_TRIO')
os.chdir('Example6_TRIO')
shutil.copyfile(mstool.TRIOSTRUCTURE, './cg.pdb')
shutil.copyfile(mstool.TRIOMAPPING, './mapping.dat')
shutil.copyfile(mstool.TRIOFF, './ff.xml')
mstool.Backmap('cg.pdb', mapping_add='mapping.dat', ff_add='ff.xml', nsteps=100)
os.chdir('../')

### Example 7
os.makedirs('Example7_MembraneProtein')
os.chdir('Example7_MembraneProtein')
shutil.copyfile(mstool.MPCG, './cg.pdb')
shutil.copyfile(mstool.MPAA, './protein_AA.pdb')
u = mstool.Universe('cg.pdb')
non_protein_bA = ~u.atoms.resname.isin(mstool.three2one.keys())
non_protein = mstool.Universe(data=u.atoms[non_protein_bA])
non_protein.dimensions = u.dimensions
non_protein.write('cg_nonprotein.pdb')
mstool.Backmap(AA='protein_AA.pdb', structure='cg_nonprotein.pdb', nsteps=100)
os.chdir('../')



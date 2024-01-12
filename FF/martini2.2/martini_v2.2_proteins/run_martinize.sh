#!/bin/bash

for aa in CYS ASP SER GLN LYS ILE PRO THR PHE ASN GLY HIS LEU ARG TRP ALA VAL GLU TYR MET HIP HSD HSE HID HIE; do
    rm -rf $aa
    mkdir  $aa

    cat > ${aa}.py << EOF
import MDAnalysis as mda
u = mda.Universe('protein.gro')
ag = u.select_atoms('resname ${aa}')
residue = ag.residues[0]
residue.atoms.write('${aa}/input.pdb')
EOF

    python ${aa}.py
    rm -rf ${aa}.py

    python2 martinize_v2.6.py -f ${aa}/input.pdb -o topol.top -x input_cg.pdb -v -ss C -name ${aa} -nt
    sed -i "s/${aa}_X/${aa}/" ${aa}_X.itp
    mv ${aa}_X.itp ${aa}/${aa}.itp
    mv input_cg.pdb ${aa}

done

cat */*.itp > proteins.itp
rm topol.top


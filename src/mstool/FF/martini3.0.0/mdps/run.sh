#!/bin/bash

gmx grompp -f step6.0.mdp -c step3.pdb -o step6.0.tpr -maxwarn 5
gmx mdrun -deffnm step6.0 -v

for i in {1..6}; do
gmx grompp -f step6.${i}.mdp -c step6.$((i-1)).gro -o step6.${i}.tpr -maxwarn 5
gmx mdrun -deffnm step6.${i} -v
done

#gmx grompp -f step7.mdp -c step6.6.gro -o step7.tpr -maxwarn 5
#gmx mdrun -deffnm step7 -v


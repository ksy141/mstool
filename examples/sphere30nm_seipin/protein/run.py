import MDAnalysis as mda
from   MDAnalysis.analysis.align import alignto

cg  = mda.Universe('input_cg_prot.pdb')

ags = []
for chain in 'ABCDEFGHIJK':
    aa = mda.Universe('input_aa_prot.pdb')
    alignto(aa, cg, select = {'mobile':    f'chainID {chain} and name CA and resid 37 41 45 49 53',
                              'reference': f'chainID {chain} and name CA and resid 25 29 33 37 41'})
    ags.append(aa.select_atoms(f'chainID {chain} and resid 23 - 49'))

for chain in 'ABCDEFGHIJK':
    aa = mda.Universe('input_aa_prot.pdb')
    alignto(aa, cg, select = {'mobile':    f'chainID {chain} and name CA and resid 233 237 241 245 249',
                              'reference': f'chainID {chain} and name CA and resid 245 249 253 257 261'})
    ags.append(aa.select_atoms(f'chainID {chain} and resid 234 - 265'))


final = mda.Merge(*ags)
final.atoms.write('input_aa_prot_aligned.pdb')

body = ' '.join(cg.select_atoms('chainID A and resid 60-219').residues.resids.astype(str))
alignto(aa, cg, select='name CA and resid ' + body)
aa.select_atoms('resid 60-219').write('input_aa_body_aligned.pdb')
   

import msprot
import pandas as pd
msprot.Map(     'input_aa_prot_aligned.pdb',         'input_aa_prot_aligned_map.pdb')
msprot.Backmap( 'input_aa_prot_aligned_map.pdb',     'input_aa_prot_aligned_backmap.pdb')
msprot.REM(     'input_aa_prot_aligned_backmap.pdb', 'input_aa_prot_aligned_backmap_rem.pdb', pbc=False)


body         = msprot.Universe('input_aa_body_aligned.pdb')
termini      = msprot.Universe('input_aa_prot_aligned_backmap_rem.pdb')
u            = msprot.Universe(data = pd.concat([body.atoms, termini.atoms], ignore_index=True))
u.dimensions = msprot.Universe('../ext_306.pdb').dimensions
u.cell       = msprot.Universe('../ext_306.pdb').cell
u.sort()
u.write('input_aa_prot_final.pdb')

msprot.LoopModeler(protein = 'input_aa_prot_final.pdb',
                   fasta   = 'seipin.fasta')



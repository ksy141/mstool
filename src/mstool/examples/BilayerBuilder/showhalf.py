from chimerax.atomic import all_atoms
atoms = all_atoms(session)
bA = atoms.scene_coords[:,1] < 0
atoms[~bA].residues.atoms.displays = False
atoms[bA].residues.atoms.displays = True


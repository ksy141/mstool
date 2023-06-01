import MDAnalysis as mda
from   MDAnalysis.analysis.dihedrals import Dihedral
import numpy as np

u = {'gromacs':     mda.Universe('../relaxed.pdb', '../gromacs/input.whole.xtc'),
     'openmm':      mda.Universe('../relaxed.pdb', '../openmm/input.dcd'),
     'desmond_MTK': mda.Universe('../desmond_MTK/trj.pdb', '../desmond_MTK/trj.xtc'),
     'desmond_MC':  mda.Universe('../desmond_MC/trj.pdb', '../desmond_MC/trj.xtc')}

#u = {'gromacs':     mda.Universe('../relaxed.pdb', '../gromacs/input.whole.xtc'),
#     'openmm':      mda.Universe('../relaxed.pdb', '../openmm/input.dcd')}


def VOL(u):
    data = []
    for ts in u.trajectory[100:]:
        xyz = u.dimensions[0] * u.dimensions[1] * u.dimensions[2]
        data.append(xyz / u.residues.n_residues)
    return data


def angle(u, resname, atom1, atom2, atom3):
    data = []
    ag1  = u.select_atoms(f'resname {resname} and name ' + atom1)
    ag2  = u.select_atoms(f'resname {resname} and name ' + atom2)
    ag3  = u.select_atoms(f'resname {resname} and name ' + atom3)

    for ts in u.trajectory[100:]:
        dr1 = ag1.positions - ag2.positions
        dr2 = ag3.positions - ag2.positions

        dr1n = np.linalg.norm(dr1, axis=-1)
        dr2n = np.linalg.norm(dr2, axis=-1)
        cos  = np.sum(dr1 * dr2, axis=-1) / dr1n / dr2n
        cos  = cos[(cos < 1.0) & (-1.0 < cos)]
        data.extend(np.arccos(cos) * 180.0/3.141592)

    return data

def dihedral(u, resname):
    ags = []
    atoms = u.select_atoms(f'resname {resname}')
    for i in range(0, atoms.n_atoms, 4):
        ags.append(atoms[i: i  + 4])
    
    d = Dihedral(ags).run(start=100)
    results = d.results.angles.flatten()
    return results[~np.isnan(results)]


def bond(u, resname, atom1, atom2):
    data = []
    ag1  = u.select_atoms(f'resname {resname} and name ' + atom1)
    ag2  = u.select_atoms(f'resname {resname} and name ' + atom2)

    for ts in u.trajectory[100:]:
        dr1 = ag1.positions - ag2.positions
        dr1n = np.linalg.norm(dr1, axis=-1)
        data.extend(dr1n)

    return data


def plot_bond(u, resname, atom1, atom2, r0):
    fig, ax = plt.subplots(figsize=(3.5,2.8))
    for key, value in u.items():
        #ax.hist(bond(value, atom1, atom2), bins=50, alpha=0.4, label=key, density=True)
        counts, bins = np.histogram(bond(value, resname, atom1, atom2), bins=50, density=True)
        ax.plot(0.5 * (bins[1:] + bins[:-1]), counts, label=key)
        ax.fill_between(0.5 * (bins[1:] + bins[:-1]), counts, 0, alpha=.3)

    ax.axvline(r0, c='black')
    ax.legend(frameon=False, fontsize=8)
    ax.set_xlabel('r (A)')
    ax.set_title(f'bond_{resname}_{atom1}_{atom2}')
    fig.tight_layout()
    fig.savefig(f'analysis_bond_{resname}_{atom1}_{atom2}.pdf')


def plot_angle(u, resname, atom1, atom2, atom3, t0):
    fig, ax = plt.subplots(figsize=(3.5,2.8))
    for key, value in u.items():
        #ax.hist(angle(value, atom1, atom2, atom3), bins=50, alpha=0.4, label=key, density=True)
        counts, bins = np.histogram(angle(value, resname, atom1, atom2, atom3), bins=50, density=True)
        ax.plot(0.5 * (bins[1:] + bins[:-1]), counts, label=key)
        ax.fill_between(0.5 * (bins[1:] + bins[:-1]), counts, 0, alpha=.3)

    ax.axvline(t0, c='black')
    ax.legend(frameon=False, fontsize=8)
    ax.set_xlabel('t (degree)')
    ax.set_title(f'angle_{resname}_{atom1}_{atom2}_{atom3}')
    fig.tight_layout()
    fig.savefig(f'analysis_angle_{resname}_{atom1}_{atom2}_{atom3}.pdf')


def plot_dihedral(u, resname, t0):
    fig, ax = plt.subplots(figsize=(3.5,2.8))
    for key, value in u.items():
        #ax.hist(dihedral(value, resname), bins=50, alpha=0.4, label=key, density=True)
        counts, bins = np.histogram(dihedral(value, resname), bins=50, density=True)
        ax.plot(0.5 * (bins[1:] + bins[:-1]), counts, label=key)
        ax.fill_between(0.5 * (bins[1:] + bins[:-1]), counts, 0, alpha=.3)

    ax.axvline(t0, c='black')
    ax.legend(frameon=False, fontsize=8)
    ax.set_xlabel('t (degree)')
    ax.set_title(f'dihedral_{resname}')
    fig.tight_layout()
    fig.savefig(f'analysis_dihedral_{resname}.pdf')




import matplotlib.pyplot as plt

# VOL
fig, ax = plt.subplots(figsize=(3.5,2.8))
for key, value in u.items():
    ax.hist(VOL(value), bins=50, alpha=0.4, label=key, density=True)
ax.set_xlabel('v (A^3)')
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig('analysis_vol.pdf')

# BOND
plot_bond(u, 'TP', 'GL1', 'PO4', 4.7)
plot_bond(u, 'TP', 'C1A', 'GL1', 4.7)
plot_bond(u, 'TP', 'C1A', 'C2A', 4.7)
plot_bond(u, 'TI', 'C1A', 'C2A', 4.7)
plot_bond(u, 'TI', 'C1A', 'C2A', 4.7)
plot_bond(u, 'TI', 'C1A', 'C2A', 4.7)

# ANGLE
plot_angle(u, 'TP', 'PO4', 'GL1', 'C1A', 120)
plot_angle(u, 'TP', 'GL1', 'C1A', 'C2A', 120)
plot_angle(u, 'TI', 'PO4', 'GL1', 'C1A', 120)
plot_angle(u, 'TI', 'GL1', 'C1A', 'C2A', 120)

# DIHEDRAL
plot_dihedral(u, 'TP', 125)
plot_dihedral(u, 'TI', 125)



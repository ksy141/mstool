import MDAnalysis as mda
import numpy as np

#u = {'gromacs':    mda.Universe('../input.pdb', '../gromacs/input.whole.xtc'),
#     'gromacs_cf': mda.Universe('../input.pdb', '../gromacs_cutoff/input.whole.xtc'),
#     'gromacsh':   mda.Universe('../input.pdb', '../gromacsh/input.whole.xtc'),
#     'gromacshm':  mda.Universe('../input.pdb', '../gromacshm/input.whole.xtc'),
#     'openmm':     mda.Universe('../input.pdb', '../openmm/input.dcd')}

u = {'gromacs':    mda.Universe('../input.pdb', '../gromacs/input.whole.xtc'),
     'openmm':     mda.Universe('../input.pdb', '../openmm/input.dcd')}

def APL(u):
    data = []
    for ts in u.trajectory:
        xy = u.dimensions[0] * u.dimensions[1] / 16
        data.append(xy)
    return data


def angle(u, atom1, atom2, atom3):
    data = []
    ag1  = u.select_atoms('name ' + atom1)
    ag2  = u.select_atoms('name ' + atom2)
    ag3  = u.select_atoms('name ' + atom3)

    for ts in u.trajectory:
        dr1 = ag1.positions - ag2.positions
        dr2 = ag3.positions - ag1.positions

        dr1n = np.linalg.norm(dr1, axis=-1)
        dr2n = np.linalg.norm(dr2, axis=-1)
        cos  = np.sum(dr1 * dr2, axis=-1) / dr1n / dr2n
        data.extend(np.arccos(cos) * 180.0/3.141592)

    return data


def bond(u, atom1, atom2):
    data = []
    ag1  = u.select_atoms('name ' + atom1)
    ag2  = u.select_atoms('name ' + atom2)

    for ts in u.trajectory:
        dr1 = ag1.positions - ag2.positions
        dr1n = np.linalg.norm(dr1, axis=-1)
        data.extend(dr1n)

    return data


def plot_bond(u, atom1, atom2, r0):
    fig, ax = plt.subplots(figsize=(3.5,2.8))
    for key, value in u.items():
        ax.hist(bond(value, atom1, atom2), bins=50, alpha=0.4, label=key)
    ax.axvline(r0, c='black')
    ax.legend(frameon=False)
    ax.set_xlabel('r (A)')
    fig.tight_layout()
    fig.savefig('analysis_bond_{:s}_{:s}.pdf'.format(atom1, atom2))


def plot_angle(u, atom1, atom2, atom3, t0):
    fig, ax = plt.subplots(figsize=(3.5,2.8))
    for key, value in u.items():
        ax.hist(angle(value, atom1, atom2, atom3), bins=50, alpha=0.4, label=key)
    ax.axvline(t0, c='black')
    ax.legend(frameon=False)
    ax.set_xlabel('t (degree)')
    fig.tight_layout()
    fig.savefig('analysis_angle_{:s}_{:s}_{:s}.pdf'.format(atom1, atom2, atom3))



import matplotlib.pyplot as plt

# APL
fig, ax = plt.subplots(figsize=(3.5,2.8))
for key, value in u.items():
    ax.hist(APL(value), bins=50, alpha=0.4, label=key)
ax.set_xlabel('APL (A^2)')
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig('analysis_apl.pdf')

# BOND
plot_bond(u, 'NC3', 'PO4', 4.7)
plot_bond(u, 'PO4', 'GL1', 4.7)
plot_bond(u, 'GL1', 'GL2', 3.7)
plot_bond(u, 'GL1', 'C1A', 4.7)
plot_bond(u, 'C1A', 'D2A', 4.7)
plot_bond(u, 'D2A', 'C3A', 4.7)
plot_bond(u, 'C3A', 'C4A', 4.7)


# ANGLE
plot_angle(u, 'PO4', 'GL1', 'GL2', 120)
plot_angle(u, 'PO4', 'GL1', 'C1A', 180)
plot_angle(u, 'GL1', 'C1A', 'D2A', 180)
plot_angle(u, 'C1A', 'D2A', 'C3A', 120)
plot_angle(u, 'D2A', 'C3A', 'C4A', 180)
plot_angle(u, 'GL2', 'C1B', 'C2B', 180)
plot_angle(u, 'GL2', 'C1B', 'C2B', 180)
plot_angle(u, 'C1B', 'C2B', 'C3B', 180)
plot_angle(u, 'C2B', 'C3B', 'C4B', 180)



import MDAnalysis as mda
import numpy as np

#u = {'gromacs':     mda.Universe('../input.pdb', '../gromacs/input.whole.xtc'),
#     'openmm':      mda.Universe('../input.pdb', '../openmm/input.dcd'),
#     'desmond_MTK': mda.Universe('../desmond_MTK/trj.pdb', '../desmond_MTK/trj.xtc'),
#     'desmond_MC':  mda.Universe('../desmond_MC/trj.pdb', '../desmond_MC/trj.xtc')}

u = {'gromacs':     mda.Universe('../input.pdb', '../gromacs/input.whole.xtc'),
     'openmm':      mda.Universe('../input.pdb', '../openmm/input.dcd')}


def VOL(u):
    data = []
    for ts in u.trajectory[100:]:
        xyz = u.dimensions[0] * u.dimensions[1] * u.dimensions[2]
        data.append(xyz / u.residues.n_residues)
    return data


def angle(u, atom1, atom2, atom3):
    data = []
    ag1  = u.select_atoms('name ' + atom1)
    ag2  = u.select_atoms('name ' + atom2)
    ag3  = u.select_atoms('name ' + atom3)

    for ts in u.trajectory[100:]:
        dr1 = ag1.positions - ag2.positions
        dr2 = ag3.positions - ag2.positions

        dr1n = np.linalg.norm(dr1, axis=-1)
        dr2n = np.linalg.norm(dr2, axis=-1)
        cos  = np.sum(dr1 * dr2, axis=-1) / dr1n / dr2n
        data.extend(np.arccos(cos) * 180.0/3.141592)

    return data


def bond(u, atom1, atom2):
    data = []
    ag1  = u.select_atoms('name ' + atom1)
    ag2  = u.select_atoms('name ' + atom2)

    for ts in u.trajectory[100:]:
        dr1 = ag1.positions - ag2.positions
        dr1n = np.linalg.norm(dr1, axis=-1)
        data.extend(dr1n)

    return data


def plot_bond(u, atom1, atom2, r0):
    fig, ax = plt.subplots(figsize=(3.5,2.8))
    for key, value in u.items():
        #ax.hist(bond(value, atom1, atom2), bins=50, alpha=0.4, label=key, density=True)
        counts, bins = np.histogram(bond(value, atom1, atom2), bins=50, density=True)
        ax.plot(0.5 * (bins[1:] + bins[:-1]), counts, label=key)
        ax.fill_between(0.5 * (bins[1:] + bins[:-1]), counts, 0, alpha=.3)

    ax.axvline(r0, c='black')
    ax.legend(frameon=False)
    ax.set_xlabel('r (A)')
    fig.tight_layout()
    fig.savefig('analysis_bond_{:s}_{:s}.pdf'.format(atom1, atom2))


def plot_angle(u, atom1, atom2, atom3, t0):
    fig, ax = plt.subplots(figsize=(3.5,2.8))
    for key, value in u.items():
        #ax.hist(angle(value, atom1, atom2, atom3), bins=50, alpha=0.4, label=key, density=True)
        counts, bins = np.histogram(angle(value, atom1, atom2, atom3), bins=50, density=True)
        ax.plot(0.5 * (bins[1:] + bins[:-1]), counts, label=key)
        ax.fill_between(0.5 * (bins[1:] + bins[:-1]), counts, 0, alpha=.3)

    ax.axvline(t0, c='black')
    ax.legend(frameon=False)
    ax.set_xlabel('t (degree)')
    fig.tight_layout()
    fig.savefig('analysis_angle_{:s}_{:s}_{:s}.pdf'.format(atom1, atom2, atom3))



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
plot_bond(u, 'GL1', 'PO4', 4.7)
plot_bond(u, 'C1A', 'GL1', 4.7)

# ANGLE
plot_angle(u, 'PO4', 'GL1', 'C1A', 120)


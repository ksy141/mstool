import MDAnalysis as mda
import matplotlib.pyplot as plt
import numpy as np
from   MDAnalysis.analysis.rdf import InterRDF

u = {'gromacs': mda.Universe('../gromacs/W.pdb', '../gromacs/input.xtc'),
     'openmm':  mda.Universe('../openmm/W.pdb', '../openmm/W.dcd'),
     'desmond_MTK': mda.Universe('../desmond_MTK/trj.pdb', '../desmond_MTK/trj.xtc'),
     'desmond_MC': mda.Universe('../desmond_MC/trj.pdb', '../desmond_MC/trj.xtc')}


def VOL(u):
    data = []
    for ts in u.trajectory[100:]:
        xyz = u.dimensions[0] * u.dimensions[1] * u.dimensions[2]
        data.append(xyz / u.atoms.n_atoms)
    print(np.average(data))
    return data

def RDF(u):
    rdf = InterRDF(u.atoms, u.atoms, range=(0.0,20.0))
    rdf.run(start=100)
    return rdf.results.bins[1:], rdf.results.rdf[1:]


def plot_vol(u):
    fig, ax = plt.subplots(figsize=(3.5,2.8))
    for key, value in u.items():
        ax.hist(VOL(value), bins=50, alpha=0.4, label=key)
    ax.legend(frameon=False)
    ax.set_xlabel('v (A^3)')
    fig.tight_layout()
    fig.savefig('analysis_vol.pdf')

def plot_rdf(u):
    fig, ax = plt.subplots(figsize=(3.5,2.8))
    for key, value in u.items():
        bins, r = RDF(value)
        ax.plot(bins, r, label=key)
    ax.legend(frameon=False)
    ax.set_xlabel('r (A)')
    ax.set_ylabel('g(r)')
    fig.tight_layout()
    fig.savefig('analysis_rdf.pdf')
        
plot_vol(u)
plot_rdf(u)


"""
openmm2psf.py
Dump a CHARMM/X-PLOR-style PSF from an OpenMM System **without ParmEd**.

Key features
------------
✓ Preserves force‑field atom‑types by re‑parsing the ForceField XML (Option 1)  
✓ Writes bonds, angles, dihedrals, **impropers**, and **CMAP cross‑terms**  
✓ No external dependencies beyond OpenMM

Example
-------
from openmm.app import *
from openmm import *
from openmm.unit import *

pdb  = PDBFile("protein.pdb")
ff   = ForceField("charmm36.xml", "tip3p.xml")
system = ff.createSystem(pdb.topology, nonbondedMethod=app.PME)
write_psf(system, pdb.topology, ff, "protein.psf")
"""

import datetime as _dt
from openmm import (
    HarmonicBondForce,
    HarmonicAngleForce,
    PeriodicTorsionForce,
    CustomTorsionForce,
    CMAPTorsionForce,
    unit as _u,
)
from ..core.readxml import ReadXML
from openmm.app import DesmondDMSFile, PDBFile, ForceField

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _wrap(seq, per_line):
    """Yield successive chunks with *per_line* tuples each."""
    for i in range(0, len(seq), per_line):
        yield seq[i : i + per_line]


def _block(fh, label, items, per_line):
    fh.write(f"\n{len(items):8d} {label}\n")
    for chunk in _wrap(items, per_line):
        fh.write("".join(f"{x:8d}" for tup in chunk for x in tup) + "\n")


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def WritePSF(structure, out, ff=None, ff_add=None):
    ff = ff or []
    ff_add = ff_add or []

    xml = ReadXML(ff=ff, ff_add=ff_add)
    forcefield = ForceField(*xml.ff)

    if structure.endswith(".pdb"):
        pdb = PDBFile(structure)
    elif structure.endswith(".dms"):
        pdb = DesmondDMSFile(structure)
    else:
        raise ValueError("Unsupported structure format. Expected .pdb or .dms")

    topology = pdb.topology

    system = forcefield.createSystem(topology)
    data = forcefield._SystemData(topology)
    templateForResidue = forcefield._matchAllResiduesToTemplates(
        data, topology, [False] * topology.getNumResidues(), False
    )

    # --- 2.  Collect atom block --------------------------------------------
    atoms, index_map = [], {}
    charges = {}

    # Grab charges from the first NonbondedForce
    from openmm import NonbondedForce

    for f in system.getForces():
        if isinstance(f, NonbondedForce):
            for i in range(system.getNumParticles()):
                q, *_ = f.getParticleParameters(i)
                charges[i] = q.value_in_unit(_u.elementary_charge)
            break

    aid = 1
    default_segid = "SYS"
    if topology is not None:
        for chain in topology.chains():
            segid = chain.id or default_segid
            for res in chain.residues():
                resid = int(res.id or 1)
                resname = (res.name or "RES")[:4]
                for atom in res.atoms():
                    mass = system.getParticleMass(atom.index).value_in_unit(_u.dalton)
                    charge = charges.get(atom.index, 0.0)

                    atype = data.atomType[atom]
                    mass = forcefield._atomTypes[atype].mass

                    atoms.append(
                        (aid, segid, resid, resname, atom.name[:4], atype, charge, mass)
                    )
                    index_map[atom.index] = aid
                    aid += 1
    else:  # no topology – invent names
        for omm_idx in range(system.getNumParticles()):
            mass = system.getParticleMass(omm_idx).value_in_unit(_u.dalton)
            charge = charges.get(omm_idx, 0.0)
            name = f"A{omm_idx+1}"
            atoms.append((
                omm_idx + 1,
                default_segid,
                1,
                "RES",
                name,
                name,
                charge,
                mass,
            ))
            index_map[omm_idx] = omm_idx + 1

    natom = len(atoms)

    # --- 3.  Connectivity ---------------------------------------------------
    bonds, angles, diheds, impropers, cmaps = [], [], [], [], []
    bonds = [[index_map[bond[0].index], index_map[bond[1].index]] for bond in topology.bonds()]

    for f in system.getForces():
        #if isinstance(f, HarmonicBondForce):
        #    for i in range(f.getNumBonds()):
        #        a, b, *_ = f.getBondParameters(i)
        #        bonds.append((index_map[a], index_map[b]))

        if isinstance(f, HarmonicAngleForce):
            for i in range(f.getNumAngles()):
                a, b, c, *_ = f.getAngleParameters(i)
                angles.append((index_map[a], index_map[b], index_map[c]))

        elif isinstance(f, PeriodicTorsionForce):
            # CHARMM proper dihedrals come in here; impropers usually DON'T
            for i in range(f.getNumTorsions()):
                a, b, c, d, *_ = f.getTorsionParameters(i)
                diheds.append((index_map[a], index_map[b], index_map[c], index_map[d]))

        elif isinstance(f, CustomTorsionForce):
            # ForceField builds CHARMM-style harmonic impropers with a
            # quadratic form: k*(phi - phi0)^2   → look for "k" and "phi0"
            expr = f.getEnergyFunction().replace(" ", "").lower()
            is_improper = ("phi0" in expr or "theta0" in expr) and "k" in expr
            if is_improper:
                for i in range(f.getNumTorsions()):
                    a, b, c, d, *_ = f.getTorsionParameters(i)
                    impropers.append(
                        (
                            index_map[a],
                            index_map[b],
                            index_map[c],
                            index_map[d],
                        )
                    )
            else:
                # some users encode extra proper dihedrals via CustomTorsionForce
                for i in range(f.getNumTorsions()):
                    a, b, c, d, *_ = f.getTorsionParameters(i)
                    diheds.append(
                        (
                            index_map[a],
                            index_map[b],
                            index_map[c],
                            index_map[d],
                        )
                    )

        elif isinstance(f, CMAPTorsionForce):
            for i in range(f.getNumTorsions()):
                _map,  a, b, c, d, e, f_, g, h = f.getTorsionParameters(i)
                cmaps.append(
                    (
                        index_map[a],
                        index_map[b],
                        index_map[c],
                        index_map[d],
                        index_map[e],
                        index_map[f_],
                        index_map[g],
                        index_map[h],
                    )
                )

    # --- 4.  Write PSF ------------------------------------------------------
    flags = ["PSF", "EXT"]

    # CMAP: include flag if any CMAP terms present
    if cmaps:
        flags.append("CMAP")
    # XPLOR: default text-based types
    flags.append("XPLOR")

    with open(out, "w") as fh:
        fh.write(" ".join(flags) + "\n\n")
        fh.write(f"{2:8d} !NTITLE\n")
        fh.write(
            f" REMARKS Generated by openmm2psf on {_dt.datetime.now():%Y-%m-%d %H:%M:%S}\n"
        )
        fh.write(" REMARKS Atom types from ForceField XML; impropers & CMAP terms included\n")

        # Atom block (EXT format)
        fh.write(f"\n{natom:8d} !NATOM\n")
        for idx, seg, resid, rname, aname, atype, q, m in atoms:
            fh.write(
                f"{idx:8d} {seg:<8s} {resid:8d} {rname:<8s} {aname:<8s} "
                f"{atype:<8s}{q:10.6f}{m:13.4f}     0\n"
            )

        _block(fh, "!NBOND: bonds", bonds, 4)
        _block(fh, "!NTHETA: angles", angles, 3)
        _block(fh, "!NPHI: dihedrals", diheds, 2)
        _block(fh, "!NIMPHI: impropers", impropers, 2)
        _block(fh, "!NCRTERM: cross-terms", cmaps, 1)

        fh.write("\n       0 !NDON\n\n")
        fh.write("       0 !NACC\n\n\n")

    print(
        f"Wrote {out}: {natom} atoms / "
        f"{len(bonds)} bonds / {len(angles)} angles / "
        f"{len(diheds)} dihedrals / {len(impropers)} impropers / "
        f"{len(cmaps)} CMAP terms"
    )

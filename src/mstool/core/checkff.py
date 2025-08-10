from   ..utils.atomic_data import anum_mass_vdw_from_name
import networkx as nx

def mol2nx(mol, resname=None):
    """Convert a molecule to a networkx graph.
    >>> molnx = mstool.mol2nx(mstool.ReadXML().RESI['POPC'], resname='POPC')
    """
    G = nx.Graph()
    for nn, tt, qq in zip(mol['names'], mol['types'], mol['charges']):
        upper_name = tt.upper()
        if upper_name in anum_mass_vdw_from_name.keys():
            anum, mass, vdw = anum_mass_vdw_from_name[upper_name]
        elif len(upper_name) > 1.5 and upper_name[:2] in anum_mass_vdw_from_name.keys():
            anum, mass, vdw = anum_mass_vdw_from_name[upper_name[:2]]
        elif upper_name[0] in anum_mass_vdw_from_name.keys():
            anum, mass, vdw = anum_mass_vdw_from_name[upper_name[0]]
        else:
            raise ValueError(f"Unknown atom name: {upper_name}")

        G.add_node(str(nn), type=str(tt), charge=float(qq), atomic_number=int(anum))

    for namei, namej in mol['bonds']:
        G.add_edge(str(namei), str(namej))

    if resname is not None:
        G.resname = resname

    return G

def node_match(n1, n2):
    return n1['atomic_number'] == n2['atomic_number']

def CheckFF(ff1, ff2):
    """Check if there are any duplicate residues in two force fields."""
    ff1nx = []
    for resname in ff1.RESI.keys():
        mol = ff1.RESI[resname]
        if len(mol['names']) < 4:
            continue
        ff1nx.append(mol2nx(mol, resname))

    ff2nx = []
    for resname in ff2.RESI.keys():
        mol = ff2.RESI[resname]
        if len(mol['names']) < 4:
            continue
        ff2nx.append(mol2nx(mol, resname))

    count = 0
    for f1nx in ff1nx:
        for f2nx in ff2nx:
            if nx.is_isomorphic(f1nx, f2nx, node_match=node_match):
                print(f"Duplicate residue found: {f1nx.resname} {f2nx.resname} in both force fields.")
                count += 1
    
    return count






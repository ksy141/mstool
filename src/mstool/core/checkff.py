from   ..utils.atomic_data import anum_mass_vdw_from_name
import networkx as nx
from rdkit import Chem

def nx2rdkit(graph: nx.Graph) -> Chem.Mol:
    """
    Converts a graph with atomic_number node attributes to an RDKit Mol object.
    All bonds are assumed to be single bonds.
    
    Parameters:
        graph (nx.Graph): Nodes must have 'atomic_number' attribute.
        
    Returns:
        Chem.Mol: RDKit molecule object.
    """
    mol = Chem.RWMol()
    atom_idx_map = {}

    # Add atoms
    for node in graph.nodes:
        atomic_number = graph.nodes[node].get('atomic_number', None)
        if atomic_number is None:
            raise ValueError(f"Node {node} missing 'atomic_number' attribute.")
        atom = Chem.Atom(atomic_number)
        idx = mol.AddAtom(atom)
        atom_idx_map[node] = idx

    # Add bonds (assume SINGLE)
    for u, v in graph.edges:
        mol.AddBond(atom_idx_map[u], atom_idx_map[v], Chem.BondType.SINGLE)

    # Finalize Mol
    final_mol = mol.GetMol()
    Chem.SanitizeMol(final_mol)  # raises exception if invalid
    return final_mol


def mol2nx(mol, resname=None):
    """Convert a molecule to a networkx graph.
    >>> molnx = mstool.mol2nx(mstool.ReadXML().RESI['POPC'], resname='POPC')
    """
    G = nx.Graph()
    for nn, tt, qq, number in zip(mol['names'], mol['types'], mol['charges'], mol['numbers']):
        if number == 0:
            upper_name = tt.upper()
            if upper_name in anum_mass_vdw_from_name.keys():
                anum, mass, vdw = anum_mass_vdw_from_name[upper_name]
            elif len(upper_name) > 1.5 and upper_name[:2] in anum_mass_vdw_from_name.keys():
                anum, mass, vdw = anum_mass_vdw_from_name[upper_name[:2]]
            elif upper_name[0] in anum_mass_vdw_from_name.keys():
                anum, mass, vdw = anum_mass_vdw_from_name[upper_name[0]]
            else:
                raise ValueError(f"Unknown atom name: {upper_name}")
        else:
            anum = number

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


def CheckMol(nxmol1, nxmol2):
    """Match two nx mols and print out the difference"""
    if not nx.is_isomorphic(nxmol1, nxmol2, node_match=node_match):
        print("two molecules are not the same")
        return
    
    GM = nx.isomorphism.GraphMatcher(nxmol1, nxmol2, node_match=node_match)
    for mapping in GM.subgraph_isomorphisms_iter():
        for nodeid1, nodeid2 in mapping.items():
            data1 = nxmol1.nodes[nodeid1]
            data2 = nxmol2.nodes[nodeid2]

            if data1['charge'] != data2['charge'] or data1['type'] != data2['type']:
                print(f"{nodeid1}: charge={data1['charge']}, type={data1['type']} "
                      f"<-> {nodeid2}: charge={data2['charge']}, type={data2['type']}")
        break


from rdkit import Chem
import networkx as nx
import os
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib.pyplot as plt
pwd = os.path.dirname(os.path.realpath(__file__))

def VisG(G, attrs=True):
    # Layout positions
    pos = nx.spring_layout(G)

    # Custom labels
    if attrs:
        labels = {
            n: f"{n} ({G.nodes[n]['type']}, q={G.nodes[n]['charge']})"
            for n in G.nodes
        }
    else:
        labels = {n: f"{n}" for n in G.nodes}

    # Define element color map
    color_map = {
        1:  'lightgray',   # Hydrogen
        6:  'gray',        # Carbon
        7:  'blue',        # Nitrogen
        8:  'red',         # Oxygen
        15: 'orange',      # Phosphorus
        16: 'yellow'       # Sulfur
    }

    # Assign colors and sizes
    node_colors = []
    node_sizes = []

    for n in G.nodes:
        atomic_number = G.nodes[n].get('atomic_number')
        color = color_map.get(atomic_number, 'green')  # Fallback color
        size = 300 + atomic_number * 50
        node_colors.append(color)
        node_sizes.append(size)

    # Edge widths based on bond order
    edge_widths = [G[u][v].get('order', 1) * 2 for u, v in G.edges()]

    # Draw
    nx.draw(
        G, pos,
        with_labels=False,
        node_color=node_colors,
        node_size=node_sizes,
        width=edge_widths
    )
    nx.draw_networkx_labels(G, pos, labels, font_size=9)

    plt.title("Graph with name, type, and charge per node")
    plt.axis('off')
    plt.show()


def DFS(G):
    nodes   = [node for node in G.nodes()]
    visited = set()
    queue   = [nodes[0]]

    while queue:
        current = queue.pop()
        if current not in visited:
            visited.add(current)
            for neighbor in G.neighbors(current):
                if neighbor not in visited:
                    queue.append(neighbor)

    return nodes, visited


def rdkit_to_nx(mol):
    G = nx.Graph()

    for atom in mol.GetAtoms():
        G.add_node(atom.GetIdx(), atomic_number=atom.GetAtomicNum())

    for bond in mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), order=bond.GetBondType())

    for idx, atom in G.nodes(data=True):
        count_C = 0; count_H = 0
        if atom['atomic_number'] == 6:
            neighbors = G.neighbors(idx)
            for neighbor in neighbors:
                if G.nodes[neighbor]['atomic_number'] == 1:
                    count_H += 1
                elif G.nodes[neighbor]['atomic_number'] == 6:
                    count_C += 1

            neighbors = G.neighbors(idx)
            if count_H == 3 and count_C == 1:
                G.nodes[idx]['type']   = 'CTL3'
                G.nodes[idx]['charge'] = -0.27
                G.nodes[idx]['z']      = -10
                for neighbor in neighbors:
                    if G.nodes[neighbor]['atomic_number'] == 1:
                        G.nodes[neighbor]['type']   = 'HAL3'
                        G.nodes[neighbor]['charge'] = 0.09
                        G.nodes[neighbor]['z']      = -10

            elif count_H == 2 and count_C == 2:
                G.nodes[idx]['type']   = 'CTL2'
                G.nodes[idx]['charge'] = -0.18
                G.nodes[idx]['z']      = -10
                for neighbor in neighbors:
                    if G.nodes[neighbor]['atomic_number'] == 1:
                        G.nodes[neighbor]['type']   = 'HAL2'
                        G.nodes[neighbor]['charge'] = 0.09
                        G.nodes[neighbor]['z']      = -10

            elif count_H == 1 and count_C == 2:
                G.nodes[idx]['type']   = 'CEL1'
                G.nodes[idx]['charge'] = -0.15
                G.nodes[idx]['z']      = -10
                for neighbor in neighbors:
                    if G.nodes[neighbor]['atomic_number'] == 1:
                        G.nodes[neighbor]['type']   = 'HEL1'
                        G.nodes[neighbor]['charge'] = 0.15
                        G.nodes[neighbor]['z']      = -10

    G.smiles = Chem.MolToSmiles(mol)
    return G

def node_match(n1, n2):
    return n1['atomic_number'] == n2['atomic_number']

def edge_match(e1, e2):
    return e1['order'] == e2['order']

def ReadBuildingBlock(ifiles):
    graphs = []

    for ifile in ifiles:
        root = ET.parse(ifile).getroot()

        for block in root.findall('BuildingBlock'):
            name = block.get("name") or "Unnamed"
            G = nx.Graph()

            for pos in block.findall("Pos"):
                if pos.get("z") is not None:
                    z = float(pos.get("z"))
                else:
                    z = 0.0

            for atom in block.findall("Atom"):
                atom_name = atom.get("name")
                G.add_node(atom_name,
                           charge=float(atom.get("charge", 0.0)),
                           type=atom.get("type"), z=z,
                           atomic_number=int(atom.get("number", 0)) if atom.get("number") else None)

            # Add edges based on Bonds
            for bond in block.findall("Bond"):
                atom1 = bond.get("atomName1")
                atom2 = bond.get("atomName2")
                order_str = bond.get("order")

                try:
                    order = int(order_str) if order_str is not None else 1
                except ValueError:
                    order = 1  # Fallback in case of malformed input

                G.add_edge(atom1, atom2, order=order)

            DFS_nodes, DFS_visited = DFS(G)
            if len(DFS_nodes) == len(DFS_visited):
                graphs.append(G)
            else:
                print(f"{name} has more than one fragment. {DFS_nodes}, {DFS_visited}")

    graphs.sort(key=lambda g: len(g.nodes))
    return graphs


def CisTrans(mol):
    cis   = []
    trans = []

    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE and bond.GetStereo() in (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ):
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()

            # Get substituents (only 1 neighbor that's not the double-bonded atom)
            begin_neighbors = [nbr.GetIdx() for nbr in begin_atom.GetNeighbors() if nbr.GetIdx() != end_atom.GetIdx()]
            end_neighbors = [nbr.GetIdx() for nbr in end_atom.GetNeighbors() if nbr.GetIdx() != begin_atom.GetIdx()]

            # Sort to ensure reproducibility
            begin_neighbors.sort()
            end_neighbors.sort()

            if len(begin_neighbors) != 1 or len(end_neighbors) != 1:
                print(f"Warning: Unexpected number of neighbors for E/Z bond {bond.GetIdx()}")
                continue

            i = begin_neighbors[0]  # atom1
            j = begin_atom.GetIdx() # atom2
            k = end_atom.GetIdx()   # atom3
            l = end_neighbors[0]    # atom4

            stereo = bond.GetStereo()

            if stereo == Chem.BondStereo.STEREOZ:
                cis.append([i, j, k, l])
            else:
                trans.append([i, j, k, l])

    return cis, trans


def Chiral(mol):
    chirals = []
    Chem.FindMolChiralCenters(mol, includeUnassigned=False)  # to update atom stereo tags
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)

    for center_idx, config in chiral_centers:
        center_atom = mol.GetAtomWithIdx(center_idx)
        neighbors = center_atom.GetNeighbors()

        if len(neighbors) != 4:
            print(f"Warning: atom {center_idx} is not bonded to 4 atoms.")
            continue

        hydrogen_idx = None
        heavy_atom_indices = []

        for nbr in neighbors:
            if nbr.GetAtomicNum() == 1:  # Hydrogen
                hydrogen_idx = nbr.GetIdx()
            else:
                heavy_atom_indices.append(nbr.GetIdx())

        if hydrogen_idx is None:
            # No hydrogen, fall back to default
            neighbor_indices = sorted([nbr.GetIdx() for nbr in neighbors])
            atom1, atom3, atom4, atom5 = neighbor_indices
        else:
            if len(heavy_atom_indices) != 3:
                print(f"Warning: unexpected number of heavy atoms at center {center_idx}")
                continue
            atom1 = hydrogen_idx
            atom3, atom4, atom5 = sorted(heavy_atom_indices)

        atom2 = center_idx  # chiral center
        chirals.append([atom1, atom2, atom3, atom4, atom5])

    return chirals


class NewLipid:
    """Create a new lipid type based on SMILES.
    >>> from mstool.membrane.newlipid import NewLipid
    >>> n = NewLipid('CCCCCCCC/C=C\CCCCCCCC(O[C@H](COC(CCCCCCCCCCCCCCC)=O)COP(OCC[N+](C)(C)C)([O-])=O)=O')
    >>> VisG(n.nxmol)
    """
    def __init__(self, smiles, bb=[], bb_add=[]):
        self.smiles    = smiles
        self.rdkitmol  = Chem.MolFromSmiles(smiles)
        self.rdkitmolH = Chem.AddHs(self.rdkitmol)
        self.nxmol     = rdkit_to_nx(self.rdkitmolH)

        if len(bb) == 0:
            bb = [pwd + '/buildingblock.xml']
        bb_final = bb + bb_add
        self.graphs = ReadBuildingBlock(bb_final)

        # Perfect Match
        self.PerfectMatch()

        # Check
        self.HasType()

        # StereoChemistry
        self.cis, self.trans = CisTrans(self.rdkitmol)
        self.chiral = Chiral(self.rdkitmolH)


    def PerfectMatch(self):
        for ffgraph in self.graphs:
            nodes    = ffgraph.nodes(data=True)
            GM       = nx.isomorphism.GraphMatcher(self.nxmol, ffgraph, node_match=node_match)
            mappings = [mapping for mapping in GM.subgraph_isomorphisms_iter()]

            for mapping in mappings:
                for index, name in mapping.items():
                    self.nxmol.nodes[index]['type']   = nodes[name]['type']
                    self.nxmol.nodes[index]['charge'] = nodes[name]['charge']
                    self.nxmol.nodes[index]['z']      = nodes[name]['z']

    def HasType(self):
        qtot = 0.0
        for nodeid, attr in self.nxmol.nodes(data=True):
            if 'type' not in attr.keys():
                print(f"{nodeid} has no assigned atomtype")
            if 'charge' in attr.keys():
                qtot += attr['charge']

        if abs(int(qtot) - qtot) > 0.0001:
            print(f'qtot = {qtot}')



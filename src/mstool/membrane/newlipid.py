from rdkit import Chem
import networkx as nx
import os
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib.pyplot as plt
from   ..core.readxml import ReadXML
pwd = os.path.dirname(os.path.realpath(__file__))

E2Z = {
    "H": 1,   "He": 2,  "Li": 3,  "Be": 4,  "B": 5,
    "C": 6,   "N": 7,   "O": 8,   "F": 9,   "Ne": 10,
    "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15,
    "S": 16,  "Cl": 17, "Ar": 18, "K": 19,  "Ca": 20,
    "Sc": 21, "Ti": 22, "V": 23,  "Cr": 24, "Mn": 25,
    "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30,
    "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35,
    "Kr": 36, "Rb": 37, "Sr": 38, "Y": 39,  "Zr": 40,
    "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44, "Rh": 45,
    "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
    "Sb": 51, "Te": 52, "I": 53,  "Xe": 54, "Cs": 55,
    "Ba": 56, "La": 57, "Ce": 58, "Pr": 59, "Nd": 60,
    "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65,
    "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70,
    "Lu": 71, "Hf": 72, "Ta": 73, "W": 74,  "Re": 75,
    "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80,
    "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85,
    "Rn": 86, "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90,
    "Pa": 91, "U": 92,  "Np": 93, "Pu": 94, "Am": 95,
    "Cm": 96, "Bk": 97, "Cf": 98, "Es": 99, "Fm": 100,
    "Md": 101,"No": 102,"Lr": 103,"Rf": 104,"Db": 105,
    "Sg": 106,"Bh": 107,"Hs": 108,"Mt": 109,"Ds": 110,
    "Rg": 111,"Cn": 112,"Nh": 113,"Fl": 114,"Mc": 115,
    "Lv": 116,"Ts": 117,"Og": 118,
}

Z2E = {value: key for key, value in E2Z.items()}

def VisG(G, attrs=True):
    # Layout positions
    pos = nx.spring_layout(G)

    # Custom labels
    if attrs:
        labels = {}
        for n in G.nodes:
            tt = G.nodes[n]['type']
            qq = G.nodes[n]['charge']
            labels[n] = f"{n} ({tt}, q={qq})"

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
    # Kekulize the molecule in-place
    # aromatic bonds into single or double
    Chem.Kekulize(mol, clearAromaticFlags=True)

    bond_order_map = {
        Chem.rdchem.BondType.SINGLE: 1,
        Chem.rdchem.BondType.DOUBLE: 2,
        Chem.rdchem.BondType.TRIPLE: 3,
        Chem.rdchem.BondType.AROMATIC: 1.5,
        Chem.rdchem.BondType.QUADRUPLE: 4,
    }

    G = nx.Graph()

    for atom in mol.GetAtoms():
        G.add_node(atom.GetIdx(), atomic_number=atom.GetAtomicNum(), z=0.0)

    for bond in mol.GetBonds():
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), order=bond_order_map[bond.GetBondType()])

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

    if isinstance(ifiles, str):
        ifiles = [ifiles]

    for ifile in ifiles:
        root = ET.parse(ifile).getroot()

        for block in root.findall('BuildingBlock'):
            name = block.get("name") or "Unnamed"
            G = nx.Graph()
            G.name = name

            z = 0.0
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
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))

    for center_idx, config in chiral_centers:
        center_atom = mol.GetAtomWithIdx(center_idx)
        neighbors = center_atom.GetNeighbors()

        if len(neighbors) != 4:
            print(f"Warning: atom {center_idx} is not bonded to 4 atoms.")
            continue

        neighbors = center_atom.GetNeighbors()
        # print(f"Chiral center at atom {center_idx}, config: {config}")

        # Collect neighbor ranks and symbols
        neighbor_info = []
        for nbr in neighbors:
            idx = nbr.GetIdx()
            rank = ranks[idx]
            neighbor_info.append((rank, idx, nbr.GetSymbol()))
        neighbor_info.sort()

        if config == 'R':
            chirals.append([neighbor_info[3][1],
                            center_idx,
                            neighbor_info[0][1],
                            neighbor_info[1][1],
                            neighbor_info[2][1]])
        elif config == 'S':
            chirals.append([neighbor_info[3][1],
                            center_idx,
                            neighbor_info[0][1],
                            neighbor_info[2][1],
                            neighbor_info[1][1]])

    return chirals


class NewLipid:
    """Create a new lipid type based on SMILES.
    >>> from mstool.membrane.newlipid import NewLipid
    >>> POPC    = NewLipid('POPC1',   'CCCCCCCC/C=C\CCCCCCCC(O[C@H](COC(CCCCCCCCCCCCCCC)=O)COP(OCC[N+](C)(C)C)([O-])=O)=O')
    >>> POPE    = NewLipid('POPE1',   '[H][C@@](COP([O-])(OCC[NH3+])=O)(OC(CCCCCCC/C=C\CCCCCCCC)=O)COC(CCCCCCCCCCCCCCC)=O')
    >>> BMPSS   = NewLipid('BMPSS',   '[H][C@](COP(OC[C@@]([H])(O)COC(CCCCCCC/C=C\CCCCCCCC)=O)([O-])=O)(O)COC(CCCCCCC/C=C\CCCCCCCC)=O')
    >>> BMPSR   = NewLipid('BMPSR',   '[H][C@@](COP(OC[C@@]([H])(O)COC(CCCCCCC/C=C\CCCCCCCC)=O)([O-])=O)(O)COC(CCCCCCC/C=C\CCCCCCCC)=O')
    >>> BMPRS   = NewLipid('BMPRS',   '[H][C@](COP(OC[C@]([H])(O)COC(CCCCCCC/C=C\CCCCCCCC)=O)([O-])=O)(O)COC(CCCCCCC/C=C\CCCCCCCC)=O')
    >>> BMPRR   = NewLipid('BMPRR',   '[H][C@@](COP(OC[C@]([H])(O)COC(CCCCCCC/C=C\CCCCCCCC)=O)([O-])=O)(O)COC(CCCCCCC/C=C\CCCCCCCC)=O')
    >>> SM102N  = NewLipid('SM102N',  'O=C(OC(CCCCCCCC)CCCCCCCC)CCCCCCC[N](CCO)CCCCCC(OCCCCCCCCCCC)=O')
    >>> SM102P  = NewLipid('SM102P',  'O=C(OC(CCCCCCCC)CCCCCCCC)CCCCCCC[NH+](CCO)CCCCCC(OCCCCCCCCCCC)=O')
    >>> PEG2000 = NewLipid('PEG2000', 'O=C(CCCCCCCCCCCCC)OCC(OC(CCCCCCCCCCCCC)=O)COCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOC')
    >>> VisG(POPC.nxmol)
    >>> POPC.WriteMapping('lipid.itp')
    >>> POPC.WriteFF('ff.xml')
    >>> POPE.WriteMapping('lipid.itp')
    >>> POPE.WriteFF('ff.xml')
    """
    def __init__(self, resname, smiles, bb=[], bb_add=[], ff=[], ff_add=[]):
        self.resname   = resname
        self.smiles    = smiles
        self.rdkitmol  = Chem.MolFromSmiles(smiles)
        # added hydrogens will be at the end of the molecule
        # therefore, it will not mess up with the original atom order
        self.rdkitmolH = Chem.AddHs(self.rdkitmol)
        self.nxmol     = rdkit_to_nx(self.rdkitmolH)
        self.xml       = ReadXML(ff=ff, ff_add=ff_add)

        if isinstance(bb, str):
            bb = [bb]
        if isinstance(bb_add, str):
            bb_add = [bb_add]
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
        Chem.AssignStereochemistry(self.rdkitmolH, force=True, cleanIt=True)
        self.chiral = Chiral(self.rdkitmolH)

        # names
        self.names = []
        for nodeid, attr in self.nxmol.nodes(data=True):
            name = Z2E[attr['atomic_number']] + str(nodeid)
            self.names.append(name)
            self.nxmol.nodes[nodeid]['name'] = name

        # Amide
        self.AmideMatch()


    def WriteMapping(self, ofile):
        """Write the mapping of the new lipid to a file."""
        data = {'NC3': [], 'PO4': [], 'GL1': [], 'C3A': []}
        for nodeid, attr in self.nxmol.nodes(data=True):
            posz = attr.get('z', 0.0)
            if posz <= -5:
                data['C3A'].append(attr['name'])
            elif posz > -5 and posz <= 0:
                data['GL1'].append(attr['name'])
            elif posz > 0 and posz <= 2.5:
                data['PO4'].append(attr['name'])
            elif posz > 2.5:
                data['NC3'].append(attr['name'])

        with open(ofile, 'a') as f:
            f.write(f'\nRESI {self.resname}\n')

            for key, atoms in data.items():
                if atoms:
                    f.write(f'[ {key} ]\n')
                    for atom in atoms:
                        f.write(f'{atom}\n')
                    f.write('\n')

            for stereo, values in [('chiral', self.chiral), ('cis', self.cis), ('trans', self.trans)]:
                if len(values) > 0:
                    f.write(f'[ {stereo} ]\n')
                    for value in values:
                        text = ''
                        for v in value:
                            text += self.names[v] + ' '
                        f.write(text.strip() + '\n')
                    f.write('\n')

    def WriteFF(self, ofile):
        """Append this lipid's force field info as a <Residue> to the XML file"""
        # Read existing content except closing tags
        realcontent = []; bondcontent = []
        if os.path.exists(ofile):
            with open(ofile, 'r') as fin:
                for line in fin:
                    if line.strip().startswith('<ForceField>'):
                        continue
                    elif line.strip().startswith('<Residues>'):
                        continue
                    elif line.strip().startswith('</ForceField>'):
                        continue
                    elif line.strip().startswith('</Residues>'):
                        continue
                    elif line.strip().startswith('<HarmonicBondForce>'):
                        continue
                    elif line.strip().startswith('</HarmonicBondForce>'):
                        continue
                    elif line.strip().startswith('<Bond k='):
                        bondcontent.append(line)
                        continue
                    realcontent.append(line)

        # Build residue block as lines
        harmonic_bonds = dict()
        lines = [f'    <Residue name="{self.resname}">']
        for nodeid, attr in self.nxmol.nodes(data=True):
            lines.append(
                f'      <Atom charge="{attr.get("charge", 0.0)}" name="{attr.get("name", f"A{nodeid}")}" type="{attr.get("type", "UNK")}"/>'
            )
        for atom1, atom2, bond_attr in self.nxmol.edges(data=True):
            name1 = self.nxmol.nodes[atom1]["name"]
            name2 = self.nxmol.nodes[atom2]["name"]
            lines.append(f'      <Bond atomName1="{name1}" atomName2="{name2}"/>')
            type1 = self.nxmol.nodes[atom1]["type"]
            type2 = self.nxmol.nodes[atom2]["type"]
            if type1 > type2: type1, type2 = type2, type1
            if (type1, type2) not in self.xml.BOND.keys():
                numb1 = self.nxmol.nodes[atom1]['atomic_number']
                numb2 = self.nxmol.nodes[atom2]['atomic_number']
                if numb1 == 1 or numb2 == 1:
                    harmonic_bonds[(type1, type2)] = [0.1111, 269449.6]
                else:
                    harmonic_bonds[(type1, type2)] = [0.1530, 186188.0]
        lines.append('    </Residue>')

        # Write everything back with blank lines between each line
        with open(ofile, 'w') as f:
            f.write('<ForceField>\n  <Residues>\n')
            if realcontent:
                for line in realcontent:
                    f.write(line)

            for line in lines:
                f.write(line + '\n')
            f.write('  </Residues>\n')

            if bondcontent or harmonic_bonds:
                f.write('  <HarmonicBondForce>\n')
                for bond in bondcontent:
                    f.write(bond)
                for key, value in harmonic_bonds.items():
                    f.write(f'    <Bond k="{value[1]}" length="{value[0]}" type1="{key[0]}" type2="{key[1]}"/>\n')
                f.write('  </HarmonicBondForce>\n')

            f.write('</ForceField>\n')

    def WriteMartini(self, ofile):
        with open(ofile, 'a') as f:
            f.write(f"""
 [ moleculetype ]
 ; molname      nrexcl
   {self.resname}          1

 [ atoms ]
 ; id    type    resnr   residu  atom    cgnr    charge
    1    Q0   1  {self.resname}    NC3      1  1.0
    2    Qa   1  {self.resname}    PO4      2  -1.0
    3    Na   1  {self.resname}    GL1      3  0
    4    Na   1  {self.resname}    GL2      4  0
    5    C1   1  {self.resname}    C1A      5  0
    6    C3   1  {self.resname}    D2A      6  0
    7    C1   1  {self.resname}    C3A      7  0
    8    C1   1  {self.resname}    C4A      8  0
    9    C1   1  {self.resname}    C1B      9  0
   10    C1   1  {self.resname}    C2B     10  0
   11    C1   1  {self.resname}    C3B     11  0
   12    C1   1  {self.resname}    C4B     12  0

 [ bonds ]
 ;  i  j     funct   length  force.c.
    1  2     1   0.47    1250
    2  3     1   0.47    1250
    3  4     1   0.37    1250
    3  5     1   0.47    1250
    5  6     1   0.47    1250
    6  7     1   0.47    1250
    7  8     1   0.47    1250
    4  9     1   0.47    1250
    9 10     1   0.47    1250
   10 11     1   0.47    1250
   11 12     1   0.47    1250

 [ angles ]
 ;  i  j  k  funct   angle   force.c.
    2  3  4  2   120.0   25.0
    2  3  5  2   180.0   25.0
    3  5  6  2   180.0   25.0
    5  6  7  2   120.0   45.0
    6  7  8  2   180.0   25.0
    4  9 10  2   180.0   25.0
    9 10 11  2   180.0   25.0
   10 11 12  2   180.0   25.0
\n""")

    def PerfectMatch(self):
        for ffgraph in self.graphs:
            nodes    = ffgraph.nodes(data=True)
            GM       = nx.isomorphism.GraphMatcher(self.nxmol, ffgraph, node_match=node_match, edge_match=edge_match)
            mappings = [mapping for mapping in GM.subgraph_isomorphisms_iter()]

            for mapping in mappings:
                #print(ffgraph.name)
                for index, name in mapping.items():
                    self.nxmol.nodes[index]['type']   = nodes[name]['type']
                    self.nxmol.nodes[index]['charge'] = nodes[name]['charge']
                    self.nxmol.nodes[index]['z']      = nodes[name]['z']

    def AmideMatch(self):
        amide = nx.Graph()
        amide.add_node('H', atomic_number=1, type='H')
        amide.add_node('N', atomic_number=7, type='N')
        amide.add_node('C', atomic_number=6, type='C')
        amide.add_node('O', atomic_number=8, type='O')
        amide.add_edge('H', 'N', order=1)
        amide.add_edge('N', 'C', order=1)
        amide.add_edge('C', 'O', order=2)

        GM = nx.isomorphism.GraphMatcher(self.nxmol, amide, node_match=node_match, edge_match=edge_match)
        for mapping in GM.subgraph_isomorphisms_iter():
            if ''.join(mapping.values()) == 'HNCO':
                self.trans.append(list(mapping.keys()))


    def HasType(self):
        qtot = 0.0
        for nodeid, attr in self.nxmol.nodes(data=True):
            if 'type' not in attr.keys():
                print(f"{nodeid} has no assigned atomtype")
            if 'charge' in attr.keys():
                qtot += attr['charge']

        self.charge = qtot
        if abs(int(qtot) - qtot) > 0.0001:
            print(f'qtot = {qtot}')



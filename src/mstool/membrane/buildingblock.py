from ..core.checkff import mol2nx
import networkx as nx
import numpy as np

def node_match(n1, n2):
    return n1['atomic_number'] == n2['atomic_number']

def FS(G, i=0, j=-1):
    nodes   = [node for node in G.nodes()]
    visited = set()
    queue   = [nodes[i]]
    subGs   = []
    q = 0.0

    while queue:
        current = queue.pop(j)
        q += G.nodes[current]['charge']
        if current not in visited:
            visited.add(current)
            for neighbor in G.neighbors(current):
                if neighbor not in visited:
                    queue.append(neighbor)

        if abs(q - int(q)) < 1e-6:
            subG = G.subgraph(visited).copy()
            try:
                subG.resname = G.resname
            except:
                pass
            subGs.append(subG)

    return subGs


def BuildingBlock(mol, resname=None):
    """ Get building blocks from a molecule.
    >>> subGs = BuildingBlock(mstool.ReadXML().RESI['POPC'], 'POPC')
    """
    molnx = mol2nx(mol, resname=resname)
    subGs = []
    for i in range(len(molnx.nodes)):
        subGs.extend(FS(molnx, i, j=-1)) #DFS
        subGs.extend(FS(molnx, i, j=0))  #BFS

    hashes = set()
    final = []
    for subG in subGs:
        if subG.number_of_nodes() == len(molnx.nodes):
            continue
        h = nx.weisfeiler_lehman_graph_hash(subG)
        if h not in hashes:
            hashes.add(h)
            final.append(subG)
    return final

def nx2xml(nxmols, ofile, lower_cutoff=3.5, upper_cutoff=49.5):
    """ Create an XML file for building blocks from a list of networkx molecules.
    >>> d = mstool.ReadToppars('/Users/siyoungkim/mstool/src/mstool/FF/charmm36.toppar/toppar.lipid.str')
    >>> mstool.nx2xml(d.GRAPH, 'bb.xml')
    """
    lines = []
    for i, nxmol in enumerate(nxmols):
        if len(nxmol.nodes) < lower_cutoff or len(nxmol.nodes) > upper_cutoff:
            continue

        lines.append(f'  <BuildingBlock name="{nxmol.resname}_{i}">')

        charges = np.array([nxmol.nodes[n]['charge'] for n in nxmol.nodes])
        if max(abs(charges)) > 0.55:
            lines.append(f'    <Pos z="2.5"/>')
        elif max(abs(charges)) < 0.28:
            lines.append(f'    <Pos z="-5.0"/>')
        else:
            lines.append(f'    <Pos z="0.0"/>')

        for nodeid, attr in nxmol.nodes(data=True):
            lines.append(
                f'    <Atom name="{nodeid}" number="{attr.get("atomic_number", "0")}" type="{attr.get("type", "UNK")}" charge="{attr.get("charge", 0.0)}"/>'
            )
        for atom1, atom2, bond_attr in nxmol.edges(data=True):
            if 'order' in bond_attr.keys():
                border = bond_attr['order']
                if border != 1:
                    lines.append(f'    <Bond atomName1="{atom1}" atomName2="{atom2}" order="{border}"/>')
                else:
                    lines.append(f'    <Bond atomName1="{atom1}" atomName2="{atom2}"/>')
            else:
                lines.append(f'    <Bond atomName1="{atom1}" atomName2="{atom2}"/>')
        lines.append('  </BuildingBlock>')

    # Write everything back with blank lines between each line
    with open(ofile, 'w') as f:
        f.write('<BuildingBlocks>\n')
        for line in lines:
            f.write(line + '\n')
        f.write('</BuildingBlocks>')

def CreateBuildingBlockXML(readxml, ofile, select=[], lower_cutoff=3.5, upper_cutoff=49.5):
    """ Read XML, make groups based on XML, and create a building block XML file.
    >>> readxml = mstool.ReadXML()
    >>> mstool.CreateBuildingBlockXML(readxml, 'bb.xml')
    """

    subGs = []
    for resname in readxml.RESI.keys():
        if select and resname not in select:
            continue
        if len(readxml.RESI[resname]['names']) < lower_cutoff:
            continue
        subGs.extend(BuildingBlock(readxml.RESI[resname], resname=resname))

    hashes = set()
    final = []
    for subG in subGs:
        h = nx.weisfeiler_lehman_graph_hash(subG)
        if h not in hashes:
            hashes.add(h)
            final.append(subG)
    final = sorted(final, key=lambda g: len(g.nodes)) 
    print("Number of unique building blocks found:", len(final))
    nx2xml(final, ofile, lower_cutoff, upper_cutoff)


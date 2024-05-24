import networkx as nx
import matplotlib.pyplot as plt

def visMol(mol, out=None):
    '''visualize a molecule.
    CG molecules:
    >>> from mstool.utils.vismol import visMol
    >>> martini = mstool.ReadMartini()
    >>> bonds = martini.martini['molecules']['POPC']
    >>> visMol(bonds)

    AA molecules:
    >>> ff = mstool.ReadXML()
    >>> bonds = ff.RESI['POPC']
    >>> visMol(bonds)
    '''

    G = nx.Graph()
    try:
        # Martini
        names = mol['atoms']['name']
    except:
        # XML
        names = mol['names']

    for name in names:
        G.add_node(name)
    
    bonds = mol['bonds']
    for bond in bonds:
        if G.has_edge(bond[0], bond[1]):
            print(f"bond between {bond[0]} and {bond[1]} is defined more than once")
        G.add_edge(bond[0], bond[1])
    
    # pos = nx.spring_layout(G)
    # pos = nx.spectral_layout(G)
    # pos = nx.circular_layout(G)
    nx.draw(G, with_labels=True)
    if out:
        plt.savefig(out)
    else:
        plt.show()
        

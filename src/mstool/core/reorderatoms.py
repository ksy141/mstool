import pandas as pd
from .universe import Universe
from .readxml  import ReadXML
#pwd = os.path.dirname(os.path.realpath(__file__))

def reorder_atoms(group, xml):
    resname = group['resname'].values[0]
    names   = group['name'].values
    ffnames = list(xml.RESI[resname]['names'])

    #if set(names) != set(ffnames):
    #    for resname, value in xml.RESI.items():
    #        if set(names) == set(value['names']):
    #            ffnames = xml.RESI[resname]['names']
    #            break
    
    H_missing_atoms = []; O_missing_atoms = []
    for name in group['name']:
        if name not in ffnames:
            if name.startswith('H'):
                H_missing_atoms.append(name)
            else:
                O_missing_atoms.append(name)
    
    # the first one is likely N
    all_names = ffnames[:1] + H_missing_atoms + ffnames[1:] + O_missing_atoms
    group['name'] = pd.Categorical(group['name'], categories=all_names, ordered=True)
    return group.sort_values('name')

def ReorderAtoms(structure, out=None, ff=[], ff_add=[]):
    if isinstance(structure, Universe):
        u = structure
    else:
        u = Universe(structure)
    
    if not isinstance(ff_add, list): ff_add = [ff_add]
    xml = ReadXML(ff=ff, ff_add=ff_add)

    u.atoms = u.atoms.groupby('resn').apply(lambda group: reorder_atoms(group, xml))
    u.atoms.reset_index(drop=True, inplace=True)
    
    if out: u.write(out)
    return u


import numpy as np
from   scipy.spatial.distance import pdist
from   ..core.universe import Universe

class Rock:
    def __init__(self, structure, out='ROCK', resname='ROCK', type='HL', rcut=1.2, Kbond=5000.0, name='HROCK'):
        
        # Kbond: kJ/mol/nm^2
        u    = Universe(structure)
        N    = len(u.atoms)
        pos  = u.atoms[['x','y','z']].to_numpy()
        dist = pdist(pos)


        ### Change atomic information
        u.atoms.resname = resname
        u.atoms.chain   = resname
        u.atoms.type    = type
        u.atoms.resid   = 1
        u.atoms.name    = [name + str(i) for i in range(len(u.atoms))]
        #u.atoms.mass    = 12.0
        #u.atoms.anum    = 6
        u.atoms.mass     = 1.00
        u.atoms.anum     = 1
        u.atoms.bfactor = 0.0
        u.atoms.charge  = 0.0
        self.dms = out + '.dms'
        u.write(self.dms)


        ### Save it to XML
        self.xml = out + '.xml'
        W = open(self.xml, 'w')
        W.write('<ForceField>\n')
        W.write('  <Residues>\n')
        W.write(f'    <Residue name="{resname}">\n')
        for name in u.atoms.name.tolist():
            W.write(f'      <Atom charge="0.00" name="{name}" type="{type}"/>\n')
        W.write(f'    </Residue>')
        W.write('  </Residues>\n')
        W.write('</ForceField>\n')
        W.close()


        ### Bonds
        bonds = []
        for i in range(N):
            for j in range(N):
                if i < j:
                    k = N * i + j - ((i + 2) * (i + 1)) // 2
                    d = dist[k] * 0.1 #A to nm
                    
                    if d < rcut:
                        bonds.append([i, j, d, Kbond])
                        
        self.bonds = np.array(bonds)


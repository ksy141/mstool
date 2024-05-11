import numpy as np
from   scipy.spatial.distance import pdist
from   ..core.universe import Universe

class Rock:
    def __init__(self, structure, out='ROCK', resname='ROCK',
                 rockHtype='HAL3', rockCtype='CTL3', 
                 rcut=1.2, Kbond=5000.0, ENM=False):
        
        # Kbond: kJ/mol/nm^2
        u    = Universe(structure)
        N    = len(u.atoms)
        pos  = u.atoms[['x','y','z']].to_numpy()

        ### Change atomic information
        u.atoms.resname = resname
        u.atoms.chain   = resname
        u.atoms.resid   = 1
        bAH = u.atoms['name'].str.startswith('H')
        u.atoms.loc[bAH,  'type'] = rockHtype
        u.atoms.loc[~bAH, 'type'] = rockCtype
        u.atoms.name = [u.atoms.iloc[i]['type'][0] + 'ROCK' + str(i) for i in range(len(u.atoms))]
        u.atoms.mass = 100.0

        #if rocktype[0] == 'H':
        #    u.atoms.mass     = 1.00
        #    u.atoms.anum     = 1
        #elif rocktype[0] == 'C':
        #    u.atoms.mass    = 12.0
        #    u.atoms.anum    = 6
        
        if ENM:
            u.atoms.bfactor = 0.0
        else:
            u.atoms.bfactor = 1.0
        u.atoms.charge  = 0.0
        self.dms = out + '.dms'
        u.write(self.dms)


        ### Save it to XML
        self.xml = out + '.xml'
        W = open(self.xml, 'w')
        W.write('<ForceField>\n')
        W.write('  <Residues>\n')
        W.write(f'    <Residue name="{resname}">\n')
        for nn, tt in zip(u.atoms.name.tolist(), u.atoms.type.tolist()):
            W.write(f'      <Atom charge="0.00" name="{nn}" type="{tt}"/>\n')
        W.write(f'    </Residue>\n')
        W.write('  </Residues>\n')
        W.write('</ForceField>\n')
        W.close()


        ### Bonds
        self.bonds = []

        if ENM:
            dist  = pdist(pos)
            bonds = []
            for i in range(N):
                for j in range(N):
                    if i < j:
                        k = N * i + j - ((i + 2) * (i + 1)) // 2
                        d = dist[k] * 0.1 #A to nm
                        
                        if d < rcut:
                            bonds.append([i, j, d, Kbond])
                            
            self.bonds = np.array(bonds)


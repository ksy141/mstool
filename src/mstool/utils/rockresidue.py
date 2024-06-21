from   ..core.universe import Universe

def generate_list(group):
    index = group.groupby('resn').cumcount()
    return [f"{'H' if group['name'].iloc[i][0] == 'H' else 'C'}{index.iloc[i]}" for i in range(len(group))]

def get_resname(group, uniques):
    resname = 'ROCK' + ''.join(group['name'])
    return [resname] * len(group['name'])


class RockResidue:
    def __init__(self, structure, out='ROCK', resname='ROCK',
                 rockHtype='HAL3', rockCtype='CTL3'):
        
        self.dms = out + '.dms'
        u = Universe(structure)
        u.bonds         = []
        u.atoms.mass    = 100.0
        u.atoms.bfactor = 1.0
        u.atoms.charge  = 0.0

        # change name to C0 H1 H2 C3 H4 H5 for each residue
        u.atoms['name'] = u.atoms.groupby('resn').apply(generate_list).explode().reset_index(level=0, drop=True)
        anum, mass, vdw = u.update_anum_mass_vdw_from_name()
        u.atoms['anum'] = anum
        
        # obtain unique residues
        uniques = {}
        for name in u.atoms.groupby('resn').name:
            namelist = name[1].tolist()
            resname  = 'ROCK' + ''.join(namelist)

            #if len(namelist) != len(set(namelist)):
            #    assert 0 == 1, print(namelist)

            if resname not in uniques.keys():
                uniques[resname] = namelist


        # change resname
        u.atoms['resname'] = u.atoms.groupby('resn').apply(get_resname, uniques=uniques).explode().reset_index(level=0, drop=True)
        u.write(self.dms)
        self.bonds = []
        

        ### Save it to XML
        self.xml = out + '.xml'
        W = open(self.xml, 'w')
        W.write('<ForceField>\n')
        W.write('  <Residues>\n')
    
        for unique_resname, unique_name in uniques.items():
            W.write(f'    <Residue name="{unique_resname}">\n')
            for nn in unique_name:
                if nn.startswith('C'):
                    tt = rockCtype
                else:
                    tt = rockHtype
                W.write(f'      <Atom charge="0.00" name="{nn}" type="{tt}"/>\n')
            W.write(f'    </Residue>\n')

        W.write('  </Residues>\n')
        W.write('</ForceField>\n')
        W.close()


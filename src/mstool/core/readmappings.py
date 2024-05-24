import os
pwd = os.path.dirname(os.path.realpath(__file__))
scheme = [pwd + '/../mapping/martini.protein.c36m.dat',
          pwd + '/../mapping/martini.lipid.c36.dat']


class ReadMappings:
    def __init__(self, mapping = [], mapping_add = []):

        if not isinstance(mapping, list):
            mapping = [mapping]

        if not isinstance(mapping_add, list): 
            mapping_add = [mapping_add]

        if mapping == []:
            mapping = scheme

        self.ifiles = [*mapping, *mapping_add]
        self.collect()


    def get_groups(self, seq, group_by):
        data = []
        for line in seq:
            line = line.strip()
            if not line: continue
            line = line.split(';')[0]
            if not line: continue
            if line.startswith(group_by):
                if data:
                    yield data, index
                    data = []

            data.append(line)

        # return the last one
        if data:
            yield data, index


    def collect(self):
        d = {}
        for ifile in self.ifiles:
            with open(ifile) as W:
                for line in W.readlines():
                    if line.startswith(';'): continue
                    if not line.strip(): continue

                    if line.startswith('RESI'):
                        sl = line.split()
                        resname = sl[1]

                        if resname in d.keys():
                            print('%s duplicate mapping scheme' %resname)
                            print('updating with the latest one\n')

                        d[resname] = {}
                        d[resname]['CGAtoms']      = {}
                        d[resname]['AAAtoms']      = []
                        d[resname]['cis']          = []
                        d[resname]['trans']        = []
                        d[resname]['chiral']       = []
                        d[resname]['dihedral']     = []
                        d[resname]['antidihedral'] = []
                        continue

                    #if line.startswith(('[ cis', '[cis')):
                    #    read = 'cis'
                    #    continue

                    #if line.startswith(('[ trans', '[trans')):
                    #    read = 'trans'
                    #    continue

                    #if line.startswith(('[ chiral', '[chiral')):
                    #    read = 'chiral'
                    #    continue

                    if line.startswith('['):
                        # isomeric information
                        read = line.split('[')[1].split(']')[0].strip().lower()
                        if read in ['cis', 'trans', 'chiral', 'chirals', 'dihedral', 'dihedrals', 'antidihedral', 'antidihedrals']:
                            continue
                        
                        # CG bead
                        CGName = line.split('[')[1].split(']')[0].strip()
                        d[resname]['CGAtoms'][CGName] = []
                        read = 'atoms'
                        continue

                    if read == 'atoms':
                        d[resname]['CGAtoms'][CGName] += line.split()
                        d[resname]['AAAtoms']         += line.split()
                        continue

                    if read == 'cis':
                        d[resname]['cis'].append(line.split())
                        continue

                    if read == 'trans':
                        d[resname]['trans'].append(line.split())
                        continue

                    if read == 'chiral' or read == 'chirals':
                        d[resname]['chiral'].append(line.split())
                        continue

                    if read == 'dihedral' or read == 'dihedrals':
                        sl = line.split()
                        if len(sl) % 5 != 0:
                            assert 0 == 1, \
                            'dihedral should be atomNameA atomNameB atomNameC atomNameD dihedral(degree)'

                        for i in range(4, len(sl), 5):
                            sl[i] = float(sl[i])

                        d[resname]['dihedral'].append(sl)
                        continue

                    if read == 'antidihedral' or read == 'antidihedrals':
                        sl = line.split()
                        if len(sl) % 5 != 0:
                            assert 0 == 1, \
                            'antidihedral should be atomNameA atomNameB atomNameC atomNameD dihedral(degree)'

                        for i in range(4, len(sl), 5):
                            sl[i] = float(sl[i])

                        d[resname]['antidihedral'].append(sl)
                        continue

        self.RESI = d


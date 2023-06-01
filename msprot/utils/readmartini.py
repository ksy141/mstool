import os
pwd = os.path.dirname(os.path.realpath(__file__))

class ReadMartini:
    def __init__(self,
        ff = [pwd + '/../../FF/martini2.2/martini_v2.2.itp',
              pwd + '/../../FF/martini2.2/martini_v2.0_ions.itp',
              pwd + '/../../FF/martini2.2/POPC.itp',
              pwd + '/../../FF/martini2.2/martini_v2.2_proteins/proteins.itp'],
              #pwd + '/../../FF/martini2.2/martini_v2.0_lipids_all_201506.itp'],
        ff_add = [],
        Kc2b   = 10000.0):

        self.Kc2b = Kc2b
        if not isinstance(ff_add, list): ff_add = [ff_add]
        self.ifiles = [*ff, *ff_add]
        self.collect()

    def collect(self):
        d = {}
        for ifile in self.ifiles:
            with open(ifile) as W:
                for line in W.readlines():
                    if line.startswith(';'): continue
                    if not line.strip(): continue
                    line = line.split(';')[0].strip()

                    if '[' in line and ']' in line:
                        read = line.split('[')[1].split(']')[0].strip()
                        continue

                    if read == 'defaults':
                        ### non-bonded function type
                        #   1 for LJ
                        #   2 for Buckingham

                        ### combination rule = 1
                        #  Vii = C6  = 4 eps sig**6 [kJ/mol/nm^6]
                        #  Wii = C12 = 4 eps sig**12 [kJ/mol/nm^12]
                        #  Vij = (C6i C6j)^(0.5)
                        #  Wij = (C12i C12j)^(0.5)

                        ### combination rule = 2
                        #  Vii = sigma [nm]
                        #  Wii = eps   [kJ/mol]
                        #  Vij = 0.5 * (sig1 + sig2)
                        #  Wij = (eps1 * eps2)^(0.5)

                        # gen-pairs - default is no
                        # fudge LJ  - default is 1 
                        # fudge QQ  - default is 1

                        sl = line.split()
                        if sl[0] == 1:
                            d['nbfunc'] = 'LJ'
                        if sl[0] == 2:
                            d['nbfunc'] = 'Buckingham'

                        if sl[1] == 1:
                            d['energy'] = 'C12/r^12 - C6/r^6'
                        if sl[1] == 2:
                            d['energy'] = '4*eps*((sigma/r)^12 - (sigma/r)^6)'

                    if read == 'atomtypes':
                        if 'atomtypes' not in d.keys(): d['atomtypes'] = {}
                        
                        sl = line.split()
                        t = sl[0]
                        m = float(sl[1])
                        q = float(sl[2])
                        V = float(sl[4])
                        W = float(sl[5])

                        d['atomtypes'][t] = {}
                        d['atomtypes'][t]['m'] = float(sl[1])
                        d['atomtypes'][t]['q'] = float(sl[2])
                        d['atomtypes'][t]['V'] = float(sl[4])
                        d['atomtypes'][t]['W'] = float(sl[5])

                    if read == 'nonbond_params':
                        if 'nb' not in d.keys(): d['nb'] = {}
                        sl = line.split()
                        t1 = sl[0]; t2 = sl[1]
                        V  = float(sl[3])
                        W  = float(sl[4])
                        d['nb'][(t1,t2)] = [V, W]

                    if read == 'moleculetype':
                        if 'molecules' not in d.keys(): d['molecules'] = {}

                        sl = line.split()
                        resname = sl[0]
                        d['molecules'][resname] = {}
                        d['molecules'][resname]['nrexcl'] = int(sl[1])
                        d['molecules'][resname]['atoms']  = {'id':[], 'type':[], 
                                                             'name':[], 'q':[]}
                        d['molecules'][resname]['bonds']          = []
                        d['molecules'][resname]['angles']         = []
                        d['molecules'][resname]['dihedrals']      = []
                        d['molecules'][resname]['exclusions']     = []
                        d['molecules'][resname]['virtual_sites3'] = []

                    if read == 'atoms':
                        sl = line.split()
                        d['molecules'][resname]['atoms']['id'].append(int(sl[0]))
                        d['molecules'][resname]['atoms']['type'].append(sl[1])
                        d['molecules'][resname]['atoms']['name'].append(sl[4])
                        d['molecules'][resname]['atoms']['q'].append(float(sl[6]))

                    if read == 'bonds':
                        # V = 0.5 * k * (r - r0)^2
                        sl = line.split()
                        id1 = int(sl[0]); id2 = int(sl[1])
                        name1 = d['molecules'][resname]['atoms']['name'][id1-1]
                        name2 = d['molecules'][resname]['atoms']['name'][id2-1]
                        #if name1 > name2: name2 = name1
                        l = float(sl[3]) #nm
                        k = float(sl[4]) #kJ/nm^2
                        d['molecules'][resname]['bonds'].append([name1, name2, l, k])

                    if read == 'angles':
                        # V = 0.5 * k * (t - t0)^2 if func = 1
                        # V = 0.5 * k * (cost - cost0)^2 if func = 2
                        sl = line.split()
                        id1 = int(sl[0]); id2 = int(sl[1]); id3 = int(sl[2])
                        name1 = d['molecules'][resname]['atoms']['name'][id1-1]
                        name2 = d['molecules'][resname]['atoms']['name'][id2-1]
                        name3 = d['molecules'][resname]['atoms']['name'][id3-1]
                        l = float(sl[4]) #degree
                        k = float(sl[5]) #kJ/mol
                        func = int(sl[3])
                        d['molecules'][resname]['angles'].append([name1, name2, name3, l, k, func])

                    if read == 'dihedrals':
                        # proper (1):   V = k * (1 + cos(nt - t0)); t, k, n
                        # improper (2): V = 0.5 * k * (t - t0)^2;   t, k (kJ/mol/rad^2)
                        sl = line.split()
                        id1 = int(sl[0]); id2 = int(sl[1]); id3 = int(sl[2]); id4 = int(sl[3])
                        name1 = d['molecules'][resname]['atoms']['name'][id1-1]
                        name2 = d['molecules'][resname]['atoms']['name'][id2-1]
                        name3 = d['molecules'][resname]['atoms']['name'][id3-1]
                        name4 = d['molecules'][resname]['atoms']['name'][id4-1]
                        func  = int(sl[4])
                        t     = float(sl[5])
                        k     = float(sl[6])

                        if len(sl) == 8:
                            n = int(sl[7])
                            d['molecules'][resname]['dihedrals'].append([
                                name1, name2, name3, name4, func, t, k, n])
                        elif len(sl) == 7:
                            d['molecules'][resname]['dihedrals'].append([
                                name1, name2, name3, name4, func, t, k, 0])

                        else:
                            assert 0 == 1, 'unsuccesful dihedral addition'

                    if read == 'constraints':
                        # V = 0.5 * k * (r - r0)^2
                        sl = line.split()
                        id1 = int(sl[0]); id2 = int(sl[1])
                        name1 = d['molecules'][resname]['atoms']['name'][id1-1]
                        name2 = d['molecules'][resname]['atoms']['name'][id2-1]
                        #if name1 > name2: name2 = name1
                        f = int(sl[2])
                        # f == 1 --> Just like bonds -> exclusions should be added
                        # f == 2 --> No connection, and so no exclusions, are generated for this interaction
                        l = float(sl[3]) #nm
                        d['molecules'][resname]['bonds'].append([name1, name2, l, self.Kc2b, f])

                    if read == 'exclusions':
                        sl = line.split()
                        for i in range(len(sl)-1):
                            id1 = int(sl[0]); id2 = int(sl[i+1])
                            name1 = d['molecules'][resname]['atoms']['name'][id1-1]
                            name2 = d['molecules'][resname]['atoms']['name'][id2-1]
                            d['molecules'][resname]['exclusions'].append([name1, name2])

                    if read == 'virtual_sites3':
                        sl = line.split()
                        id1 = int(sl[0])
                        id2 = int(sl[1])
                        id3 = int(sl[2])
                        id4 = int(sl[3])
                        f   = int(sl[4])
                        a   = float(sl[5])
                        b   = float(sl[6])

                        # f == 1 -> 3-body virtual site / virtual_lc3
                        # f == 4 -> 3-body virtual site (out)
                        
                        name1 = d['molecules'][resname]['atoms']['name'][id1-1]
                        name2 = d['molecules'][resname]['atoms']['name'][id2-1]
                        name3 = d['molecules'][resname]['atoms']['name'][id3-1]
                        name4 = d['molecules'][resname]['atoms']['name'][id4-1]

                        try:
                            # 3out has c. inverse nm to inverse A
                            d['molecules'][resname]['virtual_sites3'].append([name1, name2, name3, name4, f, a, b, 0.1 * float(sl[7])])
                        except:
                            d['molecules'][resname]['virtual_sites3'].append([name1, name2, name3, name4, f, a, b])


        self.martini = d




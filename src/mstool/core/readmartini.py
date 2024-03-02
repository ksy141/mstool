import os
import glob
pwd = os.path.dirname(os.path.realpath(__file__))

class ReadMartini:
    '''Read martini forcefields. The default forcefield files read are
    $mstool/FF/martini2.2/martini_v2.2.itp, 
    $mstool/FF/martini2.2/martini_v2.0_ions.itp, 
    $mstool/FF/martini2.2/martini_v2.2_proteins/proteins.itp, 
    $mstool/FF/martini2.2/martini_v2.2_proteins/proteins_HIS.itp,
    $mstool/FF/martini2.2/martini_v2.0_lipids_all_201506.itp,
    $mstool/FF/martini2.2/martini_lipids_add.itp,
    $mstool/FF/martini2.2/structures/*.itp.

    To add files on top of the default,
    use ff_add = [].

    To override the default files,
    use ff = [].

    Martini forcefields are native to gromacs, which uses
    kJ/mol, nm, kJ/mol/rad^2, ps, bar, unified atomic mass unit as its MD units.
    Although, for angles, please note that angles are written in degress not in radians,
    which will be internally converted to radians. This is due to user-friendliness.

    Parameters
    ----------
    ff : list
        A list of martini forcefield files. The default files read are
        $mstool/FF/martini2.2/martini_v2.2.itp, 
        $mstool/FF/martini2.2/martini_v2.0_ions.itp, 
        $mstool/FF/martini2.2/martini_v2.2_proteins/proteins.itp, 
        $mstool/FF/martini2.2/martini_v2.2_proteins/proteins_HIS.itp,
        $mstool/FF/martini2.2/martini_v2.0_lipids_all_201506.itp,
        $mstool/FF/martini2.2/martini_lipids_add.itp,
        $mstool/FF/martini2.2/structures/*.itp

        To override the default files, use this.
        e.g., ff = ['ff1.itp', 'ff2.itp']

        One can also provide the complete file.
        e.g., ff = ['topol.top']

    ff_add : list
        A list of *additional* martini forcefield files.
        If there are duplicate moleculetypes, the parameters are updated with the latest one.
        The default is [].

    define : dict
        Reads gromacs macro. Available keywords are:
        #ifdef, #ifndef, #else, #end, #define.
        e.g., define={'BILAYER_LIPIDHEAD_FC': '500.0', 'FLEXIBLE': 'True'}

    Kc2b : float
        The default is 10000.0.

    Attributes
    ----------
    ifiles : list
        List of files that were read

    martini : dict
        Contains all the information. 
        ``self.martini['molecules']`` saves ``moleculetypes`` field.
        ``self.martini['mols']`` saves ``molecules`` field.
        ``self.martini['nb']`` saves ``defaults`` field.
        ``self.martini['system']`` saves ``system`` field.
        ``self.martini['nbfunc']`` saves nonbonded potential (``LJ`` or ``Buckingham``)
        ``self.martini['energy']`` saves a potential form``


    define : dict
        Dictionary to indicate which gromacs macro was used when reading files.


    Examples
    --------
    >>> m = mstool.ReadMartini(ff_add=['add_molecules.itp'], 
    ...     define={'FLEXIBLE': 'True', 'MICELLE_LIPIDHEAD_FC: '300.0',
    ...             'VESICLE_LIPIDTAIL_R': '100.0', 
    ...             'BILAYER_LIPIDHEAD_FC': '200.0'})
    >>> m.martini['molecules']['POPC']
    >>> m.martini['define']
    >>>
    >>> m = mstool.ReadMartini('topol.top')
    >>> m.martini['molecules']
    >>> m.martini['mols']
    '''


    def __init__(self,
        ff = [], ff_add = [], 
        define = {}, Kc2b = 10000.0,
        constraint_to_bond=False):

        self.Kc2b = Kc2b
        self.ctb  = constraint_to_bond
        if not isinstance(ff, list): ff = [ff]
        if not isinstance(ff_add, list): ff_add = [ff_add]

        if ff == []:
            ff = [pwd + '/../FF/martini2.2/martini_v2.2.modified.itp',
                  pwd + '/../FF/martini2.2/martini_v2.0_ions.itp',
                  pwd + '/../FF/martini2.2/martini_v2.2_proteins/proteins.itp',
                  pwd + '/../FF/martini2.2/martini_v2.2_proteins/proteins_HIS.itp',
                  pwd + '/../FF/martini2.2/martini_v2.0_lipids_all_201506.itp',
                  pwd + '/../FF/martini2.2/martini_lipids_add.itp',
                  pwd + '/../FF/martini2.2/structures/*.itp']

        self.ifiles = []
        for ifile in [*ff, *ff_add]:
            if '*' in ifile:
                self.ifiles += glob.glob(ifile)
            else:
                self.ifiles.append(ifile)

        self.define = define
        self.collect()


    def collect(self):
        d = {}; read=None
        for ifile in self.ifiles:
            with open(ifile) as W:
                
                conditions = []
                for line in W.readlines():
                    if line.startswith(';'): continue
                    line = line.split(';')[0].strip()
                    if not line: continue


                    ### Add ifiles
                    if line.startswith('#include'):
                        self.ifiles.append(line.split()[1].replace('"', '').replace("'", ''))


                    ######### ifdef ifndef define #########
                    if line.startswith('#ifdef'):
                        key = line.split()[1]
                        if key in self.define.keys(): 
                            conditions.append([key, True])
                        else:
                            conditions.append([key, False])
                        continue

                    if line.startswith('#ifndef'):
                        key = line.split()[1]
                        if key in self.define.keys():
                            conditions.append([key, False])
                        else:
                            conditions.append([key, True])
                        continue

                    if line.startswith('#else'):
                        conditions[-1][1] = not conditions[-1][1]
                        continue

                    if line.startswith('#end'):
                        conditions.pop()
                        continue

                    if len(conditions) != 0:
                        if conditions[-1][1]:
                            pass
                        else:
                            continue

                    if line.startswith('#define'):
                        sl  = line.split()
                        self.define[sl[1]] = ' '.join(sl[2:])
                        continue

                    for key, value in self.define.items():
                        line = line.replace(key, value)
                    #######################################


                    if '[' in line and ']' in line:
                        read = line.split('[')[1].split(']')[0].strip()
                        continue

                    if read == 'system':
                        d['system'] = line.strip()

                    if read == 'molecules':
                        if 'mols' not in d.keys(): d['mols'] = []
                        sl = line.strip().split()
                        d['mols'].append([sl[0], int(sl[1])])

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
                        if int(sl[0]) == 1:
                            d['nbfunc'] = 'LJ'
                        if int(sl[0]) == 2:
                            d['nbfunc'] = 'Buckingham'

                        if int(sl[1]) == 1:
                            d['energy'] = 'C12/r^12 - C6/r^6'
                        if int(sl[1]) == 2:
                            d['energy'] = '4*eps*((sigma/r)^12 - (sigma/r)^6)'

                    if read == 'atomtypes':
                        if 'atomtypes' not in d.keys(): d['atomtypes'] = {}
                        
                        sl = line.split()
                        t = sl[0]
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
                                                             'resid': [], 'resname': [], 
                                                             'name':[], 'q':[], 'm': []}

                        d['molecules'][resname]['bonds']               = []
                        d['molecules'][resname]['angles']              = []
                        d['molecules'][resname]['dihedrals']           = []
                        d['molecules'][resname]['exclusions']          = []
                        d['molecules'][resname]['virtual_sites3']      = []
                        d['molecules'][resname]['position_restraints'] = []

                        d['molecules'][resname]['idx_bonds']               = []
                        d['molecules'][resname]['idx_angles']              = []
                        d['molecules'][resname]['idx_dihedrals']           = []
                        d['molecules'][resname]['idx_exclusions']          = []
                        d['molecules'][resname]['idx_virtual_sites3']      = []
                        d['molecules'][resname]['idx_position_restraints'] = []


                    if read == 'atoms':
                        sl = line.split()
                        d['molecules'][resname]['atoms']['id'].append(int(sl[0]))
                        d['molecules'][resname]['atoms']['type'].append(sl[1])
                        d['molecules'][resname]['atoms']['resid'].append(int(sl[2]))
                        d['molecules'][resname]['atoms']['resname'].append(sl[3])
                        d['molecules'][resname]['atoms']['name'].append(sl[4])
                        d['molecules'][resname]['atoms']['q'].append(float(sl[6]))

                        try:
                            mass = float(sl[7])
                        except:
                            mass = d['atomtypes'][sl[1]]['m']
                        d['molecules'][resname]['atoms']['m'].append(mass)
                        

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
                        d['molecules'][resname]['idx_bonds'].append([id1, id2, l, k])

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
                        d['molecules'][resname]['idx_angles'].append([id1, id2, id3, l, k, func])

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

                        try:
                            n = int(sl[7])
                        except:
                            n = 0

                        d['molecules'][resname]['dihedrals'].append([
                            name1, name2, name3, name4, func, t, k, n])
                        d['molecules'][resname]['idx_dihedrals'].append([
                            id1, id2, id3, id4, func, t, k, n])

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

                        if self.ctb:
                            d['molecules'][resname]['bonds'].append([name1, name2, l, self.Kc2b])
                            d['molecules'][resname]['idx_bonds'].append([id1, id2, l, self.Kc2b])

                        else:
                            d['molecules'][resname]['bonds'].append([name1, name2, l, self.Kc2b, f])
                            d['molecules'][resname]['idx_bonds'].append([id1, id2, l, self.Kc2b, f])

                    if read == 'exclusions':
                        sl = line.split()
                        for i in range(len(sl)-1):
                            id1 = int(sl[0]); id2 = int(sl[i+1])
                            name1 = d['molecules'][resname]['atoms']['name'][id1-1]
                            name2 = d['molecules'][resname]['atoms']['name'][id2-1]
                            d['molecules'][resname]['exclusions'].append([name1, name2])
                            d['molecules'][resname]['idx_exclusions'].append([id1, id2])

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
                            #vv = 0.1 * float(sl[7])
                            vv = float(sl[7])
                            d['molecules'][resname]['virtual_sites3'].append([name1, name2, name3, name4, f, a, b, vv])
                            d['molecules'][resname]['idx_virtual_sites3'].append([id1, id2, id3, id4, f, a, b, vv])

                        except:
                            d['molecules'][resname]['virtual_sites3'].append([name1, name2, name3, name4, f, a, b])
                            d['molecules'][resname]['idx_virtual_sites3'].append([id1, id2, id3, id4, f, a, b])


                    if read == 'position_restraints':
                        sl  = line.split()
                        id1 = int(sl[0])
                        name1 = d['molecules'][resname]['atoms']['name'][id1-1]
                        f   = int(float(sl[1]))

                        if f == 1:
                            kx = float(sl[2])
                            ky = float(sl[3])
                            kz = float(sl[4])
                            d['molecules'][resname]['position_restraints'].append([name1, f, kx, ky, kz])
                            d['molecules'][resname]['idx_position_restraints'].append([id1, f, kx, ky, kz])

                        else:
                            g = int(sl[2])
                            r = float(sl[3])
                            k = float(sl[4])
                            d['molecules'][resname]['position_restraints'].append([name1, f, g, r, k])
                            d['molecules'][resname]['idx_position_restraints'].append([id1, f, g, r, k])


        self.martini = d




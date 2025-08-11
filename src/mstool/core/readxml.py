import os
import numpy as np
import xml.etree.ElementTree as ET

pwd  = os.path.dirname(os.path.realpath(__file__))
aaff = [pwd + '/../FF/charmm36/charmm36.xml',
        pwd + '/../FF/charmm36/pip.xml',
        pwd + '/../FF/charmm36/water.xml',
        pwd + '/../FF/charmm36/chyo.xml']

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

class ReadXML:
    def __init__(self, ff=[], ff_add=[]):

        if not isinstance(ff, list):
            ff = [ ff ]

        if not isinstance(ff_add, list):
            ff_add = [ ff_add ]

        if ff == []:
            ff = aaff

        self.ff = [*ff, *ff_add]
        self.collect()


    def collect(self):
        self.RESI = {}
        self.TYPE = {}

        for toppar in self.ff:
            root = ET.parse(toppar).getroot()

            for atomtypes in root.findall('AtomTypes'):
                for atomtype in atomtypes.findall('Type'):
                    if atomtype.get("class") == 'DUM':
                        self.TYPE["DUM"] = 0
                    else:
                        self.TYPE[atomtype.get("class")] = E2Z[atomtype.get("element")]

            for residue in root.findall('Residues/Residue'):
                resname = residue.get('name')

                names   = []
                types   = []
                charges = []
                numbers = []
                for atom in residue.findall('Atom'):
                    names.append(atom.get('name'))
                    types.append(atom.get('type'))
                    charges.append(float(atom.get('charge')))
                    try:
                        numbers.append(self.TYPE[atom.get('type')])
                    except:
                        numbers.append(0)

                bonds = []
                for bond in residue.findall('Bond'):
                    b1 = bond.get('atomName1')
                    b2 = bond.get('atomName2')
                    bonds.append([b1, b2])

                data = {'names': np.array(names),
                        'types': np.array(types),
                        'numbers': np.array(numbers), 
                        'charges': np.array(charges),
                        'bonds': np.array(bonds)}

                if resname in self.RESI.keys():
                    raise Exception(f"residue {resname} already defined")
                else:
                    self.RESI[resname] = data

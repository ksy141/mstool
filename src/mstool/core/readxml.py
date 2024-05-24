import os
import numpy as np
import xml.etree.ElementTree as ET

pwd  = os.path.dirname(os.path.realpath(__file__))
aaff = [pwd + '/../FF/charmm36/charmm36.xml',
        pwd + '/../FF/charmm36/pip.xml',
        pwd + '/../FF/charmm36/water.xml',
        pwd + '/../FF/charmm36/chyo.xml']

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

        for toppar in self.ff:
            root = ET.parse(toppar).getroot()

            for residue in root.findall('Residues/Residue'):
                resname = residue.get('name')

                names   = []
                types   = []
                charges = []
                for atom in residue.findall('Atom'):
                    names.append(atom.get('name'))
                    types.append(atom.get('type'))
                    charges.append(float(atom.get('charge')))

                bonds = []
                for bond in residue.findall('Bond'):
                    b1 = bond.get('atomName1')
                    b2 = bond.get('atomName2')
                    bonds.append([b1, b2])

                data = {'names': np.array(names),
                        'types': np.array(types),
                        'charges': np.array(charges),
                        'bonds': np.array(bonds)}

                if resname in self.RESI.keys():
                    raise Exception(f"residue {resname} already defined")
                else:
                    self.RESI[resname] = data

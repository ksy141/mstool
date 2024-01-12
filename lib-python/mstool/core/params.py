from dataclasses import dataclass, fields

import os
pwd = os.path.dirname(os.path.realpath(__file__))

def realpath(ifile):
    return os.path.realpath(os.path.expanduser(ifile))

@dataclass
class Params:
    martiniff:   list[int] = field(default_factory=list, 
                                   default = [realpath(pwd + '/../../../FF/martini2.2/martini_v2.2.modified.itp'),
                                              realpath(pwd + '/../../../FF/martini2.2/martini_v2.0_ions.itp'),
                                              realpath(pwd + '/../../../FF/martini2.2/martini_v2.2_proteins/proteins.itp'),
                                              realpath(pwd + '/../../../FF/martini2.2/martini_v2.2_proteins/proteins_HIS.itp'),
                                              realpath(pwd + '/../../../FF/martini2.2/martini_v2.0_lipids_all_201506.itp'),
                                              realpath(pwd + '/../../../FF/martini2.2/martini_lipids_add.itp')])

    mapping:     list[int] = field(default_factory=list,
                                   default = [realpath(pwd + '/../../../mapping/martini.protein.c36m.dat'),
                                              realpath(pwd + '/../../../mapping/martini.lipid.c36.dat')])


    mapping_add: list[int] = field(default_factory=list)



    def showParams(self):
        for field in fields(self):
            print(field)




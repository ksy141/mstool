#!/usr/bin/env python3

import subprocess
import argparse

description = """This script dumps a DMS file"""

def parse_args():
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("dms", help='provide a dms file')
    return parser.parse_args()    

def main():
    dms = parse_args().dms

    cmdlines = ['select * from global_cell',
                'select * from particle',
                'select * from bond',
                'select * from exclusion']

    cmdlines.append("""SELECT p0, p1, r0, fc, constrained
        FROM stretch_harm_term INNER JOIN stretch_harm_param
        ON stretch_harm_term.param=stretch_harm_param.id""")
    
    cmdlines.append("""SELECT p0, p1, p2, theta0, fc, constrained
        FROM angle_harm_term INNER JOIN angle_harm_param
        ON angle_harm_term.param=angle_harm_param.id""")

    cmdlines.append("""SELECT p0, p1, p2, cos_theta0, fc
            FROM angle_harmcos_term INNER JOIN angle_harmcos_param
            ON angle_harmcos_term.param=angle_harmcos_param.id
            WHERE type='table' AND name='angle_harmcos_term'""")

    cmdlines.append("""SELECT p0, p1, p2, p3, phi0, fc0, fc1, fc2, fc3, fc4, fc5, fc6
        FROM dihedral_trig_term INNER JOIN dihedral_trig_param
        ON dihedral_trig_term.param=dihedral_trig_param.id""")

    cmdlines.append("""SELECT p0, p1, p2, p3, phi0, fc
        FROM improper_harm_term INNER JOIN improper_harm_param
        ON improper_harm_term.param=improper_harm_param.id""")

    cmdlines.append("""SELECT p0, p1, p2, p3, p4, p5, p6, p7, cmapid
        FROM torsiontorsion_cmap_term INNER JOIN torsiontorsion_cmap_param
        ON torsiontorsion_cmap_term.param=torsiontorsion_cmap_param.id""")

    cmdlines.append("""SELECT particle.id, charge, sigma, epsilon
        FROM particle INNER JOIN nonbonded_param
        ON particle.nbtype=nonbonded_param.id ORDER BY particle.id""")

    cmdlines.append("""SELECT p0, p1, aij, bij, qij
        FROM pair_12_6_es_term INNER JOIN pair_12_6_es_param
        ON pair_12_6_es_term.param=pair_12_6_es_param.id""")
    
    cmdlines.append("""SELECT * FROM nonbonded_combined_param""")

    cmdlines.append("""SELECT p0, p1, p2, c1
        FROM virtual_lc2_term INNER JOIN virtual_lc2_param
        ON virtual_lc2_term.param=virtual_lc2_param.id""")

    cmdlines.append("""SELECT p0, p1, p2, p3, c1, c2
        FROM virtual_lc3_term INNER JOIN virtual_lc3_param
        ON virtual_lc3_term.param=virtual_lc3_param.id""")

    cmdlines.append("""SELECT p0, p1, p2, p3, c1, c2, c3
        FROM virtual_out3_term INNER JOIN virtual_out3_param
        ON virtual_out3_term.param=virtual_out3_param.id""")
    
    cmdlines.append("""SELECT p0, x0, y0, z0, fcx, fcy, fcz 
        FROM posre_harm_term INNER JOIN posre_harm_param 
        ON posre_harm_term.param=posre_harm_param.id""")

    for cmdline in cmdlines:
        #subprocess.run(['sqlite3', '-header', dms, cmdline], stderr=subprocess.DEVNULL)
        subprocess.run(['sqlite3', '-header', dms, cmdline])

if __name__ == '__main__':
    main()


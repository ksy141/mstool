#!/usr/bin/python
import sys
try:
    import numpy as np
except:
    print('failed to import NumPy. Not installed?')
    sys.exit()


if '-h' in sys.argv or '--help' in sys.argv:
    print('')
    print('')   
    print('This script backmaps a coarse grained structure named NAMEcg.pdb to an all-atom one named NAMEbm.pdb.')
    print('It also looks up for a file NAME.pdb that is used as reference structure for histidine labeling and RMSD calculations.')
    print('')
    print('Usage: cg2aa.py NAME [OPTION]')
    print('')
    print('Options:')
    print('  -h or --help: shows this modest help')
    print('  -rmsd: writes RMSD comparison to NAME.rmsd (please note that reference structure is mandatory).')
    print('')
    sys.exit()

if len(sys.argv)>1:
    code = sys.argv[1].split('.')[0]
    cgfile = open(code+'cg.pdb', mode='r')
    try:
        pdbfile = open(code+'.pdb', mode='r')
    except:
        pdbfile = None
else:
    cgfile = open('cg.pdb', mode='r')
    try:
        pdbfile = open(code+'ch.pdb', mode='r')
    except:
        pdbfile = None

lines = [x.split() for x in cgfile.readlines()]
print(lines)

has_chain = int(lines[0][4][0].isalpha()) #checks if the first char of the 5th pos of the first line is a letter
#print lines[0],has_chain

if pdbfile:
    reflines = [x.split() for x in pdbfile.readlines()]
    ref_has_chain = int(reflines[0][4][0].isalpha())
    histidines = {int(y[4+ref_has_chain]): y[3] for y in reflines if y[0] == 'ATOM' and y[2] == 'CA' and y[3].startswith('HI')} 
    pdbfile.close()
else:
    histidines = None
grains = [(int(x[1]), x[2], x[3], int(x[4+has_chain])) for x in lines if len(x) > 1]
positions = np.array([map(np.float, [x[5+has_chain], x[6+has_chain], x[7+has_chain]]) for x in lines if len(x) > 1])
#starts = list(set([int(lines[0][4+has_chain])] + [int(x[4+has_chain]) for i,x in enumerate(lines[1:]) if len(lines[i-1]) == 1]))
cgfile.close()
starts = list({int(lines[0][4+has_chain])}.union({int(x[4+has_chain]) for a,x in zip(lines,lines[1:]) if len(a) == 1}))
# print lines[0]
# print lines[1:][:5]
# print len(lines[0])
#print starts
starts.sort()
ends =  [i-1 for i in starts[1:]] + [grains[-1][3]]
print( len(starts), 'chain'+('s' if len(starts) > 1 else '')+' detected.',)
print( 'They are' if len(starts) > 1 else 'It is', 'defined by the following residue numbers:')
print( ','.join([str(a)+'-'+str(object=b) for a,b in zip(starts, ends)]))

backbone = [x[0]-1 for x in grains if x[1] == 'CA']
tot_res = len(backbone) # grains[-1][-1]

formatString = 'ATOM %6i %s %s %5i     %7.3f %7.3f %7.3f %1.2f %1.2f'
translate_res={'AHV':'ALA','RHV':'ARG','NHV':'ASN','DHV':'ASP','CHV':'CYS','GHV':'GLY','QHV':'GLN','EHV':'GLU','IHV':'ILE','LHV':'LEU','KHV':'LYS','MHV':'MET','FHV':'PHE','PHV':'PRO','SHV':'SER','THV':'THR','WHV':'TRP','YHV':'TYR','VHV':'VAL','HEV':'HEM'}  #'HHV':'HIE'
mass = {'C': 12.0107,'H': 1.00794,'N': 14.0067,'O': 15.9994,'S': 32.065,'F':55.845}

cosg1 = 0.93544403088 #g1 = 20.7 CA CA CO
sing1 = 0.35347484377
cosg2 = 0.96944534989 #g2 = 14.2 CA CA N
sing2 = 0.24530738587
cos989 = 0.15471038629 #CSC
sin989 = 0.98795986576
cos120 = 0.5 #CCO
sin120 = 0.86602540378 #hexagon
cos108 = 0.30901699437 #pentagon
sin108 = 0.95105651629
cos126 = 0.58778525229 #pentagon external turn
sin126 = 0.80901699437
cos1274 = 0.60737583972 #Fe-N-C in HEM
sin1274 = 0.79441462053
cos1103 = 0.34693565157 #? in HEM
sin1103 = 0.93788893461
cos1249 = 0.57214587344 #CC-CB-CY in HEM
sin1249 = 0.82015187587
cos1166 = 0.44775908783 #CCN
sin1166 = 0.89415423683
cos1095 = 0.33380685923
sin1095 = 0.94264149109
cos111 = 0.35836794954
sin111 = 0.93295353482
cos117 = 0.45399049974 #CT - C - O2
sin117 = 0.89100652418
cosd = 0.57714519003 #d = half a tetrahedral angle
sind = 0.81664155516
cos45 = 0.70710678118

dctct = 1.526
dcc = 1.522
dcbcc = 1.44400 #HEM
dcbcy = 1.501 #HEM
dcycx = 1.34 #HEM
dco = 1.229
dcn = 1.335
dnc = 1.384 #HEM
dccna = 1.385
dcccv = 1.375
dcrna = 1.343
dcvnb = 1.394

dctn = 1.449
dcts = 1.810
dctoh = 1.410
dfech = 3.41315894736826
dfen = 2.01

def translate_residue(res, j):
    if res in translate_res.values():#i.e. it is valid name
        return res
    elif translate_res.has_key(res):
        return translate_res[res]
    else:
        return histidines[int(j)] if histidines else 'HIS'


#proy_ca = np.dot(positions[backbone[1:-1]] - positions[backbone[:-2]], co) # / np.dot(positions[backbone[1:-1]] - positions[backbone[:-2]], positions[backbone[1:-1]] - positions[backbone[:-2]])) * (positions[backbone[1:-2]] - positions[backbone[:-3]])

def append2pdb(atom_name, pos, residue, index):
    global pdb, atom_num
    atom_num += 1
    x,y,z = pos
    pdb.append(formatString % (atom_num, atom_name, translate_residue(residue,index), index, x,y,z,1.00,0.00)) #P



def completeTetra(v1,v2): #finds v3,v4 the other vertices of the tetra centered at zero
    v1 = v1 / np.linalg.norm(v1)
    v2 = v2 / np.linalg.norm(v2)
    nv = np.cross(v1,v2)
    nv = 0.8164965809277259 * nv / np.linalg.norm(nv)
    v = (-v1-v2)/2
    w1 = v + nv
    w2 = v - nv
#    w1 = w1 / np.linalg.norm(w1)
#    w2 = w2 / np.linalg.norm(w2)
#    print np.linalg.norm(v1 + v2 + w1 + w2),
    return	(w1,w2)

def companionVertex(v1,w): #finds the vertex of the tetrahedron of center zero in the direction of w
    v1 = v1 / np.linalg.norm(v1)
    w  = w - np.inner(v1, w)  * v1
    w  = w / np.linalg.norm(w)
    return - cos1095 * v1 + sin1095 * w

def turn(v1,v2,w, d, cos, sin): # turns at v2 in the direction of w
    v = v1- v2
    v = v / np.linalg.norm(v)
    w = w - v2
    w = w - np.inner(v, w)  * v
    w = w / np.linalg.norm(w)
    return v2 - d * cos * v + d * sin * w

def antibisect(v1,v2): #finds vector centered in zero in the direction of the angle given by v1,v2
    v1 = v1 / np.linalg.norm(v1)
    v2 = v2 / np.linalg.norm(v2)
    v = v1 + v2
    return - v / np.linalg.norm(v)

def normal2plane(v1,v2,v3):
    v = np.cross(v2-v1,v3-v1)
    return v / np.linalg.norm(v)

def hexagon(v,w): #constructs the hexagon ring starting from w that is bonded to v.
    pass
    

ca1 = positions[backbone[1:-1]] - positions[backbone[:-2]] # ca1[i] = CAi+1 - CAi
ca2 = positions[backbone[2:]] - positions[backbone[:-2]]   # ca2[i] = CAi+2 - CAi
p_vector = np.cross(ca1, ca2) #peptide plane vector.

#[ np.inner(ca1[i],ca2[i]) for i,_ in enumerate(ca1)]
atom_num = 0
pdb = []
for g in grains:
    i, grain_name, res, j = g #i grain number, j residue number
    if grain_name == 'CA': #New resid starts
        posCA = positions[i-1]
        if not j in ends and not j+1 in ends:
            e1 = ca1[j-1] / np.linalg.norm(ca1[j-1]) #normalized vector in the direction CA - CA
            v2 = p_vector[j-1] - np.inner(e1, p_vector[j-1])  * e1
            e2 = v2 / np.linalg.norm(v2)
            posC = posCA + dcc * cosg1 * e1 + dcc * sing1 * e2
            #posN = positions[i] - dctn * cosg2 * e1 - dctn * sing2 * e2
            if not j in starts:
                posN = posofnextN
            posofnextN = posC - dcn * cos1166 * e2 + dcn * sin1166 * e1
            posO = turn(posCA, posC, posofnextN, dco, cos120, -sin120)#posC + dco * cos120 * e2 + dco * sin120 * e1
        elif j in ends:
            posN = posofnextN
            v1 = posN - posCA
            v2 = positions[backbone[j-2]] - posCA
            posC = posCA + dctct * companionVertex(v2,v1)
            posO = turn(posCA, posC, posN, dco, cos120, -sin120)
            posOXT = turn(posCA, posC, posN, dco, cos120, sin120)
        elif j+1 in ends:
            posN = posofnextN
            posC = turn(posN,posCA,positions[backbone[j]],dctct,cos1095,sin1095) #
            posofnextN = turn(posCA,posC,positions[backbone[j]],dcn,cos1166,sin1166)#posC - dcn * cos1166 * e2 + dcn * sin1166 * e1)
            posO = turn(posCA, posC, posofnextN, dco, cos120, -sin120)
        if j in starts:
            if j != 1:
                pdb.append('TER')
            posCB = turn(posC,posCA,posCA + e2,dctct,cos1095,sin1095)
            v, _ = completeTetra(posC-posCA, posCB-posCA)
            posN = posCA + dcn * v
        else:
            v1 = posN - posCA
            v2 = posC - posCA
            v3, _ = completeTetra(v1, v2)
            posCB = posCA + dctct * v3
        append2pdb(' N  ', posN, res, j)
        append2pdb(' CA ', posCA, res, j)
        if not res in ['GHV', 'GLY']: #i.e. it has a CB
#            angle = np.dot(CBs[j] - posCA, posCB - posCA) / (np.linalg.norm(CBs[j] - posCA) * np.linalg.norm(posCB - posCA))
#            print j,translate_residue(res,j),  angle, np.linalg.norm(CBs[j]-posCB)
            append2pdb(' CB ', posCB, res, j)
        if res in ['AHV', 'GHV', 'ALA', 'GLY']: #i.e. there is not another coarse grain with CB
            append2pdb(' C  ', posC, res, j)
            append2pdb(' O  ', posO, res, j)
            if j in ends:
                append2pdb(' OXT', posOXT, res, j)

    if grain_name in ['AP','AP1','CB']:  # and j>1: ?
        if res in ['NHV','DHV','YHV', 'LHV','FHV', 'ASN', 'ASP', 'TYR', 'LEU', 'PHE']:  #AP = [' CB ','HB2','HB3',' CG '  ]
            v1 = posCA - posCB
            v2 = positions[i-1] - posCB
            posCG = posCB + dctct * companionVertex(v1,v2)
            append2pdb(' CG ', posCG, res, j)
        if res in ['HHV', 'WHV', 'HSD', 'TRP']:
            posAP = positions[i-1]
        if res in ['RHV', 'ARG']: #ARG, code neglets the effect of H in center of mass and a few suppositions more (CT-CT-CT-CT dihedral = 180).
            v1 = posCA - posCB
            v2 = positions[i-1] - posCB
            posCG = posCB + dctct * companionVertex(v1,v2)  #big aprox?
            posCD = posCG - v1 #(CT-CT-CT-CT dihedral = 180)
            append2pdb(' CG ', posCG, res, j)
            append2pdb(' CD ', posCD, res, j)
        if res in ['KHV', 'LYS']: #LYS, code neglets the effect of H in center of mass and a few suppositions more (CT-CT-CT-CT dihedral = 180).
            v1 = posCA - posCB
            v2 = positions[i-1] - posCB
            posCG = posCB + dctct * companionVertex(v1,v2)  #big aprox?
            posCD = posCG - v1 #(CT-CT-CT-CT dihedral = 180)
            append2pdb(' CG ', posCG, res, j)
            append2pdb(' CD ', posCD, res, j)
        if res in ['CHV', 'CYS']: #, AP = [' CB ','HB2','HB3',' SG ', HG  ] code neglets the effect of HG in center of mass
            v1 = posCA - posCB
            v2 = positions[i-1] - posCB
            posSG = posCB + dcts * companionVertex(v1,v2)  #params say 108.6
            append2pdb(' SG ', posSG, res, j)
            append2pdb(' C  ', posC, res, j)
            append2pdb(' O  ', posO, res, j)
            if j in ends:
                append2pdb(' OXT', posOXT, res, j)
        if res in ['QHV','EHV', 'MHV', 'GLN', 'GLU', 'MET']: #AP = ['CB','HB2','HB3','CG','HG2','HG3'] code neglets the effect of HG in center of mass
            v1 = posCA - posCB
            v2 = positions[i-1] - posCB
            posCG = posCB + dctct * companionVertex(v1,v2)
            append2pdb(' CG ', posCG, res, j)
        if res in ['VHV', 'VAL']:#posAP, CB,CA and HB are coplanar.
            v1 = posCA - posCB
            w = positions[i-1] - posCB
            v2 = companionVertex(v1,-w)
            v3, v4 = completeTetra(v1, v2)
            posCG1 = posCB + dctct * v3
            posCG2 = posCB + dctct * v4
            append2pdb(' CG1', posCG1, res, j)
            append2pdb(' CG2', posCG2, res, j)
            append2pdb(' C  ', posC, res, j)
            append2pdb(' O  ', posO, res, j)
            if j in ends:
                append2pdb(' OXT', posOXT, res, j)
        if res in ['THV', 'THR']: # simple aprox.
            v1 = posCA - posCB
            w = positions[i-1] - posCB
            v2 = companionVertex(v1,-w)
            v3, v4 = completeTetra(v1, v2)
            posOG1 = posCB + dco * v3
            posCG2 = posCB + dctct * v4
            append2pdb(' OG1', posOG1, res, j)
            append2pdb(' CG2', posCG2, res, j)
        if res in ['IHV', 'ILE']: #Sets position of CG2 #### chirality ?
            v1 = posCA - posCB
            w = positions[i-1] - posCB
            v2 = companionVertex(v1,w)
            posCG2 = posCB + dctct * v2
            append2pdb(' CG2', posCG2, res, j)
        if res in ['PHV', 'PRO']:
            posCG = turn(posCA, posCB, positions[i-1], dcc, cos108, sin108)
            posCD = turn(posCB, posCG, posN, dcc, cos108, sin108)
            append2pdb(' CD ', posCD, res, j)
            append2pdb(' CG ', posCG, res, j)
            append2pdb(' C  ', posC, res, j)
            append2pdb(' O  ', posO, res, j)
            if j in ends:
                append2pdb(' OXT', posOXT, res, j)
            
    if res in ['QHV', 'GLN'] and grain_name in ['PL', 'CG']: 
        posPL = positions[i-1]
        posCD =  posCG + dctct * companionVertex(posCB-posCG, posPL-posCG)
        v1 = posCG - posCD
        v1 = v1 / np.linalg.norm(v1)
        v2 = posPL - posCD
        v2 = v2 - np.inner(v1, v2)  * v1
        v2 = v2 / np.linalg.norm(v2)
        posOE1 = posCD - dco * cos120 * v1 - dco * sin120 * v2   #NH2 beats O #not really working? #SHOULD they be coplanar w/ CG and CD ?
        posNE2 = posCD - dcn * cos1166 * v1 + dcn * sin1166 * v2
        append2pdb(' CD ', posCD, res, j)
        append2pdb(' OE1', posOE1, res, j)
        append2pdb(' NE2', posNE2, res, j)
        append2pdb(' C  ', posC, res, j)
        append2pdb(' O  ', posO, res, j)
        if j in ends:
            append2pdb(' OXT', posOXT, res, j)
    if res in ['EHV', 'GLU'] and grain_name in ['NE', 'CG']:
        posNE = positions[i-1]
        posCD = posCG + dctct * companionVertex(posCB-posCG,posNE-posCG)
        v1 = posCG - posCD
        v1 = v1 / np.linalg.norm(v1)
        v2 = posNE - posCD
        v2 = v2 - np.inner(v1, v2)  * v1
        v2 = v2 / np.linalg.norm(v2)
        posOE1 = posCD - dco * cos120 * v1 - dco * sin120 * v2
        posOE2 = posCD - dco * cos120 * v1 + dco * sin120 * v2
        append2pdb(' CD ', posCD, res, j)
        append2pdb(' OE1', posOE1, res, j)
        append2pdb(' OE2', posOE2, res, j)
        append2pdb(' C  ', posC, res, j)
        append2pdb(' O  ', posO, res, j)
        if j in ends:
            append2pdb(' OXT', posOXT, res, j)

    if res in ['DHV', 'ASP'] and grain_name == ['NE', 'CG']:
        v1 = posCB - posCG
        v1 = v1 / np.linalg.norm(v1)
        v2 = positions[i-1] - posCG
        v2 = v2 - np.inner(v1, v2)  * v1
        v2 = v2 / np.linalg.norm(v2)
        posOD1 = posCG - dco * cos120 * v1 - dco * sin120 * v2   #In theory this shoudn't work because CB CG and NE are colinear.
        posOD2 = posCG - dco * cos120 * v1 + dco * sin120 * v2
        append2pdb(' OD1', posOD1, res, j)
        append2pdb(' OD2', posOD2, res, j)
        append2pdb(' C  ', posC, res, j)
        append2pdb(' O  ', posO, res, j)
        if j in ends:
            append2pdb(' OXT', posOXT, res, j)

    if res == ['NHV', 'ASN'] and grain_name == ['PL', 'CG']:
        v1 = posCB - posCG
        v1 = v1 / np.linalg.norm(v1)
        v2 = positions[i-1] - posCG
        v2 = v2 - np.inner(v1, v2)  * v1
        v2 = v2 / np.linalg.norm(v2)
        posOD1= posCG - dco * cos120 * v1 - dco * sin120 * v2   #NH2 beats O
        posND2 = posCG - dcn * cos1166 * v1 + dcn * sin1166 * v2
        append2pdb(' OD1', posOD1, res, j)
        append2pdb(' ND2', posND2, res, j)
        append2pdb(' C  ', posC, res, j)
        append2pdb(' O  ', posO, res, j)
        if j in ends:
            append2pdb(' OXT', posOXT, res, j)

    if res in ['SHV', 'SER'] and grain_name in  ['PL', 'CB']: #AP = [' CB ','HB2','HB3',' OG ', HG  ] code neglets the effect of HG in center of mass
        posOG = posCB + dctoh * companionVertex(posCA - posCB,positions[i-1] - posCB) #or should it be dco?
        append2pdb(' OG ', posOG, res, j)
        append2pdb(' C  ', posC, res, j)
        append2pdb(' O  ', posO, res, j)
        if j in ends:
            append2pdb(' OXT', posOXT, res, j)
        
    if res in ['IHV', 'ILE'] and grain_name in ['AP2', 'CG']:
        posCG1 = posCB + dctct * completeTetra(posCG2-posCB, posCA-posCB)[0] #check chirality
        posCD = posCG1 + dctct * companionVertex(posCB-posCG1,  positions[i-1]-posCG1)
        append2pdb(' CD1', posCD, res, j)
        append2pdb(' CG1', posCG1, res, j)
        append2pdb(' C  ', posC, res, j)
        append2pdb(' O  ', posO, res, j)
        if j in ends:
            append2pdb(' OXT', posOXT, res, j)
        
    if res in ['MHV','MET'] and grain_name in ['AP2','CG']:
        posSD = posCG + dcts * companionVertex(posCB-posCG,  positions[i-1]-posCG)
#        v1 = posCG - posSD
#        v1 = v1 / np.linalg.norm(v1)
        v2 = positions[i-1] - posSD
#        v2 = v2 - np.inner(v1, v2)  * v1
        v2 = v2 / np.linalg.norm(v2)
        posCE = posSD + dcts * v2 # cos989 * v1 + dcts * sin989 * v2   #Met angle = 98.9
        append2pdb(' SD ', posSD, res, j)
        append2pdb(' CE ', posCE, res, j)
        append2pdb(' C  ', posC, res, j)
        append2pdb(' O  ', posO, res, j)
        if j in ends:
            append2pdb(' OXT', posOXT, res, j)

    if res in ['KHV', 'LYS'] and grain_name in ['PS','CG']:
        posCE = posCD + dctct* companionVertex(posCG - posCD,  positions[i-1] - posCD)
        posNZ = posCE + dctn * companionVertex(posCD - posCE, positions[i-1] - posCE)
        append2pdb(' CE ', posCE, res, j)
        append2pdb(' NZ ', posNZ, res, j)
        append2pdb(' C  ', posC, res, j)
        append2pdb(' O  ', posO, res, j)
        if j in ends:
            append2pdb(' OXT', posOXT, res, j)
       
    if res in ['LHV','LEU'] and grain_name in ['AP2', 'CG']:
        v1, v2 = completeTetra(posCB-posCG, -companionVertex(posCB-posCG,  positions[i-1]-posCG)) #check chirality
        posCD1 = posCG + dctct * v1
        posCD2 = posCG + dctct * v2
        append2pdb(' CD1', posCD1, res, j)
        append2pdb(' CD2', posCD2, res, j)
        append2pdb(' C  ', posC, res, j)
        append2pdb(' O  ', posO, res, j)
        if j in ends:
            append2pdb(' OXT', posOXT, res, j)

    if res in ['THV', 'THR'] and grain_name in ['PL','CB']:
        v1, v2 = completeTetra(posCA-posCB, -companionVertex(posCA-posCB,  positions[i-1]-posCB)) #check chirality
        posOG1 = posCB + dco * v1
        posCG2 = posCB + dctct * v2
        append2pdb(' OG1', posOG1, res, j)
        append2pdb(' CG2', posCG2, res, j)
        append2pdb(' C  ', posC, res, j)
        append2pdb(' O  ', posO, res, j)
        if j in ends:
            append2pdb(' OXT', posOXT, res, j)

    if res in ['RHV', 'ARG'] and grain_name in ['PS', 'CG']:
        posPS = positions[i-1]
        posNE = posCD + dctn * companionVertex(posCG-posCD,  posPS-posCD) #assumes CoM points aprox in the direction of NE
        append2pdb(' NE ', posNE, res, j)
    if res in ['RHV', 'ARG'] and grain_name in ['PL', 'CD']:
        posPL = positions[i-1]
        posCZ = turn(posCD,posNE,posPL,dcn,cos120,sin120)#should it be 123.2?
        posNH2= turn(posNE,posCZ,posPL,dcn,cos120,sin120)
        posNH1= turn(posNE,posCZ,posPL,dcn,cos120,-sin120)#posNH1= turn(posNE,posCZ,posPS,dcn,cos120,sin120)
        append2pdb(' CZ ', posCZ, res, j)
        append2pdb(' NH1', posNH1, res, j)
        append2pdb(' NH2', posNH2, res, j)
        append2pdb(' C  ', posC, res, j)
        append2pdb(' O  ', posO, res, j)
        if j in ends:
            append2pdb(' OXT', posOXT, res, j)

    if res == ['FHV', 'PHE'] and grain_name in ['AP2', 'CG']:
        posAP2 = positions[i-1] #CoM points aprox in the direction of CD1
    if res == ['FHV', 'PHE'] and grain_name in ['AP3', 'CD']:
        posAP3 = positions[i-1]
        posMC = (posAP2 * (2*mass['C'] + 2*mass['H']) + posAP3 * (3*mass['C'] + 3*mass['H'])) / (5*mass['C'] + 5*mass['H'])
        v = posMC - posCG
        v = v / np.linalg.norm(v)
        posCZ  = posCG + 2 * dcc * v
        posOH = posCZ + dco * v
        posCE1 = turn(posCZ + dcc * v,posCZ, posAP2,dcc,cos120,sin120)
        posCE2 = turn(posCZ + dcc * v,posCZ, posAP3,dcc,cos120,sin120)
        posCD1 = turn(posCZ,posCE1,posCG,dcc,cos120,sin120)
        posCD2 = turn(posCZ,posCE2,posCG,dcc,cos120,sin120)
        append2pdb(' CD1', posCD1, res, j)
        append2pdb(' CE1', posCE1, res, j)
        append2pdb(' CD2', posCD2, res, j)
        append2pdb(' CE2', posCE2, res, j)
        append2pdb(' CZ ', posCZ, res, j)
        append2pdb(' C  ', posC, res, j)
        append2pdb(' O  ', posO, res, j)
        if j in ends:
            append2pdb(' OXT', posOXT, res, j)
    if res in ['YHV', 'TYR'] and grain_name in ['AP2', 'CG']:
        posAP2 = positions[i-1] #CoM points aprox in the direction of CD1
    if res in ['YHV', 'TYR'] and grain_name in ['AP3', 'CD']:
        posAP3 = positions[i-1] #assumes CoM points aprox in the direction of NE
    if res in ['YHV', 'TYR'] and grain_name in ['PL', 'CE']:        
        posPL = positions[i-1]
        v = posPL - posCG
        v = v / np.linalg.norm(v)
        posCZ  = posCG + 2 * dcc * v
        posOH = posCZ + dco * v
        posCE1 = turn(posOH,posCZ, posAP2,dcc,cos120,sin120)
        posCE2 = turn(posOH,posCZ, posAP3,dcc,cos120,sin120)
        posCD1 = turn(posCZ,posCE1,posCG,dcc,cos120,sin120)
        posCD2 = turn(posCZ,posCE2,posCG,dcc,cos120,sin120)
        append2pdb(' CD1', posCD1, res, j)
        append2pdb(' CE1', posCE1, res, j)
        append2pdb(' CD2', posCD2, res, j)
        append2pdb(' CE2', posCE2, res, j)
        append2pdb(' CZ ', posCZ, res, j)
        append2pdb(' OH ', posOH, res, j)
        append2pdb(' C  ', posC, res, j)
        append2pdb(' O  ', posO, res, j)
        if j in ends:
            append2pdb(' OXT', posOXT, res, j)

    if res in ['HHV', 'HSD'] and grain_name in ['PL1', 'CG']: #HIS
        posPL1 = positions[i-1]
    if res in ['HHV', 'HSD'] and grain_name in ['PL2', 'CD']: #HIS
        posPL2 = positions[i-1] #assumes CoM points aprox in the direction of CD2, NE2
        v1 = posCA - posCB
        v2 = (posPL1 - posCB + posPL2 - posCB)/2
#        v2 = v2 / np.linalg.norm(v2)
        posCG = posCB + dctct * companionVertex(v1,v2)
        w = posCG - posCB
        w = w / np.linalg.norm(w)
        nv = v2 - np.inner(v2,w) * w
        posND1 = turn(posCB, posCG, posPL1-nv,dccna,cos126,sin126)
        posCD2 = turn(posCB, posCG, posPL2-nv,dcccv,cos126,sin126)
        posCE1 = turn(posCG, posND1,posPL2-nv,dcrna,cos108,sin108)
        posNE2 = turn(posCG, posCD2,posPL1-nv,dcvnb,cos108,sin108)
        append2pdb(' CG ', posCG, res, j)
        append2pdb(' ND1', posND1, res, j)
        append2pdb(' CE1', posCE1, res, j)
        append2pdb(' CD2', posCD2, res, j)
        append2pdb(' NE2', posNE2, res, j)
        append2pdb(' C  ', posC, res, j)
        append2pdb(' O  ', posO, res, j)
        if j in ends:
            append2pdb(' OXT', posOXT, res, j)

    if res in ['WHV', 'TRP'] and grain_name in ['AP1', 'CB']:
        posAP1 = positions[i-1] # ['CB','HB2','HB3','CG']
    if res in ['WHV', 'TRP'] and grain_name in ['PL', 'CG']:
        posPL = positions[i-1] #['CD1','HD1','NE1','HE1','CE2'] 
    if res in ['WHV', 'TRP'] and grain_name in ['AP2', 'CD']:
        posAP2 = positions[i-1] # ['CD2','CE3','HE3','CZ3','HZ3']
    if res in ['WHV', 'TRP'] and grain_name in ['AP3', 'CE']:
        posAP3 = positions[i-1] # ['CZ2','HZ2','CH2','HH2']
        nv = normal2plane(posPL, posAP2, posAP3)
        w = posAP - posCB
        w = w - np.inner(w,nv) * nv 
        posCG = turn(posCA, posCB,posAP1,dcc,cos108,sin108)# dcc * companionVertex(posCA - posCB,  w)
        posMC = (posAP2 + posAP3 + posPL)/3
        posCD2 = turn(posCB, posCG, posMC,dcccv,cos126,sin126)
        posCD1 = turn(posCB, posCG, posMC,dcccv,cos126,-sin126)
        posNE1 = turn(posCG, posCD1,posMC,dcvnb,cos108,sin108)
        posCE2 = turn(posCD1, posNE1, posMC,dcc,cos108,sin108)
        posCE3 = turn(posCE2, posCD2, posMC,dcc,cos120,sin120)
        posCZ3 = turn(posCD2, posCE3, posMC,dcc,cos120,sin120)
        posCZ2 = turn(posCD2, posCE2, posMC,dcc,cos120,sin120)
        posCH2 = turn(posCE2, posCZ2, posCZ3,dcc,cos120,sin120)
        append2pdb(' CG ', posCG, res, j)
        append2pdb(' CD1', posCD1, res, j)
        append2pdb(' NE1', posNE1, res, j)
        append2pdb(' CD2', posCD2, res, j)
        append2pdb(' CE2', posCE2, res, j)
        append2pdb(' CE3', posCE3, res, j)
        append2pdb(' CZ3', posCZ3, res, j)
        append2pdb(' CZ2', posCZ2, res, j)
        append2pdb(' CH2', posCH2, res, j)
        append2pdb(' C  ', posC, res, j)
        append2pdb(' O  ', posO, res, j)
        if j in ends:
            append2pdb(' OXT', posOXT, res, j)

#HEM
    if res == 'HEV': #HEM
        if grain_name == 'PL':
            posFE = positions[i-1]
        if grain_name == 'AP1':
            posAP1 = positions[i-1]
        if grain_name == 'AP2':
            posAP2 = positions[i-1]
        if grain_name == 'AP3':
            posAP3 = positions[i-1]
        if grain_name == 'AP4':
            posAP4 = positions[i-1]
        if grain_name == 'AP5':
            posAP5 = positions[i-1]
        if grain_name == 'AP6':
            posAP6 = positions[i-1]
        if grain_name == 'AP7':
            posAP7 = positions[i-1]
        if grain_name == 'AP8':
            posAP8 = positions[i-1]
        if grain_name == 'AP9':
            posAP9 = positions[i-1]
        if grain_name == 'AP10':
            posAP10 = positions[i-1]
        if grain_name == 'AP11':
            posAP11 = positions[i-1]
        if grain_name == 'AP12':
            posAP12 = positions[i-1]
        if grain_name == 'AP13':
            posAP13 = positions[i-1]
        if grain_name == 'AP14':
            posAP14 = positions[i-1]
        if grain_name == 'NE1':
            posNE1 = positions[i-1]
        if grain_name == 'NE2':
            posNE2 = positions[i-1]

            vc = (posAP6-posFE) / np.linalg.norm(posAP6-posFE)
            #va = (posAP2-posFE) / np.linalg.norm(posAP2-posFE)
            #print va-vc

            posCHC = posFE + dfech * vc
            posCHA = posFE - dfech * vc

            v1 = v1 / np.linalg.norm(v1)
            w = posAP9 - posFE
            w = w - np.inner(vc, w)  * vc
            w = w / np.linalg.norm(w)

            posCHD = posFE + dfech * w
            posCHB = posFE - dfech * w

            posNC = posFE + dfen * ( cos45 * w + cos45 * vc)
            posNB = posFE + dfen * (-cos45 * w + cos45 * vc)
            posND = posFE + dfen * ( cos45 * w - cos45 * vc)
            posNA = posFE + dfen * (-cos45 * w - cos45 * vc)

            posC1A = turn(posFE,posNA, posCHA,dnc,cos1274,sin1274)
            posC4A = turn(posFE,posNA, posCHB,dnc,cos1274,sin1274)
            posC1B = turn(posFE,posNB, posCHB,dnc,cos1274,sin1274)
            posC4B = turn(posFE,posNB, posCHC,dnc,cos1274,sin1274)
            posC1C = turn(posFE,posNC, posCHC,dnc,cos1274,sin1274)
            posC4C = turn(posFE,posNC, posCHD,dnc,cos1274,sin1274)
            posC1D = turn(posFE,posND, posCHD,dnc,cos1274,sin1274)
            posC4D = turn(posFE,posND, posCHA,dnc,cos1274,sin1274)

            posC2A = turn(posNA,posC1A, posCHA,dcbcc,cos1103,-sin1103)
            posC3A = turn(posNA,posC4A, posCHB,dcbcc,cos1103,-sin1103)
            posC2B = turn(posNB,posC1B, posCHB,dcbcc,cos1103,-sin1103)
            posC3B = turn(posNB,posC4B, posCHC,dcbcc,cos1103,-sin1103)
            posC2C = turn(posNC,posC1C, posCHC,dcbcc,cos1103,-sin1103)
            posC3C = turn(posNC,posC4C, posCHD,dcbcc,cos1103,-sin1103)
            posC2D = turn(posND,posC1D, posCHD,dcbcc,cos1103,-sin1103)
            posC3D = turn(posND,posC4D, posCHA,dcbcc,cos1103,-sin1103)

            posCAB = posC3B + dcbcy * antibisect(posC4B-posC3B,posC2B-posC3B) #turn(posC4B,posC3B, posAP8,dcbcy,cos1249,sin1249)
            posCBB = turn(posC3B,posCAB, posAP8,dcycx,cos120,sin120)
            posCMB = turn(posC1B,posC2B, posAP7,dcbcy,cos1249,sin1249)

            posCAC = posC3C + dcbcy * antibisect(posC4C-posC3C,posC2C-posC3C) #posC3C + dcbcy * antibisect(posC4C-posC3C,posC2C-posC3C) #turn(posC4C,posC3C,posAP11,dcbcy,cos1249,sin1249)
            posCBC = turn(posC3C,posCAC,posAP11,dcycx,cos120,sin120)
            posCMC = turn(posC1C,posC2C,posAP10,dcbcy,cos1249,sin1249)

            posCAA = posC2A + dcbcy * antibisect(posC3A-posC2A,posC1A-posC2A) #turn(posC1A,posC2A, posAP4,dcbcy,cos1249,sin1249)
            posCMA = turn(posC4A,posC3A, posAP3,dcbcy,cos1249,sin1249)
            posCBA = posCAA + dctct * companionVertex(posC2A - posCAA, posAP4 - posCAA)
            posCGA = posCBA + dcc * companionVertex(posCAA - posCBA, posNE1 - posCBA)
            posO1A = turn(posCBA,posCGA,posAP4,dco,cos117,-sin117)
            posO2A = turn(posCBA,posCGA,posAP4,dco,cos117,sin117)
            posCAD = posC3D + dcbcy * antibisect(posC4D-posC3D,posC2D-posC3D) #turn(posC4D,posC3D,posAP14,dcbcy,cos1249,sin1249)
            posCMD = turn(posC1D,posC2D,posAP13,dcbcy,cos1249,sin1249)
            posCBD = posCAD + dctct * companionVertex(posC3D - posCAD, posAP14 - posCAD)
            posCGD = posCBD + dcc * companionVertex(posCAD - posCBD, posNE2 - posCBD)
            posO1D = turn(posCBD,posCGD,posAP14,dco,cos117,-sin117)
            posO2D = turn(posCBD,posCGD,posAP14,dco,cos117,sin117)

            pdb.append('TER') #is it mandatory?
            append2pdb(' FE ', posFE, res, j)
            append2pdb(' NA ', posNA, res, j)
            append2pdb(' C1A', posC1A, res, j)
            append2pdb(' C2A', posC2A, res, j)
            append2pdb(' C3A', posC3A, res, j)
            append2pdb(' C4A', posC4A, res, j)
            append2pdb(' CHA', posCHA, res, j)
            append2pdb(' CHB', posCHB, res, j)
            append2pdb(' NB ', posNB, res, j)
            append2pdb(' C1B', posC1B, res, j)
            append2pdb(' C2B', posC2B, res, j)
            append2pdb(' C3B', posC3B, res, j)
            append2pdb(' C4B', posC4B, res, j)
            append2pdb(' CHC', posCHC, res, j)
            append2pdb(' NC ', posNC, res, j)
            append2pdb(' C1C', posC1C, res, j)
            append2pdb(' C2C', posC2C, res, j)
            append2pdb(' C3C', posC3C, res, j)
            append2pdb(' C4C', posC4C, res, j)
            append2pdb(' CHD', posCHD, res, j)
            append2pdb(' ND ', posND, res, j)
            append2pdb(' C1D', posC1D, res, j)
            append2pdb(' C2D', posC2D, res, j)
            append2pdb(' C3D', posC3D, res, j)
            append2pdb(' C4D', posC4D, res, j)

            append2pdb(' CAA', posCAA, res, j)
            append2pdb(' CBA', posCBA, res, j)
            append2pdb(' CGA', posCGA, res, j)
            append2pdb(' O1A', posO1A, res, j)
            append2pdb(' O2A', posO2A, res, j)
            append2pdb(' CMA', posCMA, res, j)
            append2pdb(' CAD', posCAD, res, j)
            append2pdb(' CAB', posCAB, res, j)
            append2pdb(' CBB', posCBB, res, j)
            append2pdb(' CMB', posCMB, res, j)
            append2pdb(' CAC', posCAC, res, j)
            append2pdb(' CBC', posCBC, res, j)
            append2pdb(' CMC', posCMC, res, j)
            append2pdb(' CBD', posCBD, res, j)
            append2pdb(' CGD', posCGD, res, j)
            append2pdb(' O1D', posO1D, res, j)
            append2pdb(' O2D', posO2D, res, j)
            append2pdb(' CMD', posCMD, res, j)


pdb.append('END')


ofile = open(code+'bm.pdb', 'w')
ofile.write('\n'.join(pdb))
ofile.close()


if '-rmsd' in sys.argv:
    if pdbfile:
        bmplines = [x.split() for x in pdb if not x.startswith('TER') and not x.startswith('END')]
        pdblines = [x for x in reflines if x != ['TER'] and x != ['END']]
    
        residues = 'ALA ARG ASN ASP CYS GLY GLN GLU ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL HI HEM'.split()
    
    
        cmpfile = open(code+'.rmsd', mode='w')
    
    #    tot_res = max([int(x[4]) for x in pdblines if x[3] != 'HEM'])
        tot_res = len([1 for _, grain_name, _, _ in grains if grain_name == 'CA'])
        pdbAtoms = {tuple(x[2:5]) : map(np.float, [x[5], x[6], x[7]]) for x in pdblines if not x[2].startswith('H')} # and x[2] != 'OXT']
        bmpAtoms = {tuple(x[2:5]) : map(np.float, [x[5], x[6], x[7]]) for x in bmplines if not x[2].startswith('H')}  # and x[2] != 'OXT']
        commonAtoms = list(set(pdbAtoms.keys()).intersection(set(bmpAtoms.keys())))
        Ns = list({int(x[2]) for x in commonAtoms if x[1] != 'HEM'})
        Ns.sort()
        
        pdbpositions = np.array([pdbAtoms[x] for x in commonAtoms])
        bmppositions = np.array([bmpAtoms[x] for x in commonAtoms])
        pdbpositionsBB = np.array([pdbAtoms[x] for x in commonAtoms if x[0] in ['C', 'CA', 'O', 'N']])
        bmppositionsBB = np.array([bmpAtoms[x] for x in commonAtoms if x[0] in ['C', 'CA', 'O', 'N']])
        pdbpositionsbyR = {res : np.array([pdbAtoms[x] for x in commonAtoms if not x[0] in ['C', 'CA', 'O', 'N','CB'] and x[1].startswith(res)]) for res in residues}
        bmppositionsbyR = {res : np.array([bmpAtoms[x] for x in commonAtoms if not x[0] in ['C', 'CA', 'O', 'N','CB'] and x[1].startswith(res)]) for res in residues}
        pdbpositionsBBbyN = {n : np.array([pdbAtoms[x] for x in commonAtoms if x[0] in ['C', 'CA', 'O', 'N'] and int(x[2]) == n]) for n in Ns}
        bmppositionsBBbyN = {n : np.array([bmpAtoms[x] for x in commonAtoms if x[0] in ['C', 'CA', 'O', 'N'] and int(x[2]) == n]) for n in Ns}
        pdbpositionsSCbyN = {n : np.array([pdbAtoms[x] for x in commonAtoms if not x[0] in ['C', 'CA', 'O', 'N'] and int(x[2]) == n]) for n in Ns}
        bmppositionsSCbyN = {n : np.array([bmpAtoms[x] for x in commonAtoms if not x[0] in ['C', 'CA', 'O', 'N'] and int(x[2]) == n]) for n in Ns}
        
        
        #pdbpositionsBB = np.array([map(np.float, [x[5], x[6], x[7]]) for x in pdblines if x[2] in ['C', 'CA', 'O', 'N']])
        #bmppositionsBB = np.array([map(np.float, [x[5], x[6], x[7]]) for x in bmplines if x[2] in ['C', 'CA', 'O', 'N']])
        #CBs = {int(y[4]) : np.array([float(y[5]),float(y[6]),float(y[7])]) for y in pdblines if y[2] == 'CB'}
        #print pdbpositionsbyR['HI']

        print('RMSD calculation between backmapped structue ('+code+'bm.pdb) and reference structure ('+code+'.pdb) (in A).',)
        print('Backbone: '+str(np.sqrt(np.sum(np.square(pdbpositionsBB - bmppositionsBB)) / len(pdbpositionsBB)))+', total: '+str(np.sqrt(np.sum(np.square(pdbpositions - bmppositions)) / len(pdbpositions)))+'.')
        print('More detailed information (i.e. by residue type and number) is written into file '+code+'.rmsd.')
        
        cmpfile.write('RMSD between backmapped structue ('+code+'bm.pdb) and  reference structure ('+code+'.pdb) in An\n')
        cmpfile.write('Backbone:'+str(np.sqrt(np.sum(np.square(pdbpositionsBB - bmppositionsBB)) / len(pdbpositionsBB)))+'\n')
        cmpfile.write('Total:'+str(np.sqrt(np.sum(np.square(pdbpositions - bmppositions)) / len(pdbpositions)))+'\n')
        cmpfile.write('By residue type:\n')
        for res in residues:
            cmpfile.write(res+':'+str(np.sqrt(np.sum(np.square(pdbpositionsbyR[res] - bmppositionsbyR[res])) / len(pdbpositionsbyR[res]) if len(pdbpositionsbyR[res]) else 0))+'\n')
        cmpfile.write('By residue number (backbone and total):\n')    
        for n in Ns:
            cmpfile.write(str(n)+':'+str(np.sqrt(np.sum(np.square(pdbpositionsBBbyN[n] - bmppositionsBBbyN[n])) / len(pdbpositionsBBbyN[n])))+':'+str(np.sqrt(np.sum(np.square(pdbpositionsSCbyN[n] - bmppositionsSCbyN[n])) / len(pdbpositionsSCbyN[n]) if len(pdbpositionsSCbyN[n]) else 0))+'\n')
    else:
        print('Cannot provide RMSD comparison without a reference structure.')



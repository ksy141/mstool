# https://iopscience.iop.org/article/10.1088/1367-2630/14/2/023063#:~:text=In%20this%20plot%2C%20the%20circular,be%20H%20%3D%202πR%2Fω.
# x = (R + a cos(wt)) cos t
# y = (R + a cos(wt)) sin t
# z = a sin(wt)
# w = number of helix turns = 1 turn by 3.6 a.a. = 5 turns by 18 a.a. = 10 pi by 18 a.a.
# H = translation per turn = 2 pi R / w = 3.6 residues per turn * 1.5 translation = 5.4 A
# R = 0.8594 w = 0.2387 * # of amino acid

### alpha helix
# i'th C=O -> (i+4)'th N-H hydrogen bond

import mstool
import numpy as np

### the only parameter
N = 90
ifile = 'step1.pdb'


### others will be calculated from N
R = 0.2387 * N #A
t = np.linspace(0, 2 * np.pi, N + 1)[:-1]
dt = t[1] - t[0]
a = 2.3
w = int(N / 3.6)
print(f"{N} a.a.; {int(N/3.6)} turns; {R} A radius")


def circularhelix(R, a, w, t):
    x = (R + a * np.cos(w * t)) * np.cos(t)
    y = (R + a * np.cos(w * t)) * np.sin(t)
    z = a * np.sin(w * t)
    return x, y, z


def completeTetra(v1,v2): 
    #finds v3,v4 the other vertices of the tetra centered at zero
    v1 = v1 / np.linalg.norm(v1)
    v2 = v2 / np.linalg.norm(v2)
    nv = np.cross(v1,v2)
    nv = 0.8164965809277259 * nv / np.linalg.norm(nv)
    v = (-v1-v2)/2
    w1 = v + nv
    w2 = v - nv
    w1 /= np.linalg.norm(w1)
    w2 /= np.linalg.norm(w2)
    return	w1, w2


### dict to be saved
data = {'name': [], 'resname': 'ALA', 'resid': [], 'x': [], 'y': [], 'z':[], 'chain': 'A'}


### CA atoms
x, y, z = circularhelix(R, a, w, t)
for i in range(N):
    data['resid'].append(i+1)
    data['name'].append('CA')
    data['x'].append(x[i])
    data['y'].append(y[i])
    data['z'].append(z[i])


### N atoms
xN, yN, zN = circularhelix(R, a, w, t-dt/3)
for i in range(N):
    data['resid'].append(i+1)
    data['name'].append('N')
    data['x'].append(xN[i])
    data['y'].append(yN[i])
    data['z'].append(zN[i])


### C atoms
xC, yC, zC = circularhelix(R, a, w, t+dt/3)
for i in range(N):
    data['resid'].append(i+1)
    data['name'].append('C')
    data['x'].append(xC[i])
    data['y'].append(yC[i])
    data['z'].append(zC[i])


### O and HN atoms
dx = list(xN[4:] - xC[:-4])
dy = list(yN[4:] - yC[:-4])
dz = list(zN[4:] - zC[:-4])

for i in range(4):
    j = i - 4
    dx.append(xN[i] - xC[j])
    dy.append(yN[i] - yC[j])
    dz.append(zN[i] - zC[j])

for i in range(N):
    data['resid'].append(i+1)
    data['name'].append('O')
    norm = np.sqrt(dx[i]**2 + dy[i]**2 + dz[i]**2)
    data['x'].append(xC[i] + 1.35 * dx[i] / norm)
    data['y'].append(yC[i] + 1.35 * dy[i] / norm)
    data['z'].append(zC[i] + 1.35 * dz[i] / norm)

for i in range(N):
    data['resid'].append((i+4) % N + 1)
    data['name'].append('HN')
    norm = np.sqrt(dx[i]**2 + dy[i]**2 + dz[i]**2)
    data['x'].append(xN[ (i+4) % N] - 1.10 * dx[i] / norm)
    data['y'].append(yN[ (i+4) % N] - 1.10 * dy[i] / norm)
    data['z'].append(zN[ (i+4) % N] - 1.10 * dz[i] / norm)


#### CB and HA
#v1 = np.array([x,y,z]).T - np.array([xN,yN,zN]).T
#v1 /= np.linalg.norm(v1, axis=-1)[:,None]
#
#v2 = np.array([xC,yC,zC]).T - np.array([x,y,z]).T
#v2 /= np.linalg.norm(v2, axis=-1)[:,None]
#
#for i in range(N):
#    data['resid'].append(i+1)
#    data['name'].append('CB')
#    w1, w2 = completeTetra(v1[i], v2[i])
#    print(np.linalg.norm(w1), np.linalg.norm(w2))
#
#    data['x'].append(x[i] + 1.8 * w1[0])
#    data['y'].append(y[i] + 1.8 * w1[1])
#    data['z'].append(z[i] + 1.8 * w1[2])
#
#    data['resid'].append(i+1)
#    data['name'].append('HA')
#    data['x'].append(x[i] + 1.35 * w2[0])
#    data['y'].append(y[i] + 1.35 * w2[1])
#    data['z'].append(z[i] + 1.35 * w2[2])

for i in range(N):
    pos = np.repeat([[x[i], y[i], z[i]]], 5, axis=0) + np.random.rand(5, 3) - 0.5
    for j, name in enumerate(['HA', 'CB', 'HB1', 'HB2', 'HB3']):
        data['resid'].append(i+1)
        data['name'].append(name)
        data['x'].append(pos[j][0])
        data['y'].append(pos[j][1])
        data['z'].append(pos[j][2])


### Save
u = mstool.Universe(data=data)
u.sort()
u.write(ifile)


### Append bond
index1 = u.atoms[(u.atoms.name == 'N') & (u.atoms.resid == 1)].index[0] + 1
index2 = u.atoms[(u.atoms.name == 'C') & (u.atoms.resid == N)].index[0] + 1
pdb = open(ifile, 'a')
pdb.write(f"CONECT{index1:5d}{index2:5d}")
pdb.close()


### REM
mstool.REM(ifile, 'final.pdb', pbc=False)
mstool.CheckStructure('final.pdb')

import msprot
import numpy as np

N = 0; dr = 10; Nx=8
x = []; y = []; z = [];
resid = [];
name  = [];

for i in range(0, Nx):
    for j in range(0, Nx):
        for k in range(0, Nx):
            N += 1
            resid.extend([N] * 3)
            name.extend(['PO4', 'GL1', 'C1A'])
            for l in range(3):
                x.append(i * dr + l*2.5)
                y.append(j * dr + l*2.5)
                z.append(k * dr + l*2.5)


u = msprot.Universe(data = {'x': x, 'y': y, 'z': z, 
                            'resid': resid,
                            'resname': 'TOY',
                            'name': name,
                            'chain': 'X'})

u.dimensions = [Nx * dr] * 3 + [90] * 3
u.write('input.pdb')


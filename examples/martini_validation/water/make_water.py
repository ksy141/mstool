import msprot

N = 1; dr = 5; Nx=8
x = []; y = []; z = [];
for i in range(0, Nx):
    for j in range(0, Nx):
        for k in range(0, Nx):
            x.append(i * dr)
            y.append(j * dr)
            z.append(k * dr)

u = msprot.Universe(data = {'x': x, 'y': y, 'z': z, 
                            'resid': [i for i in range(1, Nx**3+1)],
                            'resname': 'W',
                            'name': 'W',
                            'anum': -1,
                            'chain': 'X'})

u.dimensions = [Nx * dr] * 3 + [90] * 3
u.write('W.pdb')


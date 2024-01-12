
W = open('nonbond_params_LJ', 'w')

with open('nonbond_params') as fp:
    for line in fp:
        line = line.split(';')[0]
        line = line.strip()
        if not line: continue

        sl = line.split()
        atom1 = sl[0]
        atom2 = sl[1]
        C6    = float(sl[3])
        C12   = float(sl[4])
    
        sigma6 = C12 / C6
        sigma  = sigma6 ** (1/6)
        epsilon = C6 / 4 / sigma6

        W.write(f'{atom1:10s} {atom2:10s} {1:5d} {sigma:10.5f} {epsilon:10.5f}\n')
W.close()


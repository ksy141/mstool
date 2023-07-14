
def gmx_itp(u, itp, segname=None, nrexcl=1):

    if segname:
        segname = segname.upper()
    else:
        segname = u.atoms['resname'][0]

    itpfile = open(itp + segname + '.itp', 'w')
    itpfile.write('[ moleculetype ]\n')
    itpfile.write('; name  nrexcl\n')
    itpfile.write(f'{segname:s} {nrexcl}\n\n')

    itpfile.write('[ atoms ]\n')
    itpfile.write('; id    type    resnr   residu  atom    cgnr    charge\n')

    for index, atom in u.atoms.iterrows():
        aid   = index + 1
        atype = atom['type']
        aname = atom['name']
        resn  = atom['resname']
        q     = atom['charge']
        itpfile.write(f'{aid:10d} {atype:10s} {1:10d} {resn:10s} {aname:10s} {aid:10d} {q:10.3f}\n')




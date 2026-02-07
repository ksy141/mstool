import glob
import subprocess
import os
import mstool

ifiles = glob.glob('AA_*.pdb')
for ifile in ifiles:
    name = os.path.basename(ifile).split('_')[1].split('.dms')[0]
    if len(name) < 2:
        continue
    u = mstool.Universe(ifile)
    u.atoms = u.atoms[~u.select(':ACE,NMA,NME', returnbA=True)]
    u.write(name + '.pdb')
    subprocess.run(f"martinize2 -f {name}.pdb -dssp -o PROA.itp -nt -maxwarn 10 -x {name}_CG.pdb", shell=True)
    os.rename("molecule_0.itp", f"{name}.itp")
    resname = 'Z' + name
    change = False
    with open(f"{name}.itp") as fin, open(f"{name}_final.itp", "w") as fout:
        for line in fin:
        
            if line.startswith('['):
                read = line.split('[')[1].split(']')[0].strip()
                if read == 'atoms':
                    change = True
                else:
                    change = False
                fout.write(line)
                continue
            
            if change:
                sl = line.split()
                if len(sl) > 3:
                    print(sl)
                    sl[4] = sl[4] + str(int(sl[2].strip())) + ' '
                    line = ' '.join(sl)
                    line += '\n'

            line = line.replace("molecule_0", resname)
            for three in mstool.three2one.keys():
                line = line.replace(three, resname)
            fout.write(line)

    os.remove("PROA.itp")
    u = mstool.Universe(f'{name}_CG.pdb')
    u.atoms.chain = 'X'
    u.atoms.resname = resname
    u.atoms.resid = 1
    u.write(f'{resname}.pdb')
    print(name)

subprocess.run('cat *_final.itp > dipeptide.itp', shell=True)


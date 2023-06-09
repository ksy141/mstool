import mstool
from   mstool.membrane.sphere import Sphere

##### Read Martini FF
#martini = mstool.ReadMartini(define={'FLEXIBLE': 'True'})
#
#### Make a spherical bilayer of a radius of 60 A
#upper = Sphere(r=60, dr=+15, inverse=+1, chain='U', monolayer={'POPC':550, 'DOPC': 550}, martini=martini)
#lower = Sphere(r=60, dr=-15, inverse=-1, chain='L', monolayer={'POPC':200, 'DOPC': 200}, martini=martini)
#u = mstool.Merge(upper.universe.atoms, lower.universe.atoms)
#u.write('sphere.pdb')

### Solvate
mstool.SolvateMartini(structure='sphere.pdb', out='sphere.sol.pdb', t=15, conc=0.15)

#### Martinize
#mstool.MartinizeDMS(dms='sphere.sol.pdb', output='sphere.martini.dms', martini=martini)
#mstool.runMartiniEM( dms_in='sphere.martini.dms', out='sphere.EM.pdb',  soft=True)
#mstool.runMartiniNPT(dms_in='sphere.martini.dms', out='sphere.NPT.pdb', pos_in='sphere.EM.pdb', nsteps=10000)


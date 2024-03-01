import subprocess
import os

def make_mdp(epsilon_r=15.0, define=None):
    with open('input.mdp', 'w') as W:
        if define:
            W.write(f'define                   = {define}\n')

        W.write(f'''
;
; STANDARD MD INPUT OPTIONS FOR MARTINI 2.x
; Updated 15 Jul 2015 by DdJ
;
; for use with GROMACS 5
; For a thorough comparison of different mdp options in combination with the Martini force field, see:
; D.H. de Jong et al., Martini straight: boosting performance using a shorter cutoff and GPUs, submitted.

title                    = Martini

; TIMESTEP IN MARTINI 
; Most simulations are numerically stable with dt=40 fs, 
; however better energy conservation is achieved using a 
; 20-30 fs timestep. 
; Time steps smaller than 20 fs are not required unless specifically stated in the itp file.

integrator               = md
dt                       = 0.002
nsteps                   = 0
nstcomm                  = 100

nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
nstlog                   = 1000
nstenergy                = 100
nstxout-compressed       = 1000
compressed-x-precision   = 100
compressed-x-grps        = 

; NEIGHBOURLIST and MARTINI 
; To achieve faster simulations in combination with the Verlet-neighborlist
; scheme, Martini can be simulated with a straight cutoff. In order to 
; do so, the cutoff distance is reduced 1.1 nm. 
; Neighborlist length should be optimized depending on your hardware setup:
; updating ever 20 steps should be fine for classic systems, while updating
; every 30-40 steps might be better for GPU based systems.
; The Verlet neighborlist scheme will automatically choose a proper neighborlist
; length, based on a energy drift tolerance.
;
; Coulomb interactions can alternatively be treated using a reaction-field,
; giving slightly better properties.
; Please realize that electrostVatic interactions in the Martini model are 
; not considered to be very accurate to begin with, especially as the 
; screening in the system is set to be uniform across the system with 
; a screening constant of 15. When using PME, please make sure your 
; system properties are still reasonable.
;
; With the polarizable water model, the relative electrostatic screening 
; (epsilon_r) should have a value of 2.5, representative of a low-dielectric
; apolar solvent. The polarizable water itself will perform the explicit screening
; in aqueous environment. In this case, the use of PME is more realistic.


cutoff-scheme            = Verlet
nstlist                  = 20
ns_type                  = grid
pbc                      = xyz
verlet-buffer-tolerance  = 0.005

coulombtype              = Cut-off
coulomb-modifier         = None
rcoulomb                 = 1.1
epsilon_r                = {epsilon_r}
vdwtype                  = Cut-off
vdw-modifier             = None
rvdw                     = 1.1

; MARTINI and TEMPERATURE/PRESSURE
; normal temperature and pressure coupling schemes can be used. 
; It is recommended to couple individual groups in your system separately.
; Good temperature control can be achieved with the velocity rescale (V-rescale)
; thermostat using a coupling constant of the order of 1 ps. Even better 
; temperature control can be achieved by reducing the temperature coupling 
; constant to 0.1 ps, although with such tight coupling (approaching 
; the time step) one can no longer speak of a weak-coupling scheme.
; We therefore recommend a coupling time constant of at least 0.5 ps.
; The Berendsen thermostat is less suited since it does not give
; a well described thermodynamic ensemble.
; 
; Pressure can be controlled with the Parrinello-Rahman barostat, 
; with a coupling constant in the range 4-8 ps and typical compressibility 
; in the order of 10e-4 - 10e-5 bar-1. Note that, for equilibration purposes, 
; the Berendsen barostat probably gives better results, as the Parrinello-
; Rahman is prone to oscillating behaviour. For bilayer systems the pressure 
; coupling should be done semiisotropic.

; tcoupl                   = v-rescale 
; tc-grps                  = system
; tau_t                    = 1.0
; ref_t                    = 310
; Pcoupl                   = berendsen
; Pcoupltype               = isotropic
; tau_p                    = 12.0
; compressibility          = 3e-4
; ref_p                    = 1.0

gen_vel                  = yes
gen_temp                 = 310
gen_seed                 = -1

; MARTINI and CONSTRAINTS 
; for ring systems and stiff bonds constraints are defined
; which are best handled using Lincs. 

constraints              = none 
constraint_algorithm     = Lincs
''')


def getGMXEnergy(c, p='topol.top', add_mdrun='', add_grompp='', epsilon_r=15.0, define=None):
    """If fails, use add_mdurn='-rcon 1' or add_mdrun='-ntomp 1' (depending on versions)"""
    make_mdp(epsilon_r, define)

    subprocess.run(f'gmx grompp -f input.mdp -c {c} -o input.tpr -p {p} {add_grompp} -maxwarn -1'.split(), 
        stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT)

    subprocess.run(f'gmx mdrun -deffnm input {add_mdrun} -v'.split(),
        stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT)

    subprocess.run(f'echo Potential | gmx energy -f input.edr',
        shell=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT)

    output = subprocess.run("""cat energy.xvg | tail -n 1 | awk '{print $2}'""", shell=True, stdout=subprocess.PIPE)
    PE = float(output.stdout.decode('utf-8').strip())
    subprocess.run('rm -rf energy.xvg input.mdp input.edr input.gro input.log input.tpr input.xtc mdout.mdp'.split())
    print('{:7s}: {:10.3f} kJ/mol'.format('GMX', PE))
    return PE

def change_dimensions(c, dimensions=[10000.0]*3 + [90.0]*3):
    import MDAnalysis as mda
    u = mda.Universe(c)
    u.dimensions = dimensions
    u.atoms.write(c)
    

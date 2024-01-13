from .version import __version__
from .authors import __authors__

from .core.universe         import Universe, Merge, RemoveOverlappedResidues
from .core.rem              import REM
from .core.readmappings     import ReadMappings
from .core.readmartini      import ReadMartini
from .core.readxml          import ReadXML
from .core.ungroup          import Ungroup
from .core.map              import Map
from .core.gmx2dms          import GMX2DMS
from .core.dms2openmm       import DMS2openmm
from .core.martinizedms     import MartinizeDMS
from .core.loopmodeler      import LoopModeler
from .core.truncateprotein  import TruncateProtein
from .core.solvate_martini  import *
from .core.checkstructure   import CheckStructure
from .core.checktetrahedron import CheckTetrahedron
from .core.seq              import Seq
from .core.dmsfile          import DMSFile
from .core.backmap          import Backmap
from .utils.protein_sel     import three2one, one2three
from .utils.dump            import dumpsql, dumpdf
from .utils.gmx_energy      import getGMXEnergy
from .utils.openmmutils     import *
from .utils.amberselection  import amberSelection
from .utils.add_termini_atoms import addTerminiAtoms
from .utils.cap_termini_residues import capTerminiResidues
from .utils.change_his import changeHIS
from .utils.performance import Performance

from .membrane.bilayerbuilder import BilayerBuilder


# examples
import os
pwd = os.path.dirname(os.path.realpath(__file__))
ONECGBEADSTRUCTURE = pwd + '/examples/Backmapping/Example1_methane/cg.pdb'
POPCSTRUCTURE      = pwd + '/examples/Backmapping/Example4_POPC/cg.pdb'
MULTISTRUCTURE     = pwd + '/examples/Backmapping/Example5_Sphere/cg.pdb'
TRIOSTRUCTURE      = pwd + '/examples/Backmapping/Example6_TRIO/cg.pdb'
TRIOMAPPING        = pwd + '/examples/Backmapping/Example6_TRIO/mapping.dat'
TRIOFF             = pwd + '/examples/Backmapping/Example6_TRIO/ff.xml'
MPCG               = pwd + '/examples/Backmapping/Example7_MembraneProtein/cg.pdb'
MPAA               = pwd + '/examples/Backmapping/Example7_MembraneProtein/protein_AA.pdb'


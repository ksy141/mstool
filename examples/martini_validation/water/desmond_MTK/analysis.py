from usemda import *
import numpy as np
u = jobid2mda('output.dtr', system='input.dms', view_frames='*')
toxtc(u, 'trj')

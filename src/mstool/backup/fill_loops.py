from .utils.universe import Universe
from .utils.seq      import Seq
import numpy as np
import pandas as pd

class FillLoops:
	def __init__(self, mapping, universe, sequence):
		self.universe = universe
		self.mapping  = mapping
		self.sequence = sequence

		# python 3.7+
		self.chains = list(dict.fromkeys(self.universe.atoms.chain))
		assert isinstance(universe, Universe), "universe is not Universe instance"
		assert isinstance(sequence, Seq),      "sequence is not Seq instance"

		seqfromU = Seq(universe=universe)
		seqfromF = sequence

		assert len(seqfromU.seq.keys()) == len(seqfromF.seq.keys()), 'len(chains)'
		assert len(self.chains) == len(seqfromF.seq.keys()), 'len(chains)'

		nchains = len(seqfromF.seq.keys())
		for chainnum in range(nchains):
			s = seqfromU.seq[chainnum]['one']
			missing = np.where(np.array([*s]) == '-')[0]
			grouped = self.consecutive(missing)

			# disregard the missing N-termini residues
			if grouped[0][0] == 0 and len(grouped) != 1:
				grouped = grouped[1:]

			elif grouped[0][0] == 0 and len(grouped) == 1:
				assert 0 != 1, 'there is nothing to fill except for the N-termini'

			for group in grouped:
				self.fill(chainnum, group)

		self.universe.sort()


	def fill(self, chainnum, group):
		dg = {'name': [], 'resid': [], 'resname': [],'chain':[],
              'x': [],'y': [],'z': [], 'bfactor': 0.0}

		N = len(group)
		chain = self.chains[chainnum]

		bA1 = self.universe.atoms.chain  == chain
		bA2 = (self.universe.atoms.name  == 'CA') | (self.universe.atoms.name == 'BB')
		bA3 = self.universe.atoms.resid == group[0] - 1 + 1
		bA4 = self.universe.atoms.resid == group[-1] + 1 + 1

		start = self.universe.atoms[bA1 & bA2 & bA3][['x', 'y', 'z']].to_numpy()[0]
		end   = self.universe.atoms[bA1 & bA2 & bA4][['x', 'y', 'z']].to_numpy()[0]
		dr    = end - start
		step  = np.linspace(0, 1, N+2)

		for i in range(len(group)):
			resid = group[i] + 1
			posv  = start + step[i+1] * dr
			resname = self.sequence.resid(resid=resid, chainnum=chainnum)[1]
			names = self.mapping.RESI[resname]['CGAtoms'].keys()

			for name in names:
				dg['name'].append(name)
				dg['resid'].append(resid)
				dg['resname'].append(resname)
				dg['chain'].append(chain)
				pos = posv + np.random.rand(3) - 0.5
				dg['x'].append(pos[0])
				dg['y'].append(pos[1])
				dg['z'].append(pos[2])

		newu = Universe(data=dg)
		self.universe.atoms = pd.concat([self.universe.atoms, newu.atoms])



	def consecutive(self, data, stepsize=1):
	    return np.split(data, np.where(np.diff(data) != stepsize)[0]+1)
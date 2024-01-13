from  ..utils.protein_sel import one2three, three2one
from  .universe           import Universe
import numpy as np

class Seq:
	def __init__(self, fasta=None, structure=None):
		"""Get amino acid sequence from fasta or structure.
		
		Parameters
		----------
		fasta : str
			Filename of a fasta file
		structure : str
			Filename of a structure file

		Attributes
		----------
		seq : dict
			Contains all the sequence information.
			``seq[0]['three']`` is an array of three letter code of the first chain.
			``seq[0]['one']`` is an array of one letter code of the first chain.
			--- or - represents a missing residue.
	

		Examples
		--------
		>>> s = mstool.Seq(fasta='seq.fasta')
		>>> s.seq[0]['one'] # one letter code of the first chain
		>>> s.seq[1]['three'] # three letter code of the second chain
		>>> s.pretty(s.seq[0]['one'], split=40) # pretty print of one letter code of the first chain
		"""
		
		self.seq = {}

		### Read FASTA
		if fasta:
			chainnum = -1
			with open(fasta) as W:
				for line in W.readlines():
					line = line.strip()
					if not line: continue

					if line.startswith('>'):
						chainnum += 1
						self.seq[chainnum] = {'one': ''}
						continue

					self.seq[chainnum]['one'] += line


			### Assign resid and three-letter amino residues
			for key, value in self.seq.items():
				value['resid'] = [i for i in range(1, 1+len(value['one']))]
				value['three'] = []

				for one in value['one']:
					value['three'].append(one2three[one])

		### Read UNIVERSE
		if structure:
			u = Universe(structure)

			# python 3.7+
			chains = list(dict.fromkeys(u.atoms.chain))

			for chainnum, chain in enumerate(chains):
				self.seq[chainnum] = {'one': '',   'three': [], 
				                      'resid': [], 'chain': chain}

				# you only select resnames shown in three2one (protein)
				df = u.atoms.loc[(u.atoms.chain == chain) &
					(u.atoms.resname.isin(three2one.keys()))]

				if len(df) == 0:
					continue

				for i in range(1, 1+max(df.resid)):
					dg = df.loc[df.resid == i]
					self.seq[chainnum]['resid'].append(i)

					if len(dg) != 0:
						assert len(set(dg.resname)) == 1, 'chain-resid-resname dup?'
						three = dg.resname.iloc[0]
						self.seq[chainnum]['one'] += three2one[three]
						self.seq[chainnum]['three'].append(three)

					else:
						self.seq[chainnum]['one'] += '-'
						self.seq[chainnum]['three'].append('---')

		### numpy
		for key in self.seq.keys():
			self.seq[key]['three'] = np.array(self.seq[key]['three'])
			self.seq[key]['resid'] = np.array(self.seq[key]['resid'])



	def __repr__(self):
		ret = ''
		for key, value in self.seq.items():
			ret += '{:d} {:s}\n'.format(key, value['one'])
		return ret

	def resid(self, resid, chainnum=0):
		'''Get one and three letter codes of ``resid`` of ``chainnum``.
		
		Parameters
		----------
		resid : int
		chainnum : int

		Examples
		--------
		>>> s = mstool.Seq('seq.fasta')
		>>> s.resid(resid=39, chainnum=0)
		'''

		n_res = len(self.seq[chainnum]['three'])
		assert resid - 1 > -0.5,  'resid should be equal or higher than 1.'
		assert resid - 1 < n_res, f'resid should be lower than {n_res}'
		return self.seq[chainnum]['one'][resid-1], self.seq[chainnum]['three'][resid-1]

	def sequence(self, chainnum=0, pretty=True):
		if pretty:
			return self.pretty(self.seq[chainnum]['one'])
		else:
			return self.seq[chainnum]['one']

	def pretty(self, seq, split=40):
		ret = ''
		for i, s in enumerate(seq, 1):
			if i % split == 1:
				ret += '\n'
				ret += f'{i:4d} '
			if i % (split // 4) == 1:
				ret += ' '
			ret += s

		print(ret)
		#return ret





#!/usr/bin/env python2
'''
Which TMS is so and so HET binding? Let's find out!
'''
from __future__ import print_function, division

import subprocess, re, os
import numpy as np
import xml.etree.ElementTree as ET
import argparse
import Bio.PDB

MAX_INTERACTION_DISTANCE = 7

try: CODE = Bio.PDB.protein_letters_3to1
except AttributeError: CODE = Bio.PDB.to_one_letter_code

def dict2formula(atoms):
	def subscript(x):
		if x <= 1: return ''
		else: return str(x)
	out = ''

	for e in ('C', 'H', 'N', 'O'):
		if e in atoms: out += e + subscript(atoms.pop(e))
	for e in sorted(atoms):
		out += e + subscript(atoms.pop(e))
	return out

class Structure:
	def __init__(self, fn):
		self.fn = fn
		if '.cif' in fn.lower(): parser = Bio.PDB.MMCIFParser()
		else: parser = Bio.PDB.PDBParser()

		self.structure = parser.get_structure(fn, fn)
		self.hetres = []
		self.formulas = {}

	def get_het(self, water=False, modres=False):
		for residue in self.structure.get_residues():
			if not modres:
				if residue.get_resname() in CODE: continue
			if residue.get_id()[0].startswith('H_'): self.hetres.append(residue)
			if water and residue.get_id()[0] == 'W': self.hetres.append(residue)
		self.formulas = {}
		for res in self.hetres:
			atoms = {}
			for atom in res:
				elem = atom.get_fullname()[:2].strip().title()
				try: atoms[elem] += 1
				except KeyError: atoms[elem] = 1
			self.formulas[res] = dict2formula(atoms)
		return self.hetres

	def get_seq(self, outfmt='fasta', wrap=80):
		sequences = {}

		out = ''
		for model in self.structure:
			for chain in model:
				#sequences[chain] = {}
				out += '>%s_%s\n' % (self.fn, chain.get_id())
				n = 1
				for residue in chain:
					resn = residue.get_resname()
					#try: sequences[chain][residue.get_id()[1]] = CODE[resn]
					#except KeyError: pass
					try: out += CODE[resn]
					except KeyError: continue
					if wrap and not (n % wrap): out += '\n'
					n += 1
				out += '\n\n'
		return out
class SubstrateLocator:
	def __init__(self, pdbfn, water=False, modres=False, db='pdbtm'):
		self.structure = Structure(pdbfn)
		self.structure.get_het(water=water, modres=modres)
		self.pdbtm = PDBTM(os.path.basename(pdbfn)[:4], db=db)

	def get_het_contacts(self, distance=4, cathreshold=10):
		'''
distance: Maximum minimum distance allowed for contacts
cathreshold: Maximum minimum distance to C-alpha considered for plausible contact. This helps reduce the computational load involved
		'''
		hetatms = []
		contacts = {}
		for residue in self.structure.hetres:
			contacts = []

			for atom in residue: 
				hetatms.append(atom)

		for model in self.structure.structure:
			for chain in model:
				for residue in chain:
					if residue.get_resname() not in CODE: continue

					mindist = None
					for atom in residue:
						if atom.name == 'CA':
							for hetatm in hetatms:
								dist = np.linalg.norm(hetatm.get_coord()-atom.get_coord())
								if mindist is None or dist < mindist: mindist = dist
								if dist < cathreshold: break

					if mindist > cathreshold: continue

					#if mindist <= distance: 
					#	contacts[hetatm.get_parent()].add(residue)
					#	#print('contact found between %s and %s: %s-%s (%0.2f Angstroms)' % (atom.get_parent().get_full_id(), hetatm.get_parent().get_full_id(), atom.name, atom.name, dist))
					#	pass

					for atom in residue:
						found = 0
						for hetatm in hetatms:
							dist = np.linalg.norm(hetatm.get_coord()-atom.get_coord())
							if dist < distance:
								tms = self.pdbtm.get_tms(atom)

								contacts.append(Connection(atom, hetatm, tms))
								#print('%s%d and %s%d: %s-%s (%0.2f Angstroms)' % (\
								#	atom.get_parent().get_resname(), \
								#	atom.get_parent().get_id()[1], \
								#	hetatm.get_parent().get_resname(), \
								#	hetatm.get_parent().get_id()[1], \
								#	atom.name, \
								#	hetatm.name, \
								#	dist))
								found = 1
								break
						if found: break
			return contacts


	def get_tmss(self):
		self.structure.get_seq()

class Connection(object):
	'''
Container for atom-atom distances
	'''
	def __init__(self, protatom, ligatom, tms):
		self.protatom = protatom
		self.ligatom = ligatom
		self.tms = tms

	def get_length(self): return np.linalg.norm(self.ligatom.get_coord() - self.protatom.get_coord())

	def __lt__(self, other):
		if len(self) < len(other): return True
		else: return False

	def to_tabular(self):
		out = ''
		out += '\t%s' % self.ligatom.get_parent().get_resname()
		out += '\t%s' % self.ligatom.get_parent().get_id()[1]
		out += '\t%s' % self.protatom.get_parent().get_resname()
		out += '\t%s' % self.protatom.get_parent().get_parent().id
		out += '\t%s' % self.protatom.get_parent().get_id()[1]
		out += '\t%s-%s' % (self.ligatom.get_fullname()[:2].title().strip(), self.protatom.get_fullname()[:2].title().strip())
		out += '\t%s-%s' % (self.ligatom.get_fullname(), self.protatom.get_fullname())
		out += '\t%s' % self.tms
		out += '\t%0.2f' % self.get_length()
		return out.strip()

	def __str__(self): return self.to_tabular()

#def print_contactlist(l, pdb, db='pdbtm', distance=-1):
	##for hetres in sorted(l):
	##	print('### POSSIBLE INTERACTIONS FOR LIGAND %s %d ###' % (hetres.get_resname(), hetres.get_id()[1]))
	##	print('#Amino acid\tAA index\tligand\tligand id\tatom-atom\tdistance')
	##	for contact in l[hetres]: print(contact)
	##print('#Amino acid\tAA index\tligand\tligand id\tatom-atom\tdistance')
	#pdbtm = PDBTM(pdb, db)
	#print('######## %s: residues within %0.2f Angstroms of ligands ########' % (pdb, distance))
	#print('#lign\tligi\tresn\tchain\tresi\tixn\ttms\tdist')
	#for hetatm in l:
	#	for atom in l[hetatm]:
	#		tms = pdbtm.get_tms(atom.get_parent().get_id()[1], chain=atom.get_parent().get_parent().get_id())
	#		if tms == 'geom': continue
	#		line = ('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%0.2f' % (\
	#			hetatm.get_parent().get_resname(), \
	#			hetatm.get_parent().get_id()[1], \
	#			atom.get_parent().get_resname(), \
	#			atom.get_parent().get_parent().id, \
	#			atom.get_parent().get_id()[1], \
	#			'%s-%s' % (hetatm.get_fullname(), atom.get_fullname()), \
	#			tms, \
	#			l[hetatm][atom], \
	#		))
	#		print(line)

class PDBTM:
	def __init__(self, pdb, db='pdbtm'):
		self.regions = {}
		self.delete = []
		tree = ET.parse('%s/%s.xml' % (db, pdb))
		root = tree.getroot()
		for x in root:
			if x.tag.endswith('CHAIN'):
				chainid = x.attrib['CHAINID']
				self.regions[chainid] = []
				if x.attrib['TYPE'] == 'non_tm': continue
				for y in x:
					if y.tag.endswith('REGION'):
						if y.attrib['type'] in ('H', 'C', 'B'):
							self.regions[chainid].append(\
								(\
									int(y.attrib['pdb_beg']),
									int(y.attrib['pdb_end']),
								)\
							)
			elif x.tag.endswith('BIOMATRIX'): 
				for y in x:
					if y.tag.endswith('DELETE'): self.delete.append(y.attrib['CHAINID'])

	#def get_tms(self, i, chain='A'):
	def get_tms(self, query):

		if type(query) is Bio.PDB.Atom.Atom: residue = query.get_parent()
		elif type(query) is Bio.PDB.Residue.Residue: residue = query
		else: print(type(query))

		i = residue.get_id()[1]
		chain = residue.get_parent().id

		try: self.regions[chain]
		except KeyError: 
			if chain in self.delete: return 'geom'
			else: return 'solu'
		n = 1
		for reg in self.regions[chain]:
			if reg[0] <= i <= reg[1]: return str(n)
			elif i < reg[0]: 
				if n == 1: return 'LN-%d' % n
				else: return 'L%d-%d' % (n-1, n)
			n += 1
		return 'L%d-C' % n

def tabulate_contacts(pdbs, max_distance=5, include=[], exclude=[], min_distance=None, subtitle=None):
	out = ''
	for pdb in pdbs:
		sl = SubstrateLocator(pdb)
		contacts = sl.get_het_contacts(distance=max_distance)
		out += '\n#### Ligands of %s ####' % os.path.basename(pdb)

		if subtitle is not None: out += '\n# %s' % subtitle

		out += '\n# Max. distance: %0.2f Angstroms' % max_distance
		if min_distance is None: out += '\n# Min. distance: (none)'
		else: out += '\n# Min. distance: %0.2f Angstroms' % min_distance
		out += '\n# Must include: '

		if include: 
			inclstr = ''
			for x in include: inclstr += '%s, ' % x
			out += inclstr[:-2]
		else: out += '(none)'

		out +='\n# Must exclude: '
		if exclude: 
			exclstr = ''
			for x in exclude: exclstr += '%s, ' % x
			out += exclstr[:-2]
		else: out += '(none)'

		out += '\n#l_atom\
\tl_resi\
\tp_resn\
\tp_chn\
\tp_resi\
\tatoms\
\tnames\
\ttms\
\tdistance'

		n = 0
		for c in contacts: 
			if include and c.ligatom.get_parent().get_resname() not in include: continue
			if exclude and c.ligatom.get_parent().get_resname() in include: continue
			if min_distance is not None and c.get_length() < min_distance: continue
			out += '\n%s' % c
			n += 1
		if n == 0: out += '\n#\n# No contacts found\n#'
		out += '\n'
	return out.strip()

if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument('-d', metavar='distance', type=float, default=5., help='maximum distance')
	parser.add_argument('-i', nargs='+', help='must include these ligands')
	parser.add_argument('-e', nargs='+', help='must exclude these ligands')
	parser.add_argument('--pdbtm', default='pdbtm', help='PDBTM database location')
	parser.add_argument('pdb', nargs='+', help='PDB/mmCIF inputs. CIF files must contain .cif')

	parser.add_argument('-w', action='store_true', help='get TMSs for water contacts')
	parser.add_argument('-m', action='store_true', help='get TMSs for modified residue contacts. Using this argument is recommended only when proteins of interest interact with peptides; most results obtained using this option are likely to be selenomethionines (MSE), which are used to improve resolution and are not generally physiologically relevant')

	args = parser.parse_args()

	print(tabulate_contacts(args.pdb, max_distance=args.d, include=args.i, exclude=args.e))

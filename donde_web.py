#!/usr/bin/env python
from __future__ import print_function
print('Content-type: text/html\r\n\r\n')

import cgi
import cgitb
import os, sys, re

file_dir = os.path.dirname(__file__)
sys.path.append(file_dir)

DB='reduced_donde.tsv'

class Hit:
	def __init__(self, src):
		self.src = src
		self.tcdb = []
		self.maxdist = -1
		self.mindist = -1
		self.include = []
		self.exclude = []
		self.contacts = []
	def __str__(self):
		out = '#### Ligands of %s ####' % self.src

		for t in self.tcdb: out += '\n# %s' % t.strip()

		if self.maxdist == -1 or self.maxdist is None: maxdist = '(none)'
		else: maxdist = '%0.2f Angstroms' % self.maxdist
		out += '\n# Max. distance: %s' % maxdist

		if self.mindist == -1 or self.mindist is None: mindist = '(none)'
		else: mindist = '%0.2f Angstroms' % self.mindist
		out += '\n# Min. distance: %s' % mindist

		if not self.include: include = '(none)'
		else: 
			include = ''
			for i in self.include: include += i + ', '
			include = include[:-2]

		out += '\n# Must include: %s' % include

		if not self.exclude: exclude = '(none)'
		else: 
			exclude = ''
			for i in self.exclude: exclude += i + ', '

		out += '\n# Must exclude: %s' % exclude

		out += '\n#l_atom\tl_resi\tp_resn\tp_chn\tp_resi\tatoms\tnames\ttms\tdistance'

		contact = ''
		for c in self.contacts: 
			contact += '\n'
			for f in c:
				contact += '%s\t' % f
			contact = contact[:-1]
		out += contact
		return out

def get_queries(fams, max_distance=5., min_distance=None, include=None, exclude=None):
	hits = []
	with open('reduced_donde.tsv') as f:
		current = Hit('')
		for l in f:
			if l.startswith('#### Ligands'):
				if current.src and current.contacts: 
					hits.append(current)
					current = Hit(l.split()[3])
					current.include = include
					current.exclude = exclude
					current.maxdist = max_distance
					current.mindist = min_distance
				else: 
					current = Hit(l.split()[3])
					current.include = include
					current.exclude = exclude
					current.maxdist = max_distance
					current.mindist = min_distance
			elif re.match('# [0-9]\.', l): current.tcdb.append(l[2:])
			elif not l.startswith('#'):
				sl = l.strip().split('\t')
				dist = float(sl[8])
				hetres = sl[0]
				if (max_distance is not None) and (dist > max_distance): continue
				if (min_distance is not None) and (dist < min_distance): continue
				if (include is not None) and (hetres not in include): continue
				if (exclude is not None) and (hetres in include): continue
				found = 0
				for fam in fams:
					for tcid in current.tcdb:
						if tcid.startswith(fam):
							found = 1
							break
				if found: current.contacts.append(sl)
	return hits

def main(fams, max_distance=5., min_distance=None, include=None, exclude=None):
	hits = get_queries(fams=fams)
	for h in hits: print(h)
	#for pdb in sorted(pdbs): 

if __name__ == '__main__':

	form = cgi.FieldStorage()

	fams = form.getvalue('families')
	if fams is None: 
		print('[ERROR]: Please specify at least one TCID<br/><a href="javascript:history.back()">Return to form</a><br/>')
		exit()
	elif not fams: 
		print('[ERROR]: Please specify at least one TCID<br/><a href="javascript:history.back()">Return to form</a><br/>')
		exit()
	else: fams = re.split('\s*,\s*', fams)

	maxdist = form.getvalue('maxdistance')
	if maxdist is None: maxdist = 5.0
	elif not maxdist: maxdist = 5.0
	else: 
		try: maxdist = float(maxdist)
		except ValueError: maxdist = None

	mindist = form.getvalue('mindistance')
	if mindist is None: mindist = None
	elif not mindist.strip(): mindist = None
	else:
		try: mindist = float(mindist)
		except ValueError: mindist = None

	include = form.getvalue('include')
	if include is None: include = None
	elif not include.strip(): include = None
	else: include = re.split('\s*,\s*', include)

	exclude = form.getvalue('exclude')
	if exclude is None: exclude = None
	elif not exclude.strip(): exclude = None
	else: exclude = re.split('\s*,\s*', exclude)

	print('<pre>')
	main(fams=fams, max_distance=maxdist, min_distance=mindist, include=include, exclude=exclude)
	print('</pre>')


#!/usr/bin/env python
from __future__ import print_function
print('Content-type: text/html\r\n\r\n')

import cgi
import cgitb
import donde
import os, sys, re

file_dir = os.path.dirname(__file__)
sys.path.append(file_dir)
import Deuterocol1

def get_queries(fams, p1dir='deuterocol1', pdbtmdir='pdbtm'):
	p1 = Deuterocol1.Protocol1(pdbtmdir=pdbtmdir, outdir=p1dir, offline=True)
	p1.blast_pdbs()
	for fam in fams: p1.get_queries(startswith=fam)
	pdbs = {}
	#for hit in p1.hits: 
	for fam in fams: 
		#print(p1.blast.by_target(fam))
		pdbs.update(p1.blast.by_target(fam))

	return pdbs

def main(fams, p1dir='deuterocol1', pdbtmdir='pdbtm', max_distance=5., min_distance=None, include=None, exclude=None):
	pdbs = get_queries(fams=fams, p1dir=p1dir, pdbtmdir=pdbtmdir)
	for pdb in sorted(pdbs):
		fn = ('%s/pdbs_raw/%s.pdb' % (p1dir, pdb[:4]),)

		print(donde.tabulate_contacts(fn, max_distance=max_distance, min_distance=min_distance, include=include, exclude=exclude, subtitle=sorted(pdbs[pdb])[0]))
	#pdbfns = []
	#for pdb in sorted(pdbs): pdbfns.append('%s/pdbs_raw/%s.pdb' % (p1dir, pdb[:4]))
	#print(tabulate_contacts(pdbfns, max_distance=max_distance, min_distance=min_distance, include=include, exclude=exclude))

#main(fams=['1.H.1','8.A.16'])

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

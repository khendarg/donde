#!/usr/bin/env python
print('Content-type: text/html\r\n\r\n')
from __future__ import print_function

import cgi, cgitb


import Deuterocol1, donde

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
	#form = cgi.FieldStorage()
	print('TEST')

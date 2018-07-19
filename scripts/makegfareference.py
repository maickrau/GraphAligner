#!/usr/bin/python

import sys

graphfile = sys.argv[1]
nodeidfile = sys.argv[2]
referencefastafile = sys.argv[3]

with open(graphfile) as gf:
	with open(nodeidfile, 'w') as nf:
		with open(referencefastafile, 'w') as rff:
			rff.write(">ref\n")
			for line in gf:
				parts = line.split('\t')
				if parts[0] != 'S': continue
				number = parts[1]
				nf.write(number + "\n")
				seq = parts[2].rstrip()
				rff.write('NN')
				rff.write(seq)

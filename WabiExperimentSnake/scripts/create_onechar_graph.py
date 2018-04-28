#!/usr/bin/python

import fileinput

next_node = 1
for line in fileinput.input():
	l = line.strip()
	for c in l:
		print("S\t" + str(next_node) + "\t" + c)
		if next_node > 1: print("L\t" + str(next_node-1) + "\t+\t" + str(next_node) + "\t+\t0M")
		next_node += 1

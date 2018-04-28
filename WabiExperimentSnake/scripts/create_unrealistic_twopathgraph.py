#!/usr/bin/python

import fileinput
import random

next_node = 3
# add two fake nodes to simplify the loop
print("S\t1\tA")
print("S\t2\tT")
for line in fileinput.input():
	l = line.strip()
	for c in l:
		for src in [next_node-2, next_node-1]:
			for dst in [next_node, next_node+1]:
				print("L\t" + str(src) + "\t+\t" + str(dst) + "\t+\t0M")
		print("S\t" + str(next_node) + "\t" + c)
		print("S\t" + str(next_node+1) + "\t" + "ATCG"[random.randint(0, 3)])
		next_node += 2

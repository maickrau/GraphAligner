#!/usr/bin/python

import fileinput
import random

snp_rate = 0.1

seq_so_far = ""
last_nodes = []
next_node = 1
for line in fileinput.input():
	l = line.strip()
	last_cut = 0
	for i in range(0, len(l)):
		if random.uniform(0, 1) <= snp_rate:
			for last_node in last_nodes:
				print("L\t" + str(last_node) + "\t+\t" + str(next_node) + "\t+\t0M")
			seq_so_far += l[last_cut:i]
			# add an artificial A because the nodes can't be zero length
			if seq_so_far == "": seq_so_far = "A"
			print("S\t" + str(next_node) + "\t" + seq_so_far)
			print("S\t" + str(next_node+1) + "\t" + l[i])
			print("S\t" + str(next_node+2) + "\t" + "ATCG"[random.randint(0, 3)])
			print("L\t" + str(next_node) + "\t+\t" + str(next_node+1) + "\t+\t0M")
			print("L\t" + str(next_node) + "\t+\t" + str(next_node+2) + "\t+\t0M")
			last_nodes = [next_node+1, next_node+2]
			next_node += 3
			seq_so_far = ""
			last_cut = i
	seq_so_far += l[last_cut+1:]

if len(seq_so_far) > 0:
	for last_node in last_nodes:
		print("L\t" + str(last_node) + "\t+\t" + str(next_node) + "\t+\t0M")
		print("S\t" + str(next_node) + "\t" + seq_so_far)

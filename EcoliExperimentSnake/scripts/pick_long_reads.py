#!/usr/bin/python

import fileinput

read_length_cutoff = 1000

current_name = ""
current_seq = ""
for line in fileinput.input():
	if line[0] == '>':
		if len(current_seq) > read_length_cutoff:
			print(">" + current_name)
			print(current_seq)
		current_name = line[1:].strip()
		current_seq = ""
	else:
		current_seq += line.strip()

if len(current_seq) > read_length_cutoff:
	print(">" + current_name)
	print(current_seq)
current_name = line[1:].strip()
current_seq = ""

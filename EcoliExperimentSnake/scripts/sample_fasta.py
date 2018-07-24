#!/usr/bin/python

import fileinput
import random

sample_rate = 0.01

printing_this = False
for line in fileinput.input():
	if line[0] == '>':
		printing_this = random.uniform(0, 1) < sample_rate
	if printing_this:
		print(line.strip())

#!/usr/bin/python

import sys

input_graph_file = sys.argv[1]
output_graph_file = sys.argv[2]
output_mapping_file = sys.argv[3]

nodenumbers = {}
nextnumber = 1

with open(output_graph_file, 'w') as f:
	with open(input_graph_file) as input_graph:
		for line in input_graph:
			l = line.strip()
			if len(l) == 0: f.write(l + '\n')
			if l[0] == 'S':
				parts = l.split('\t')
				if parts[1] not in nodenumbers:
					nodenumbers[parts[1]] = nextnumber
					nextnumber += 1
				parts[1] = str(nodenumbers[parts[1]])
				f.write('\t'.join(parts) + '\n')
			if l[0] == 'L':
				parts = l.split('\t')
				if parts[1] not in nodenumbers:
					nodenumbers[parts[1]] = nextnumber
					nextnumber += 1
				if parts[3] not in nodenumbers:
					nodenumbers[parts[3]] = nextnumber
					nextnumber += 1
				parts[1] = str(nodenumbers[parts[1]])
				parts[3] = str(nodenumbers[parts[3]])
				f.write('\t'.join(parts) + '\n')

outputmapping = [(node_id, nodenumbers[node_id]) for node_id in nodenumbers]
outputmapping.sort(key= lambda x: x[1])
with open(output_mapping_file, 'w') as f:
	for mapping in outputmapping:
		f.write(mapping[0] + '\n')

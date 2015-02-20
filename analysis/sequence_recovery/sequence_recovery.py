# example usage:
# python sequence_recovery.py ../../output/fixbb_sequences/ ../../output/backrub_sequences/ > sequence_recovery.txt

import sys
import os

native_sequences = {}

for file in os.listdir("native_sequences"):
	if ".fasta.txt" in file:
		f = open("native_sequences/"+file, 'r')
		seq = ''
		for line in f:
			if line[0] != ">":
				seq = line.strip()
		native_sequences[file.split('.')[0]] = seq

structs = {}
f = open('../domain_sequences.txt', 'r')
for line in f:
	line = line.strip()
	structs[line.split()[0]] = line.split()[1]

input_dirs = sys.argv[1:]
dir_values = {}

for input_dir in input_dirs:
	values = []
	for file in os.listdir(input_dir):
		f = open(input_dir+"/"+file)
		native = native_sequences[structs[file.split('_')[0]]]
		sum = 0
		count = 0
		for line in f:
			if line[0] != ">":
				match = 0
				total = 0
				line = line.strip()
				for i in range(0, len(line.strip())):
					if line[i] == native[i]:
						match += 1
					total += 1
				sum += float(match) / float(total)
				count += 1
		values.append(float(sum) / float(count))
	dir_values[input_dir] = values

length = len(dir_values[input_dirs[0]])
for k in range(0,length):
	for input_dir in input_dirs:
		print dir_values[input_dir][k],
	print


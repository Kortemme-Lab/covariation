# The MIT License (MIT)
#
# Copyright (c) 2015 Noah Ollikainen
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

"""
example usage:
python sequence_recovery.py ../../output/fixbb_sequences/ ../../output/backrub_sequences/ > sequence_recovery.txt
"""


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


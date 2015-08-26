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

# Modifications by Shane O'Connor

import os
import sys

dir1 = "backrub_entropy"
dir2 = "native_entropy"
temps = ["0.0","0.3","0.6","0.9","1.2","1.8","2.4"]
domain_positions = {}
natural_entropy = {}
design_entropy = {}
valid_positions = set()
for file in os.listdir(dir2):
    domain = file.split('.')[0]
    f = open(dir2+'/'+file)
    positions = set()
    for line in f:
        positions.add(line.split()[0])
        natural_entropy[domain+"_"+line.split()[0]] = line.split()[1]
    domain_positions[domain] = positions
for kt in temps:
    design_entropy[kt] = {}
    for file in os.listdir(dir1):
        if kt in file:
            #if kt == "0.9" and "relax" not in file:
            domain = file.split('_')[0]
            f = open(dir1+'/'+file, 'r')
            for line in f:
                if line.split()[0] in domain_positions[domain]:
                    design_entropy[kt][domain+"_"+line.split()[0]] = line.split()[1]
                    #if we add valid_positions here, there should be 3009 total
                    #valid_positions.add(domain+"_"+line.split()[0])

#if we add valid_positions here, there should be 2778 total 
f = open('valid_positions', 'r')
for line in f:
    valid_positions.add(line.strip())

pos_keys = []
for p in valid_positions:
    ts = p.split('_')
    pos_keys.append((ts[0], int(ts[1]), p))

print('\t'.join(['Domain', 'Position', 'Native', 'Fixed'] + map(str, temps[1:])))
for p in sorted(pos_keys):
    s = [p[0], p[1]]
    s.append(natural_entropy[p[2]])
    for kt in temps:
        s.append(design_entropy[kt][p[2]])
    print('\t'.join(map(str, s)))

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
python profile_similarity.py ../../output/fixed_backbone/ ../../output/backrub/ > profile_similarity.txt
"""


import os
import sys
import math
import operator
import numpy
from scipy.stats import hypergeom

def get_natural_sequences(fastafile):

    f = open(fastafile, 'r')
    gapped_sequences = []
    gapped_length = 0
    length = 0
    aa = 'ACDEFGHIKLMNPQRSTVWY'
    name = ''
    seq = ''
    test = []
    for line in f:
        if line[0] == ">":
            if name != '':
                gapped_sequences.append(seq)
                #print seq[:10]
            name = line[1:].strip()
            seq = ''
        else:
            seq += line.strip()
    gapped_sequences.append(seq)
    gapped_length = len(seq)

    return gapped_sequences


def get_rosetta_sequence(fastafile, indicefile):

    f = open(indicefile, 'r')
    pfam_indices = {}
    for line in f:
        pfam_indices[int(line.split()[2])-1] = line.split()[0]
    f.close()

    f = open(fastafile, 'r')
    gapped_sequences = set()
    gapped_length = 0
    length = 0
    aa = 'ACDEFGHIKLMNPQRSTVWY'
    name = ''
    seq = ''
    test = []
    for line in f:
        if line[0] == ">":
            if name != '':
                gapped_sequences.add(seq)
                #print seq[:10]
            name = line[1:].strip()
            seq = ''
        else:
            seq += line.strip()
    gapped_sequences.add(seq)
    gapped_length = len(seq)

    sequences = []
    ungapped = []
    to_remove = set()
    for seq in gapped_sequences:
        if len(seq) != gapped_length:
            to_remove.add(seq)
    for seq in to_remove:
        gapped_sequences.remove(seq)

    for i in range(0, gapped_length):
        counter = 0
        for seq in gapped_sequences:
            char = seq[i]
            if char not in aa:
                counter += 1
        if counter < 1:
            ungapped.append(i)

    length = len(ungapped)

    ungapped_indices = {}
    designed_indices = {}
    for seq in gapped_sequences:
        ungapped_seq = ''
        count = 0
        for i in ungapped:
            if i in pfam_indices:
                ungapped_seq += seq[i]
                ungapped_indices[count] = pfam_indices[i]
                designed_indices[int(pfam_indices[i])] = count
                count += 1
        sequences.append(ungapped_seq)
        length = len(ungapped_seq)

    return sequences, designed_indices


def get_covarying_pairs(nature_mi, design_mi):
    pairs1 = {}
    pairs2 = {}
    f = open(nature_mi)
    threshold1 = 0
    to_sort = []
    dists = {}
    natural_scores = []
    for line in f:
        i = line.split()[3]
        j = line.split()[4]
        z = float(line.split()[7])
        dist = float(line.split()[8])
        if z > 0:
            z = math.sqrt(z)
        else:
            z = math.sqrt(z*-1)*-1
        #if dist <= 6:
        dists[i+" "+j] = dist
        pairs1[i+" "+j] = z
        natural_scores.append(z)
    f = open(design_mi, 'r')
    threshold2 = 0
    designed_positions = set()
    for line in f:
        i = line.split()[0]
        j = line.split()[1]
        z = float(line.split()[4])
        designed_positions.add(i)
        designed_positions.add(j)
        if z > 0:
            z = math.sqrt(z)
        else:
            z = math.sqrt(z*-1)*-1
        pairs2[i+" "+j] = z
    a_count = 0
    b_count = 0
    overlap = 0
    total = 0
    success = set()
    failure = set()
    to_sort1 = []
    to_sort2 = []
    positions = set()
    for p in pairs1:
        if p in pairs2:
            i = p.split()[0]
            j = p.split()[1]
            to_sort1.append(pairs1[p])
            to_sort2.append(pairs2[p])
    threshold1 = numpy.mean(to_sort1) + numpy.std(to_sort1)*2
    threshold2 = numpy.mean(to_sort2) + numpy.std(to_sort2)*2

    designed = []
    natural = []
    overlapping = []
    for p in pairs1:
        if p in pairs2:
            i = p.split()[0]
            j = p.split()[1]
            if pairs1[p] > threshold1:
                a_count += 1
                natural.append(p)
                positions.add(int(i))
                positions.add(int(j))
            if pairs2[p] > threshold2:
                b_count += 1
                designed.append(p)
            if pairs1[p] > threshold1 and pairs2[p] > threshold2:
                overlap += 1
                overlapping.append(p)
            total += 1
    return natural, positions


background = {
    'A' : 0.0853130414059,
    'C' : 0.0145091808885,
    'E' : 0.0697042211031,
    'D' : 0.0576610517405,
    'G' : 0.0677683836625,
    'F' : 0.0368651011894,
    'I' : 0.0658157819481,
    'H' : 0.0211289643495,
    'K' : 0.0581917850968,
    'M' : 0.0190262038642,
    'L' : 0.0958638899794,
    'N' : 0.0369395202374,
    'Q' : 0.036293485414,
    'P' : 0.0391082335344,
    'S' : 0.0594265039867,
    'R' : 0.0562652852139,
    'T' : 0.0541996845528,
    'W' : 0.0108604669712,
    'V' : 0.0866667775459,
    'Y' : 0.0283924373158
}


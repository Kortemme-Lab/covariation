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
Example usage:
python covariation_similarity.py ../../output/fixed_backbone/ ../../output/backrub/ > covariation_similarity.txt
"""


import os
import sys
import math
import operator
import numpy
from scipy.stats import hypergeom

def compute_overlap(nature_mi, design_mi, indicefile):
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
        if threshold1 > z:
            threshold1 = z
        #if dist <= 6:
        dists[i+" "+j] = dist
        pairs1[i+" "+j] = z
        natural_scores.append(z)
    threshold1 = -1 * threshold1

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
        if threshold2 > z:
            threshold2 = z
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
            positions.add(i)
            positions.add(j)
            to_sort1.append(pairs1[p])
            to_sort2.append(pairs2[p])

    threshold1 = numpy.mean(to_sort1) + numpy.std(to_sort1)*2
    threshold2 = numpy.mean(to_sort2) + numpy.std(to_sort2)*2

    designed = []
    natural = []
    overlapping = []
    indices = {}
    f = open(indicefile, 'r')
    for line in f:
        indices[line.split()[0]] = line.split()[3]
    f.close()

    for p in pairs1:
        if p in pairs2:
            i = p.split()[0]
            j = p.split()[1]
            if pairs1[p] > threshold1:
                a_count += 1
                natural.append(p)
            if pairs2[p] > threshold2:
                b_count += 1
                designed.append(p)
            if pairs1[p] > threshold1 and pairs2[p] > threshold2:
                overlap += 1
                overlapping.append(p)
            total += 1

    return float(overlap)*2 / float(a_count+b_count)


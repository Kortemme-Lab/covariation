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

import sys
import os
import math
import operator
from fsio import read_file

def get_MI(length, entropies, joint_entropies):

    MI = {}
    MI_sums = {}
    MI_means = {}
    MI_sums_correct = {}
    for i in range(0, length):
        MI_sums[i] = 0
        MI_sums_correct[i] = 0
    sum_MI = 0

    for i in range(0, length):
        for j in range(i+1, length):
            ij = str(i)+" "+str(j)
            MI[ij] = (entropies[i] + entropies[j] - joint_entropies[ij])
            MI_sums[i] += MI[ij]
            MI_sums[j] += MI[ij]
            sum_MI += MI[ij]

    mean_MI = sum_MI / float(len(MI))
    #print mean_MI

    std_sum = 0
    for ij in MI:
        std_sum += (MI[ij] - mean_MI)**2
    std_MI = math.sqrt(std_sum / float(len(MI)))
    Z = {}
    for ij in MI:
        Z[ij] = (MI[ij] - mean_MI) / std_MI

    for i in range(0, length):
        MI_means[i] = MI_sums[i] / float(length - 1)

    sum_corrected_MI = 0
    APC = {}
    corrected_MI = {}

    for ij in MI:
        i = int(ij.split()[0])
        j = int(ij.split()[1])
        APC[ij] = (MI_means[i] * MI_means[j]) / mean_MI
        corrected_MI[ij] = MI[ij] - (MI_means[i] * MI_means[j]) / mean_MI
        sum_corrected_MI += corrected_MI[ij]
    mean_corrected_MI = sum_corrected_MI / float(len(MI))
    sum_corrected_MI = 0

    for ij in corrected_MI:
        sum_corrected_MI += (corrected_MI[ij] - mean_corrected_MI)**2

    stddev_corrected_MI = math.sqrt(sum_corrected_MI / float(len(MI)))

    Zp = {}
    for ij in corrected_MI:
        Zp[ij] = (corrected_MI[ij] - mean_corrected_MI) / stddev_corrected_MI

    mean_corrected_MI_i = {}
    for i in range(0, length):
        sum = float(0.0)
        count = float(0.0)
        for j in range(0, length):
            if i != j:
                ij = str(i)+" "+str(j)
                ji = str(j)+" "+str(i)
                if ij in corrected_MI:
                    count += 1.0
                    sum += corrected_MI[ij]
                else:
                    count += 1.0
                    sum += corrected_MI[ji]
        mean_corrected_MI_i[i] = sum / count

    stddev_corrected_MI_i = {}
    for i in range(0, length):
        sum = float(0.0)
        count = float(0.0)
        for j in range(0, length):
            if i != j:
                ij = str(i)+" "+str(j)
                ji = str(j)+" "+str(i)
                if ij in corrected_MI:
                    count += 1.0
                    sum += (mean_corrected_MI_i[i] - corrected_MI[ij])**2
                else:
                    count += 1.0
                    sum += (mean_corrected_MI_i[i] - corrected_MI[ji])**2

        stddev_corrected_MI_i[i] = math.sqrt(sum / count)

    Zpx = {}
    for ij in corrected_MI:
        i = int(ij.split()[0])
        j = int(ij.split()[1])
        if stddev_corrected_MI_i[i] > 0 and stddev_corrected_MI_i[j] > 0:
            Zpx[ij] = ((corrected_MI[ij] - mean_corrected_MI_i[i]) / stddev_corrected_MI_i[i]) * ((corrected_MI[ij] - mean_corrected_MI_i[j]) / stddev_corrected_MI_i[j])
        else:
            Zpx[ij] = 0

    return MI, Z, Zp, Zpx


def get_domain_sequences(domain_sequences_file_content):
    structs = {}
    for line in domain_sequences_file_content.split('\n'):
        line = line.strip()
        if line:
            structs[line.split()[0]] = line.split()[1]
    return structs


def create_mi_file(domain, fasta_file_contents, domain_sequences, indices_directory, expectn = None):

    structure = domain_sequences[domain]
    domain_indices = read_file(os.path.join(indices_directory, '{0}.indices'.format(structure)))
    pfam_indices = {}
    for line in domain_indices.split('\n'):
        if line.strip():
            pfam_indices[int(line.split()[2])-1] = line.split()[0]
    return compute_mi(pfam_indices, fasta_file_contents, expectn = expectn)


def compute_mi(pfam_indices, fasta_file_contents, expectn = None):

    gapped_sequences = set()
    gapped_length = 0
    length = 0
    aa = 'ACDEFGHIKLMNPQRSTVWY'
    name = ''
    seq = ''
    test = []
    for line in fasta_file_contents.split('\n'):
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
    sequences = set()
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
    for seq in gapped_sequences:
        ungapped_seq = ''
        count = 0
        for i in ungapped:
            if i in pfam_indices:
                ungapped_seq += seq[i]
                ungapped_indices[count] = pfam_indices[i]
                count += 1
        sequences.add(ungapped_seq)
        length = len(ungapped_seq)

    single = {}
    for a in aa:
        single[a] = 0

    pair = {}
    for a in aa:
        for b in aa:
            pair[a+b] = 0

    num_sequences = float(len(sequences))
    if expectn and len(sequences) != expectn:
        raise Exception('Expected {0} records in but read {2}.'.format(expectn, len(sequences)))

    entropies = {}
    counts = {}
    frequencies = {}

    for i in range(0, length):
        counts[i] = single.copy()
        for seq in sequences:
            counts[i][seq[i]] += 1
        frequencies[i] = single.copy()
        sum = 0
        for char in frequencies[i]:
            freq = float(counts[i][char]) / num_sequences
            frequencies[i][char] = freq
            if freq > 0:
                sum += -1 * freq * math.log(freq,20)
        entropies[i] = sum

    joint_entropies = {}

    for i in range(0, length):
        for j in range(i+1, length):
            joint_counts = pair.copy()
            ab = str(i)+" "+str(j)
            for seq in sequences:
                joint_counts[seq[i]+seq[j]] += 1
            sum_joint_entropy = 0
            for char in joint_counts:
                freq = float(joint_counts[char]) / num_sequences
                if freq > 0:
                    sum_joint_entropy += -1 * freq * math.log(freq,20)
            joint_entropies[ab] = sum_joint_entropy

    MI, Z, Zp, Zpx = get_MI(length, entropies, joint_entropies)

    s = []
    for a,b in sorted(Z.items(), key = operator.itemgetter(1), reverse = False):
        i = int(a.split()[0])
        j = int(a.split()[1])
        if Zp[a] < 0:
            Zpx[a] = -1 * Zpx[a]
        s.append(ungapped_indices[i]+"\t"+ungapped_indices[j]+"\t"+str(Z[a])+"\t"+str(Zp[a])+"\t"+str(Zpx[a])+"\t"+str(joint_entropies[a])+"\t"+str(entropies[i])+"\t"+str(entropies[j])+"\n")
    return ''.join(s)


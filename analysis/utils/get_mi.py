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
import numpy
from fsio import read_file

canonical_aas = 'ACDEFGHIKLMNPQRSTVWY'

class SequenceMatrix(object):
    '''A class to contain a list of sequences and present them as a numpy matrix for processing.'''

    def __init__(self, sequences, expectn = None, remove_gap_indices = False):
        '''Checks to make sure that we can compute values over the sequences.'''

        # Create a indexing from amino acid type to columns
        aa_matrix_column_indices = {}
        aas_ = sorted(list(canonical_aas))
        for x in range(len(aas_)):
            aa_matrix_column_indices[aas_[x]] = x
        self.aa_matrix_column_indices = aa_matrix_column_indices
        self.remove_gap_indices = remove_gap_indices
        self.expectn = expectn
        self.sequences = sequences
        self.matrix = None
        self.gap_indices = []
        self.update()


    def update(self):
        '''Converts the sequences into a numpy matrix.'''

        # Make sure that all sequences have the same length
        sequences = self.sequences
        length = set([len(s) for s in sequences])
        assert(len(length) == 1)
        length = length.pop()
        positions = range(length)

        # If there are gaps, make sure that the gaps always occur in the same positions
        gap_frequencies = dict.fromkeys(positions, 0)
        for s in sequences:
            for i in positions:
                if s[i] == '-':
                    gap_frequencies[i] += 1
                    gaps_exist = True
        gap_indices = [i for i, v in gap_frequencies.iteritems() if v]
        for i in positions:
            assert(gap_frequencies[i] == 0 or gap_frequencies[i] == len(sequences))
        if gap_indices and self.remove_gap_indices:
            new_sequences = [s.replace('-', '') for s in sequences]
            sequences = new_sequences
            length = set([len(s) for s in sequences])
            assert(len(length) == 1)
            length = length.pop()
            gap_indices = [] # we remove the indexed positions from the sequences

        # Make sure that we have the correct number of sequences
        if self.expectn and len(sequences) != self.expectn:
            raise Exception('Expected {0} records in but read {2}.'.format(self.expectn, len(sequences)))

        self.gap_indices = gap_indices
        self.sequences = sequences
        self.matrix = numpy.array([list(seq) for seq in sequences])


    def get_sequence_entropy(self):
        '''Calculates the sequence entropy for the list of sequences. If gaps exists then they must exist at the same positions
           in all sequences. If remove_gap_indices is True then gap positions will be ignored and the indices adjusted accordingly.
           If remove_gap_indices is False then the sequence entropy for the gapped positions will be included and set to None.'''

        sequence_matrix, gap_indices, aa_matrix_column_indices = self.matrix, self.gap_indices, self.aa_matrix_column_indices
        num_sequences, sequence_length = numpy.shape(sequence_matrix)
        positions = range(sequence_length)

        # Calculate sequence entropy
        with numpy.errstate(divide='ignore', invalid='ignore'):
            count_matrix = numpy.zeros((sequence_length, 20))
            for i in positions:
                if i not in gap_indices:
                    position_aas = sequence_matrix[:,i]
                    unique, counts = numpy.unique(position_aas, return_counts=True) # this returns the SORTED unique elements of the array so we can use counts directly as
                    for cx in range(len(unique)):
                         # todo: It may be quicker to add the row in one operation. Benchmark this.
                         count_matrix[i][aa_matrix_column_indices[unique[cx]]] = counts[cx]
            count_matrix = count_matrix / num_sequences
            count_matrix = -1.0 * count_matrix * (numpy.log(count_matrix)/numpy.log(20))
            entropies = numpy.nansum(count_matrix, axis = 1) # sum up the rows (sequence positions), treating NaN as zero
            for i in positions:
                if i in gap_indices:
                    entropies[i] = None
            entropies = dict(zip(positions, entropies.tolist())) # convert the vector to a dict indexed by position. Note: we use numpy.nan instead of None as before so we may hit some errors
        return entropies


def get_MI(length, entropies, joint_entropies):
    '''length is the number of positions in the sequences.
       entropies is mapping from positions i to the Shannon entropy of the positions (H_i in the paper)
       joint_entropies is a mapping from pairs of positions i, j where i < j to their joint_entropy.  (H_i,j in the paper)
    '''

    assert(len(entropies) == length)

    MI = {}
    MI_sums = {}
    MI_means = {}
    for i in range(0, length):
        MI_sums[i] = 0
    sum_MI = 0

    for i in range(0, length):
        for j in range(i+1, length):
            ij = str(i)+" "+str(j)
            MI[ij] = (entropies[i] + entropies[j] - joint_entropies[ij]) # This computes MI_i,j as described in the paper
            MI_sums[i] += MI[ij]
            MI_sums[j] += MI[ij]
            sum_MI += MI[ij]

    # Compute \bar{MI}, the mean MI over all pairs of positions
    mean_MI = sum_MI / float(len(MI))

    # Compute \bar{MI_i}, the mean MI of position i with all other positions
    for i in range(0, length):
        MI_means[i] = MI_sums[i] / float(length - 1)

    # Compute corrected_MI = MIp_i,j for all i and j and \bar{MIp}, the average MIp over all pairs i and j
    sum_corrected_MI = 0
    APC = {}
    corrected_MI = {}
    for ij in MI:
        i = int(ij.split()[0])
        j = int(ij.split()[1])
        APC[ij] = (MI_means[i] * MI_means[j]) / mean_MI
        corrected_MI[ij] = MI[ij] - (MI_means[i] * MI_means[j]) / mean_MI # this computes MIp_i,j, the product corrected mutual information
        sum_corrected_MI += corrected_MI[ij]
    mean_corrected_MI = sum_corrected_MI / float(len(MI)) # this computes the average MIp over all pairs i and j, \bar{MIp}

    # Compute the standard deviation of MIp
    sum_corrected_MI = 0
    for ij in corrected_MI:
        sum_corrected_MI += (corrected_MI[ij] - mean_corrected_MI)**2
    stddev_corrected_MI = math.sqrt(sum_corrected_MI / float(len(MI)))

    # todo: resume commenting from here

    # Compute Z_i,j
    std_sum = 0
    for ij in MI:
        std_sum += (MI[ij] - mean_MI)**2
    std_MI = math.sqrt(std_sum / float(len(MI)))
    Z = {}
    for ij in MI:
        Z[ij] = (MI[ij] - mean_MI) / std_MI


    Zp = {}
    for ij in corrected_MI:
        Zp[ij] = (corrected_MI[ij] - mean_corrected_MI) / stddev_corrected_MI

    # Compute mean_corrected_MI_i = \bar{MIp_i}, the mean of MIp_i,j where j ranges over all positions where i != j
    mean_corrected_MI_i = {}
    for i in range(0, length):
        xsum = float(0.0)
        count = float(0.0)
        for j in range(0, length):
            if i != j:
                ij = str(i)+" "+str(j)
                ji = str(j)+" "+str(i)
                if ij in corrected_MI:
                    count += 1.0
                    xsum += corrected_MI[ij]
                else:
                    count += 1.0
                    xsum += corrected_MI[ji]
        mean_corrected_MI_i[i] = xsum / count

    stddev_corrected_MI_i = {}
    for i in range(0, length):
        xsum = float(0.0)
        count = float(0.0)
        for j in range(0, length):
            if i != j:
                ij = str(i)+" "+str(j)
                ji = str(j)+" "+str(i)
                if ij in corrected_MI:
                    count += 1.0
                    xsum += (mean_corrected_MI_i[i] - corrected_MI[ij])**2
                else:
                    count += 1.0
                    xsum += (mean_corrected_MI_i[i] - corrected_MI[ji])**2

        stddev_corrected_MI_i[i] = math.sqrt(xsum / count)

    Zpx = {}
    for ij in corrected_MI:
        i = int(ij.split()[0])
        j = int(ij.split()[1])
        if stddev_corrected_MI_i[i] > 0 and stddev_corrected_MI_i[j] > 0:
            # corrected_MI[ij] is MIp_i,j
            Zpx[ij] = ((corrected_MI[ij] - mean_corrected_MI_i[i]) / stddev_corrected_MI_i[i]) * ((corrected_MI[ij] - mean_corrected_MI_i[j]) / stddev_corrected_MI_i[j])
        else:
            Zpx[ij] = 0

    return MI, Z, Zp, Zpx


def create_mi_file(domain, fasta_file_contents, domain_sequences, indices_directory, expectn = None):

    structure = domain_sequences[domain]
    domain_indices = read_file(os.path.join(indices_directory, '{0}.indices'.format(structure)))
    pfam_indices = {}
    for line in domain_indices.split('\n'):
        if line.strip():
            # Column 3 in the indices file is the index into the PDB sequence but this is 1-indexed so we convert it to 0-indexed here
            pfam_indices[int(line.split()[2])-1] = line.split()[0]
    return compute_mi(pfam_indices, fasta_file_contents, expectn = expectn, domain = domain)


def compute_mi(pfam_indices, fasta_file_contents, expectn = None, domain = None):

    # pfam_indices is a mapping from PDB sequence indices (zero-indexed) to natural sequence indices (zero-indexed)

    gapped_sequences = set()
    gapped_length = 0
    length = 0
    aa = 'ACDEFGHIKLMNPQRSTVWY'
    name = ''
    seq = ''
    test = []

    # Read in the fasta file contents. gapped_sequences should contain expectn sequences of the same length (typically
    # the number of residues in the original PDB file unless some residues were dropped
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

    # Makes sure that all sequences are the same length and remove any that are not
    # Note: This is somewhat arbitrary since the length of the last sequence is taken. It may be better to raise an
    #       exception in this case
    gapped_length = len(seq) # gapped_length is the length of the Rosetta sequences (usually the same as the length of the input PDB sequence)
    ungapped = []
    to_remove = set()
    for seq in gapped_sequences:
        if len(seq) != gapped_length:
            to_remove.add(seq)
    for seq in to_remove:
        gapped_sequences.remove(seq)

    # Compile a zero-indexed list of positions ("ungapped") which only contain canonical residues
    for i in range(0, gapped_length):
        counter = 0
        for seq in gapped_sequences:
            char = seq[i]
            if char not in aa:
                counter += 1
        if counter < 1:
            ungapped.append(i)
    length = len(ungapped)

    # Iterate through the list of sequences and, for all positions containing only canonical residues, extract the sequence if there is a mapping to the natural sequences.
    # This creates a set of truncated sequences (the sequences variable) from the designs with corresponding resiudes in the natural sequences.
    # Note: this seems dangerous if there are any positions in the FASTA file which contain canonical residues which were removed above
    # ungapped_seq will contain sequences with only canonical residues which can be mapped to the native sequences
    # Assuming all sequences are valid at the same positions, ungapped_indices will be a mapping from indices in the (truncated) sequences set to indices in the natural sequence MSA
    sequences = set()
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

    # A counter for occurrences of an canonical amino acid
    single = {}
    for a in aa:
        single[a] = 0

    # A counter for paired-occurrences of canonical amino acids
    pair = {}
    for a in aa:
        for b in aa:
            pair[a+b] = 0

    # Make sure that we have the expected number of sequences
    num_sequences = float(len(sequences))
    if expectn and len(sequences) != expectn:
        raise Exception('Expected {0} records in but read {2}.'.format(expectn, len(sequences)))

    # This is the calculation of mutual information

    # Compute the Shannon entropy at each position
    entropies = {}
    counts = {}
    frequencies = {}
    for i in range(0, length):
        counts[i] = single.copy()
        for seq in sequences:
            counts[i][seq[i]] += 1
        frequencies[i] = single.copy()
        xsum = 0
        for char in frequencies[i]:
            freq = float(counts[i][char]) / num_sequences
            frequencies[i][char] = freq
            if freq > 0:
                xsum += -1 * freq * math.log(freq,20)
        entropies[i] = xsum

    # The joint entropy calculated at all unique pairs of distinct positions
    # If there are n positions then there will be ((n - 1) + (n - 2) + ... + (n - n)) pairs of positions
    # The keys of joint_entropies are internally sorted i.e. the first word (sequence index) in the key will always be lower than the second word
    joint_entropies = {}
    for i in range(0, length):
        for j in range(i+1, length):
            joint_counts = pair.copy()
            ab = str(i)+" "+str(j)
            for seq in sequences:
                joint_counts[seq[i]+seq[j]] += 1
            sum_joint_entropy = 0
            for char_pair in joint_counts:
                freq = float(joint_counts[char_pair]) / num_sequences
                if freq > 0:
                    sum_joint_entropy += -1 * freq * math.log(freq,20)
            joint_entropies[ab] = sum_joint_entropy

    # Z has an entry for each pair in joint_entropies (each unique pair of distinct positions) and is a mapping
    #    'i j' -> float where i and j are positions in the sequence and i < j
    MI, Z, Zp, Zpx = get_MI(length, entropies, joint_entropies)
    assert(sorted(joint_entropies) == sorted(Z.keys()) == sorted(Zp.keys()) == sorted(Zpx.keys()))

    s = []
    for a,b in sorted(Z.items(), key = operator.itemgetter(1), reverse = False):
        i = int(a.split()[0])
        j = int(a.split()[1])
        if Zp[a] < 0:
            Zpx[a] = -1 * Zpx[a]
        s.append(ungapped_indices[i]+"\t"+ungapped_indices[j]+"\t"+str(Z[a])+"\t"+str(Zp[a])+"\t"+str(Zpx[a])+"\t"+str(joint_entropies[a])+"\t"+str(entropies[i])+"\t"+str(entropies[j])+"\n")

    natural_indexed_residue_entropies = {}
    designed_sequence_matrix = SequenceMatrix(list(sequences))
    rosetta_indexed_residue_entropies = designed_sequence_matrix.get_sequence_entropy()
    for rosetta_res_id, entropy in rosetta_indexed_residue_entropies.iteritems():
        assert(abs(entropy - entropies[rosetta_res_id]) < 0.0001) # todo: sanity check - we can remove this later once we remove the double calculation of entropy
        natural_indexed_residue_entropies[int(ungapped_indices[rosetta_res_id])] = entropy

    return ''.join(s), entropies, natural_indexed_residue_entropies


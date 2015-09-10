====================================
Analysis scripts
====================================

This benchmark includes Python scripts to calculate the following metrics:

- Covariation similarity (analysis/covariation_similarity/covariation_similarity.py)
- Profile similarity (analysis/profile_similarity/profile_similarity.py)
- Native sequence recovery (analysis/sequence_recovery/sequence_recovery.py)

The usage of these scripts is as follows:

::

  python covariation_similarity.py ../../output/fixed_backbone/ ../../output/backrub/ > covariation_similarity.txt
  python profile_similarity.py ../../output/fixed_backbone/ ../../output/backrub/ > profile_similarity.txt
  python sequence_recovery.py ../../output/fixed_backbone/ ../../output/backrub/ > sequence_recovery.txt

R scripts are provided to visualize the output. The R scripts can be run as follows:

::

  R CMD BATCH covariation_similarity/covariation_similarity.r
  R CMD BATCH profile_similarity/profile_similarity.r
  R CMD BATCH sequence_recovery/sequence_recovery.r


====================================
Mutual information files (.mi)
====================================

The analysis scripts output .fasta files containing the predicted sequences and .mi files which contain the mutual
information data. These values are calculated as in the publication. The .mi files are created based on .indices files
contained in the indices directory. There is one .indices file per domain which is used to map the predicted sequence to
the native sequence. Warning: this mapping is based on the residues created by Rosetta. Other methods may require different
.indices files depending on which residues are retained and discarded.

There are two inputs used to create a mutual information file (utils/get_mi.py::compute_mi):
 - the set of predicted sequences (the .fasta file);
 - a mapping from PDB sequence indices (zero-indexed) to natural sequence indices (zero-indexed) which is extracted from the .indices file.

The creation of a mutual information file is as follows:
 - a new set of sequences - of equal length - is created which only contain those positions which correspond to positions in the .indices file (ungapped positions with canonical amino acids);
 - a mapping from indices in the (truncated) predicted sequences to indices in the natural sequence MSA is created (ungapped_indices);
 - the Shannon entropy is calculated at each position, considering all truncated sequences;
 - the joint entropy calculated at all unique pairs of distinct positions - if there are n positions then there will be ((n - 1) + (n - 2) + ... + (n - n)) pairs of positions;
 - the Zpx value described in the paper and other calculations used to determine Zpx.

The mutual information (.mi) files are structured as follows. Each line describes values for a pair of positions in the sequence. There are eight columns:
 - column 1: the left position, i;
 - column 2: the right position, j (i < j);
 - column 3: Z_i,j
 - column 4: Zp_i,j
 - column 5: Zpx_i,j
 - column 6: the joint entropy for i, j
 - column 7: the entropy for position i
 - column 8: the entropy for position j

e.g.

15      55      -1.67634845166  -0.00661831781105       -0.0434951573329        0.602062381127  0.0     0.602062381127
2       15      -1.67634845166  -0.00661831781105       -0.0697620569844        0.727410761725  0.727410761725  0.0
15      54      -1.67634845166  -0.00661831781103       -0.0229734438946        0.68562437401   0.0     0.68562437401


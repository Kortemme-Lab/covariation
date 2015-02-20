====================================
Analysis scripts
====================================

This benchmark includes Python scripts to calculate the following metrics:

- Covariation similarity (analysis/covariation_similarity/covariation_similarity.py)
- Profile similarity (analysis/profile_similarity/profile_similarity.py)
- Native sequence recovery (analysis/sequence_recovery/sequence_recovery.py)

The usage of these scripts is as follows:

python covariation_similarity.py ../../output/fixed_backbone/ ../../output/backrub/ > covariation_similarity.txt

python profile_similarity.py ../../output/fixed_backbone/ ../../output/backrub/ > profile_similarity.txt

python sequence_recovery.py ../../output/fixbb_sequences/ ../../output/backrub_sequences/ > sequence_recovery.txt

R scripts are provided to visualize the output. The R scripts can be run as follows:

R CMD BATCH covariation_similarity/covariation_similarity.r
R CMD BATCH profile_similarity/profile_similarity.r
R CMD BATCH sequence_recovery/sequence_recovery.r

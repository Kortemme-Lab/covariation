#!/bin/bash
#$ -cwd
#$ -r y
#$ -j y
#$ -l h_rt=24:00:00
#$ -t 1-500
#$ -l mem_free=2G
#$ -l arch=linux-x64

date
hostname

rosetta/backrub.linuxgccrelease -database rosetta/database -s $1.pdb -resfile NATAA.res -ex1 -ex2 -extrachi_cutoff 0 -backrub:mc_kt $2 -backrub:ntrials 10000 -nstruct 1 -out::suffix _$2_$SGE_TASK_ID -backrub:initial_pack

rosetta/fixbb.linuxgccrelease -database rosetta/database -s $1_$2_${SGE_TASK_ID}_0001_last.pdb -resfile ALLAA.res -ex1 -ex2 -extrachi_cutoff 0 -nstruct 1 -overwrite -linmem_ig 10 -no_his_his_pairE -minimize_sidechains

date

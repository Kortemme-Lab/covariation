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

rosetta/fixbb.linuxgccrelease -database rosetta/database -s $1 –resfile ALLAA.res –ex1 –ex2 –extrachi_cutoff 0 –nstruct 500 –linmem_ig 10 –no_his_his_pairE –minimize_sidechains

date

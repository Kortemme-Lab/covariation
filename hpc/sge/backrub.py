#!/usr/bin/env python2

#$ -S /usr/bin/python
#$ -cwd
#$ -r y
#$ -j y
#$ -l h_rt=24:00:00
#$ -t 1-500
#$ -l mem_free=2G
#$ -l arch=linux-x64

import sys
import os
import socket
import datetime
import shlex
import subprocess

ntrials = 10000

print "Python version:", sys.version
print "Hostname:", socket.gethostname()
print "Time:", datetime.datetime.now()
sge_task_id = long(os.environ["SGE_TASK_ID"])
print "Task:", sge_task_id

rosetta_env = os.environ.copy()
mysql_lib = '/netapp/home/kbarlow/lib/mysql-connector-c-6.1.2-linux-glibc2.5-x86_64/lib:'
try:
    rosetta_env['LD_LIBRARY_PATH'] = mysql_lib + ':' + rosetta_env['LD_LIBRARY_PATH']
except KeyError:
    rosetta_env['LD_LIBRARY_PATH'] = mysql_lib

args = shlex.split('/netapp/home/shaneoconner/compilation/r57934/main/source/bin/backrub.linuxgccrelease -database /netapp/home/shaneoconner/compilation/r57934/main/database -s %(input_pdb)s.pdb -resfile /netapp/home/shaneoconner/t14benchmarking/covariation/input/NATAA.res -ex1 -ex2 -extrachi_cutoff 0 -backrub:mc_kt %(temperature)s -backrub:ntrials ' + str(ntrials) + ' -nstruct 1 -out::suffix _%(temperature)s_' + str(sge_task_id) + ' -backrub:initial_pack')

rosetta_process = subprocess.Popen(args, stdout=subprocess.PIPE, cwd="%(job_directory)s", env=rosetta_env)
out, err = rosetta_process.communicate()
sys.stdout.write(out or '')
sys.stdout.flush()
if err:
    sys.stderr.write(err)
return_code = rosetta_process.returncode
print "Backrub return code", return_code
print ""

args = shlex.split('/netapp/home/shaneoconner/compilation/r57934/main/source/bin/fixbb.linuxgccrelease -database /netapp/home/shaneoconner/compilation/r57934/main/database -s %(input_pdb)s_%(temperature)s_' + str(sge_task_id) + '_0001_last.pdb -resfile /netapp/home/shaneoconner/t14benchmarking/covariation/input/ALLAA.res -ex1 -ex2 -extrachi_cutoff 0 -nstruct 1 -overwrite -linmem_ig 10 -no_his_his_pairE -minimize_sidechains')

rosetta_process = subprocess.Popen(args, stdout=subprocess.PIPE, cwd="%(job_directory)s", env=rosetta_env)
out, err = rosetta_process.communicate()
sys.stdout.write(out or '')
sys.stdout.flush()
if err:
    sys.stderr.write(err)
return_code = rosetta_process.returncode
print "Fixbb return code", return_code
print ""

print "Time:", datetime.datetime.now()


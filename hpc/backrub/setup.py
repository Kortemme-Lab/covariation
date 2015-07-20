#!/usr/bin/python

import sys
import os
import glob
import shutil
import subprocess
import shlex

def write_file(filepath, contents, ftype = 'w'):
    output_handle = open(filepath, ftype)
    output_handle.write(contents)
    output_handle.close()


input_pdbs = []
input_directory = os.path.join('..', '..', 'input', 'domains')
output_directory = os.path.join('job_output')

job_template = open('backrub.py').read()
for domain in os.listdir(input_directory):
    for pdb in glob.glob(os.path.join(input_directory, domain, '*.pdb')):
        pdb_name = os.path.splitext(os.path.split(pdb)[1])[0]
        job_name = '%s_%s' % (domain, pdb_name)
        input_pdbs.append(os.path.abspath(pdb))
        if not os.path.exists(job_name):
            os.mkdir(job_name)
        job_details = dict(
            temperature = '0.9',
            input_pdb = pdb_name,
            job_directory = (os.path.abspath(job_name))
        )
        shutil.copyfile(pdb, os.path.join(job_name, pdb_name + '.pdb'))
        job_script_path = os.path.join(job_name, 'backrub.py')
        write_file(job_script_path, job_template % job_details)
        subprocess.call(shlex.split('qsub ' + job_script_path))

if len(input_pdbs) == 0:
    sys.exit('No input PDB files were found.')

print('\nSetting up jobs for {0} domains.\n'.format(len(input_pdbs)))

# | grep pdb


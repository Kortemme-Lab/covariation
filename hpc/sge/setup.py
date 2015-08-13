#!/usr/bin/env python2

# The MIT License (MIT)
#
# Copyright (c) 2015 Shane O'Connor
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

methods = dict(
    backrub = dict(
        method_type = 'Flexible backbone',
        args = dict(
            backrub_step = '"%(backrub_executable)s -database %(rosetta_database_path)s -s %(input_pdb)s.pdb -resfile %(NATAA_resfile)s -ex1 -ex2 -extrachi_cutoff 0 -backrub:mc_kt %(temperature)s -backrub:ntrials %(ntrials)d -nstruct 1 -backrub:initial_pack -out::suffix _%(temperature)s_" + str(sge_task_id) + " %(extra_flags)s"',
            fixbb_step   = '"%(fixbb_executable)s -database %(rosetta_database_path)s -s %(input_pdb)s_%(temperature)s_" + str(sge_task_id) + "_0001_last.pdb -resfile %(ALLAA_resfile)s -ex1 -ex2 -extrachi_cutoff 0 -nstruct 1 -overwrite -linmem_ig 10 -no_his_his_pairE -minimize_sidechains %(extra_flags)s"'
        ),
        dependent_binaries = ['backrub', 'fixbb'],
    ),
    fixbb = dict(
        method_type = 'Fixed backbone',
        args = dict(
            fixbb_step   = '"%(fixbb_executable)s -database %(rosetta_database_path)s -s %(input_pdb)s.pdb -resfile %(ALLAA_resfile)s -ex1 -ex2 -extrachi_cutoff 0 -nstruct 1 -overwrite -linmem_ig 10 -no_his_his_pairE -minimize_sidechains -out::suffix _" + str(sge_task_id) + " %(extra_flags)s"',
        ),
        dependent_binaries = ['fixbb'],
    ),
)


__doc__ = """\
This script can be used to setup, and optionally run, the covariation benchmark. The option below
allow you to specify different design methods, the amount of sampling, temperatures (for flexible
backbone sampling), and additional Rosetta options.
Example command lines:
    setup.py backrub
    setup.py backrub --temperature 1.2 --output_directory testrun --domains PF00013 --domains PF00018


Usage:
    setup.py <method> [--domains=DOMAIN1 ...] [options]

Arguments:

    <method>
        A Rosetta method to execute. The methods currently available are:
          """ + ', '.join(methods.keys()) + """
        This file can be edited to add more methods.

    <domains>
        A list of domains to use for the benchmark. By default, all directories in input/domains
        will be searched for files with the extension '*.pdb' and these files will be used in the
        benchmark run. This option allows the user to select a subset of these domains.

Options:

    -o --output_directory OUTPUT_DIR
        The path where output data will be created. Output will be created inside a time-stamped subfolder of this
        directory. [default: ./job_output]

    --test
        When this option is set, a shorter version of the benchmark will run with fewer input domains and fewer
        generated structures. This should be used to test the scripts rather than for benchmarking.

    --talaris2014
        When this option is set, the talaris2014 score function will be used rather than the default score function.
        Warning: This option may break when talaris2014 becomes the default Rosetta score function.

    --flags OPT
        Specify a rosetta flag file containing extra options for this run.

    --nstruct NUM -n NUM
        Specify how many simulations to do for each structure in the benchmark.
        The default value is 500. [default: 500]

    --ntrials NUM
        Specify how many simulations to do for each structure in the benchmark.
        The default value is 10000. [default: 10000]

    --temperature TEMP
        Specify the value of kT for flexible backbone modeling. [default: 0.9]

    --mysql_lib PATH_TO_MYSQL_LIBRARY
        If Rosetta is built by linking against MySQL libraries, these may need to be added to the LD library path. This
        option allows you to specify the location of the libraries.

    --setup_only
        If this option is selected, the jobs are not submitted to the cluster. The files and execution script are created
        so that the job may be run later.

    --settings SETTINGS_FILE
        By default, protocols/rosetta/settings.json is used the load the benchmark settings (paths to Rosetta etc.). This
        location can be overridden with this option.
"""

import sys
import os
import glob
import shutil
import subprocess
import shlex
import copy
import json
import pprint
from libraries.utilities import *
from libraries import docopt


script_preamble = '''#!/usr/bin/env python2

#$ -S /usr/bin/python
#$ -cwd
#$ -r y
#$ -j y
#$ -l h_rt=24:00:00
#$ -t 1-%(nstruct)s
#$ -l mem_free=2G
#$ -l arch=linux-x64

import sys
import os
import socket
import datetime
import shlex
import subprocess

print "Python version:", sys.version
print "Hostname:", socket.gethostname()
print "Time:", datetime.datetime.now()
sge_task_id = long(os.environ["SGE_TASK_ID"])
print "Task:", sge_task_id

rosetta_env = os.environ.copy()
%(ld_path_extension)s

# Benchmarking the %(method_type)s method
'''

mysql_path_extension = '''
# Add the MySQL library path to the LD path
mysql_lib = '{0}'
try:
    rosetta_env['LD_LIBRARY_PATH'] = mysql_lib + ':' + rosetta_env['LD_LIBRARY_PATH']
except KeyError:
    rosetta_env['LD_LIBRARY_PATH'] = mysql_lib'''

simple_step_commands = '''
# %(step)s step
args = shlex.split(%(command_line)s)
rosetta_process = subprocess.Popen(args, stdout=subprocess.PIPE, cwd="%(job_directory)s", env=rosetta_env)
out, err = rosetta_process.communicate()
sys.stdout.write(out or '')
sys.stdout.flush()
if err:
    sys.stderr.write(err)
return_code = rosetta_process.returncode
print "%(method)s return code", return_code
print ""
'''

script_epilog = '''
print "Time:", datetime.datetime.now()
'''


def create_simple_step_command_line(method, step, command_line, job_directory):
    return simple_step_commands % locals()


def create_script(job_script_path, job_parameters):

    job_parameters['method_type'] = methods[job_parameters['method']]['method_type']
    method_details = methods[job_parameters['method']]
    script = script_preamble % job_parameters
    for step in method_details['dependent_binaries']:
        command_line = method_details['args']['{0}_step'.format(step)] % job_parameters
        script += create_simple_step_command_line(job_parameters['method'], step, command_line, job_parameters['job_directory'])
    script += script_epilog
    write_file(job_script_path, script)


if __name__ == '__main__':
    print('')
    benchmark_root = os.path.abspath(os.path.join('..', '..'))
    arguments = docopt.docopt(__doc__)

    method = arguments['<method>']
    restricted_to_domains = arguments['--domains']

    if not method in methods:
        die('Invalid method. Possible choices are: {0}.'.format(', '.join(sorted(methods.keys()))))
    method_details = methods[method]

    # Read the configuration file
    default_settings_filepath = os.path.join(benchmark_root, 'protocols', 'rosetta', 'settings.json')
    if arguments['--settings']:
        settings_filepath = os.path.abspath(os.path.normpath(arguments['--settings']))
    else:
        settings_filepath = default_settings_filepath
    if not os.path.exists(settings_filepath):
        die('The settings file {0} has not been created. See {0}.example for an example configuration.'.format(settings_filepath))
    try:
        settings = json.loads(read_file(settings_filepath))
    except:
        die('The settings file {0} is not a valid JSON file. See {1}.example for an example configuration.'.format(settings_filepath, default_settings_filepath))

    # Check to see whether this is being run from an SGE submission host or a workstation
    run_jobs_on_cluster = False
    rosetta_path = None
    p = subprocess.Popen(['qstat'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    if p.returncode == 0:
        run_jobs_on_cluster = arguments['--setup_only']
        rosetta_path = settings['cluster_rosetta_installation_path']
    else:
        rosetta_path = settings['local_rosetta_installation_path']
    rosetta_bin_path = os.path.join(rosetta_path, 'source', 'bin')
    rosetta_database_path = os.path.join(rosetta_path, 'database')
    if not os.path.exists(rosetta_bin_path):
        die('The Rosetta binary path was expected to be found at {0} but is missing.'.format(rosetta_bin_path))
    if not os.path.exists(rosetta_database_path):
        die('The Rosetta database path was expected to be found at {0} but is missing.'.format(rosetta_database_path))
    rosetta_binary_type = settings['rosetta_binary_type']

    # If Rosetta is built against the MySQL library, the user may need to add the MySQL libraries to the LD library path
    ld_path_extension = ''
    mysql_lib_path = arguments['--mysql_lib'] or settings.get('mysql_lib')
    if mysql_lib_path:
        mysql_lib_path = os.path.abspath(os.path.normpath(mysql_lib_path))
        if not os.path.exists(mysql_lib_path):
            if run_jobs_on_cluster:
                die('The MySQL library path "{0}" could not be found.'.format(mysql_lib_path))
            else:
                print('Warning: The MySQL library path "{0}" could not be found.'.format(mysql_lib_path))
        ld_path_extension = mysql_path_extension.format(mysql_lib_path)

    # Set up run parameters
    benchmark_parameters = dict(
        rosetta_database_path = rosetta_database_path,
        NATAA_resfile = os.path.join(benchmark_root, 'input/NATAA.res'),
        ALLAA_resfile = os.path.join(benchmark_root, 'input/ALLAA.res'),
        ntrials = int(arguments['--ntrials']),
        temperature = float(arguments['--temperature']),
        nstruct = int(arguments['--nstruct']),
        extra_flags = [],
        ld_path_extension = ld_path_extension,
        method = method,
    )
    if arguments['--talaris2014']:
        benchmark_parameters['extra_flags'].append('-talaris2014 true')

    # Set up test parameters
    if arguments['--test']:
        restricted_to_domains = ['PF00013', 'PF00018']
        benchmark_parameters['nstruct'] = 5
        benchmark_parameters['ntrials'] = 100

    # Add the flags file parameter
    if arguments['--flags']:
        try:
            flags_file = os.path.abspath(arguments['--flags'])
            assert(os.path.exists(flags_file))
            benchmark_parameters['extra_flags'].append('@{0}'.format(flags_file))
        except:
            die('The flags file {0} could not be located.'.format(flags_file))

    # Make sure that the necessary input files and binaries exist
    assert(os.path.exists(benchmark_parameters['NATAA_resfile']))
    assert(os.path.exists(benchmark_parameters['ALLAA_resfile']))
    for dependent_binary in method_details['dependent_binaries']:
        binary_path = os.path.join(rosetta_bin_path, '{0}{1}'.format(dependent_binary, rosetta_binary_type))
        if not os.path.exists(binary_path):
            die('The Rosetta binary {0} is missing and needs to be built.'.format(binary_path))
        benchmark_parameters['{0}_executable'.format(dependent_binary)] = binary_path

    # Create the output directory
    output_directory = arguments['--output_directory'] or os.path.join('.', 'job_output')
    try:
        output_directory = os.path.abspath(output_directory)
        os.makedirs(output_directory)
    except: pass
    if not os.path.exists(output_directory):
        die('The output directory {0} could not be created.'.format(output_directory))

    # Create the array jobs and submit them if appropriate
    input_pdbs = []
    input_directory = os.path.join(benchmark_root, 'input', 'domains')
    benchmark_parameters['extra_flags'] = ' '.join([s.strip() for s in benchmark_parameters['extra_flags']])
    run_script = '#!/bin/bash\n'
    for domain in os.listdir(input_directory):
        if (restricted_to_domains and (domain in restricted_to_domains)) or not(restricted_to_domains):
            for pdb in glob.glob(os.path.join(input_directory, domain, '*.pdb')):
                pdb_name = os.path.splitext(os.path.split(pdb)[1])[0]
                job_name = '%s_%s' % (domain, pdb_name)
                job_directory = os.path.join(output_directory, job_name)
                input_pdbs.append(os.path.abspath(pdb))
                if not os.path.exists(job_directory):
                    os.mkdir(job_directory)
                job_parameters = copy.deepcopy(benchmark_parameters)
                job_parameters.update(dict(
                    input_pdb = pdb_name,
                    job_directory = job_directory,
                ))
                shutil.copyfile(pdb, os.path.join(job_directory, pdb_name + '.pdb'))
                job_script_path = os.path.join(job_directory, '{0}.py'.format(method))
                script = create_script(job_script_path, job_parameters)
                run_script += ('qsub {0}\n'.format(job_script_path))
    run_script += '\n'

    if len(input_pdbs) == 0:
        die('No input PDB files were found.')
    else:
        benchmark_run_script = os.path.join(output_directory, 'run_benchmark.sh')
        write_file(benchmark_run_script, run_script)
        os.chmod(benchmark_run_script, 0755)  
        if run_jobs_on_cluster:
            print('\nSubmitting jobs:')
            subprocess.call(benchmark_run_script)
            print('Jobs for {0} domains have been submitted.\n'.format(len(input_pdbs)))
        else:
            print('\nJobs for {0} domains have been set up:\nExecute {1} to run the benchmark.\n'.format(len(input_pdbs), benchmark_run_script))




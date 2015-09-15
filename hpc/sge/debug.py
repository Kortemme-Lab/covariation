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

__doc__ = """\
This script can be used to debug covariation benchmark runs and let the user know if the expected
number of structures was created. If structures are missing, it will create new submission scripts
and print out the SGE submission commands to queue the missing jobs.

Usage:
    debug.py <output_directory> [options]

Arguments:

    <output_directory>
        The directory containing the benchmark run. The metadata file (benchmarks.json) will be used to
        determine whether or not the run succeeded.

Options:

    --nstruct NUM -n NUM
        Specify how many simulations to do for each structure in the benchmark.
        The default value is 500. [default: 500]

    --submission_script SUBMISSION_SCRIPT -s SUBMISSION_SCRIPT
        The name of the submission script. This will be used to create a new script with the suffix 
        .resume.py which can be used to resubmit missing/failed jobs.
"""

import sys
import os
import glob
import shutil
import subprocess
import shlex
import copy
import json
import re
import pprint
from libraries.utilities import *
from libraries import docopt
from libraries import colortext


def check_output(output_directory, nstruct_, submission_script):
    config = json.loads(open(os.path.join(output_directory, 'benchmarks.json')).read())
    for benchmark, methods in sorted(config.iteritems()):
        print('')
        qsub_commands = []
        for method, method_details in sorted(methods.iteritems()):
            missing_count = 0
            domains_missing_count = 0
            for d in glob.glob(os.path.join(output_directory, 'PF*')):
                nstruct = nstruct_
                header, remainder = None, None
                if submission_script:
                    script_path = os.path.join(d, submission_script)
                    assert(os.path.exists(script_path))
                    old_script = open(script_path).read()
                    mtchs = re.match('(.*)[#][$] -t 1-(\d+)(.*)', old_script, re.DOTALL)
                    header = mtchs.group(1)
                    nstruct = int(mtchs.group(2))
                    remainder = mtchs.group(3)
                    
                expected_files = range(1, nstruct + 1)
                file_filter = method_details['file_filter']              
                dirname = os.path.split(d)[1]
                domain = dirname.split('_')[0]
                pdb_id = dirname.split('_')[1]
                input_files = [(os.path.abspath(os.path.join(d, f)), re.match(file_filter, f).group('nstruct')) for f in sorted(os.listdir(d)) if re.match(file_filter, f)]
                found_files = [int(i[1]) for i in input_files]
                next_index = sorted(found_files)[-1] + 1
                missing_files = sorted(set(expected_files).difference(set(found_files)))
                if len(found_files) < nstruct: # in the first run, found_files will be a subset of expected_files but this does not hold in future runs since we number subsequent jobs starting at next_index. len(found_files) < nstruct is the more general test 
                    num_missing = len(missing_files)
                    domains_missing_count += 1
                    if next_index - 1 <= nstruct:
                        colortext.warning('Domain {0} is missing {1} files: {2}.'.format(domain, num_missing, missing_files)) # this only makes sense in the first run
                    else:
                        colortext.warning('Domain {0} is missing {1} files.'.format(domain, num_missing))
                    missing_count += num_missing
                    if submission_script:
                        F = open(os.path.join(d, submission_script + '.resume.py'), 'w')
                        F.write(header + '#$ -t {0}-{1}'.format(next_index, next_index + (num_missing*2)) + remainder) # queue up even more jobs in case the new ones fail
                        F.close
                        qsub_commands.append('qsub {0}'.format(os.path.join(d, submission_script + '.resume.py')))
            if missing_count:
                colortext.error('{0}, {1} is missing {2} files in {3} domains.\n'.format(benchmark, method, missing_count, domains_missing_count))
    if qsub_commands:
        colortext.message('qsub commands to finish the run:')
        print('\n'.join(qsub_commands))

if __name__ == '__main__':
    benchmark_root = os.path.abspath(os.path.join('..', '..'))
    arguments = docopt.docopt(__doc__)
    output_directory = arguments['<output_directory>']
    nstruct = int(arguments['--nstruct'])
    submission_script = arguments.get('--submission_script')
    check_output(output_directory, nstruct, submission_script)


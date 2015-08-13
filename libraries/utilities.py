#!/usr/bin/env python2

# The MIT License (MIT)
#
# Copyright (c) 2015 Roland Pache
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
import gzip


def die(msg, errorcode = None):
    print('\n{0}\n'.format(msg))
    sys.exit(errorcode or 1)


def read_file(filepath, binary = False):
    if binary:
        output_handle = open(filepath, 'rb')
    elif filepath.endswith('.gz'):
        output_handle = gzip.open(filepath, 'r')
    else:
        output_handle = open(filepath, 'r')
    contents = output_handle.read()
    output_handle.close()
    return contents


def write_file(filepath, contents, ftype = 'w'):
    output_handle = open(filepath, ftype)
    output_handle.write(contents)
    output_handle.close()

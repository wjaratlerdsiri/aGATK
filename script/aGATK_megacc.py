#!/usr/bin/env python
# encoding: utf-8

"""
MEGA

"""

import os,sys
import tempfile
import subprocess
import ConfigParser
import time
import argparse

from collections import namedtuple
from collections import defaultdict

import re
import glob
import shutil
import pdb

class SimpleMega():
    '''for lastz'''
    def __init__(self, target, out):
        self.output = out
        self.cli = '{2} -d {0} \
                -o {1} \
                -a data/megacc_distance.mao'.format(target, self.output, os.path.abspath("megacc"))

    def run(self):
        mega_stdout, mega_stderr = subprocess.Popen(
            self.cli,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE).communicate(None)
        return mega_stdout, mega_stderr


class Reader():
    """read a lastz file and return an iterator"""
    def __init__(self, lastz_file, long_format = False):
        self.file = open(lastz_file, 'rU')
        self.long_format = long_format

    def __del__(self):
        """close"""
        self.file.close()

    def __iter__( self ):
        """iterator"""
        while True:
            yield self.next()

class FullPaths(argparse.Action):
    """Expand paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))


def get_args():
    """Get arguments"""
    parser = argparse.ArgumentParser(
        description="""Run megacc in an easy way"""
    )
    parser.add_argument(
        "--target",
        required=True,
        type=str,
        action=FullPaths,
        help="""The path to the target file (fasta)"""
    )
    parser.add_argument(
        "--output",
        required=True,
        type=str,
        action=FullPaths,
        help="""The path to the output file"""
    )
    return parser.parse_args()


def main():
    start_time = time.time()
    print 'Started: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(start_time))
    args = get_args()
    d = SimpleMega(args.target, args.output)
    mgstdout, mgstderr = d.run()
    if mgstderr:
        raise IOError(mgstderr)
    end_time = time.time()
    print 'Ended: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(end_time))
    print 'Time for execution: ', (end_time - start_time)/60, 'minutes'

if __name__ == '__main__':
    main()

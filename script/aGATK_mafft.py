#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

class FullPaths(argparse.Action):
    """Expand paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))


def is_file(filename):
    if not os.path.isfile:
        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename


#mafft

def get_args():
    parser = argparse.ArgumentParser(
        description="""Align locally in FASTA file with MAFFT""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    #input is big fasta of many taxa and many loci
    parser.add_argument(
        "--fasta",
        required=True,
        action=FullPaths,
        type=str,
        help="""The file containing FASTA reads associated with target loci"""
    )
    parser.add_argument(
        "--output",
        required=True,
        action=FullPaths,
        type=str,
        help="""The file in which to store the resulting alignments."""
    )
    parser.add_argument(
        "--cores",
        type=int,
        default=1,
        help="""Process alignments in parallel using --cores for alignment. """ +
        """The number of PHYSICAL CPUs"""
    )
    return parser.parse_args()

class SimpleAlign():
    '''for lastz'''
    def __init__(self, fasta, cores, output):

        self.output = output
        self.cli = '{3} --localpair --thread {1} \
                --maxiterate 1000 {0} > {2}'.format(fasta, cores, self.output, os.path.abspath("mafft"))
        
    def run(self):
        mafft_stdout, mafft_stderr = subprocess.Popen(
            self.cli,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE).communicate(None)
        return mafft_stdout, mafft_stderr


def main():
    start_time = time.time()
    print 'Started: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(start_time))
    args = get_args()
    d = SimpleAlign(args.fasta, args.cores, args.output)
    mfstdout, mfstderr = d.run()
    if mfstderr:
        raise IOError(mfstderr)
    end_time = time.time()
    print 'Ended: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(end_time))
    print 'Time for execution: ', (end_time - start_time)/60, 'minutes'

if __name__ == '__main__':
    main()

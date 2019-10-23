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
#from operator import itemgetter


#from phyluce.pth import get_user_path

import pdb

class SimpleMega():
    '''docstring for lastz'''
    def __init__(self, target, out):
        # if not an output file, create a temp file to hold output
        #if not out:
            #fd, self.output = tempfile.mkstemp(suffix='.lastz')
            #os.close(fd)
        #else:
        self.output = out
        self.cli = '{2} -d {0} \
                -o {1} \
                -a /Users/w.jaratlerdsirigarvan.org.au/miniconda2/config/distance_estimation_overall_mean_coding.pdistance.pairwisedel.mao'.format(target, self.output, get_user_path("binaries", "megacc"))
#SimpleAlign() will get program binary from config file

    def run(self):
        mega_stdout, mega_stderr = subprocess.Popen(
            self.cli,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE).communicate(None)
        return mega_stdout, mega_stderr


class Reader():
    """read a lastz file and return an iterator over that file"""
    def __init__(self, lastz_file, long_format = False):
        self.file = open(lastz_file, 'rU')
        self.long_format = long_format

    def __del__(self):
        """close files"""
        self.file.close()

    def __iter__( self ):
        """iterator"""
        while True:
            yield self.next()

    def next(self):
        """read next fastq sequence and return as named tuple"""
        lastz_result = self.file.readline()
        if not lastz_result:
            raise StopIteration
        if not self.long_format:
            Lastz = namedtuple('Lastz', 'score,name1,strand1,zstart1,end1,length1,name2,'+
            'strand2,zstart2,end2,length2,diff,cigar,identity,percent_identity,'+
            'continuity,percent_continuity')
        else:
            Lastz = namedtuple('Lastz', 'score,name1,strand1,zstart1,end1,length1,name2,'+
            'strand2,zstart2,end2,length2,diff,cigar,identity,percent_identity,'+
            'continuity,percent_continuity,coverage,percent_coverage')
        aligns = defaultdict(lambda: defaultdict(list))
        lastz_result_split = lastz_result.strip('\n').split('\t')
        for k,v in enumerate(lastz_result_split):
            if k in [3,4,5,8,9,10]:
                lastz_result_split[k] = int(v)
            elif '%' in v:
                lastz_result_split[k] = float(v.strip('%'))
        lastz_result_split[1] = lastz_result_split[1].lstrip('>')
        lastz_result_split[6] = lastz_result_split[6].lstrip('>')
        return Lastz._make(lastz_result_split)

#if __name__ == '__main__':
#    print 'pass'

"""
helpers.py

Created by Brant Faircloth on 15 July 2011.
Copyright 2011 Brant C. Faircloth. All rights reserved.
"""
#import pdb
#it about PATH
class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))


#import time
#import argparse
#from phyluce import lastz
#from phyluce.helpers import FullPaths

#import pdb

def get_user_path(program, binary, package_only=False):
    config = ConfigParser.ConfigParser()
    # make case sensitive
    config.optionxform = str
    if package_only:
        config.read(os.path.join(sys.prefix, 'config/phyluce.conf'))
    else:
        config.read([os.path.join(sys.prefix, 'config/phyluce.conf'), os.path.expanduser('~/.phyluce.conf')])
    # ensure program is in list
    pth = config.get(program, binary)
    # expand path as necessary - replace CONDA variable placeholder
    # with sys.prefix, otherwise default to normal path expansion
    if pth.startswith("$CONDA"):
        expand_pth = pth.replace("$CONDA", sys.prefix)
    else:
        expand_pth = os.path.abspath(os.path.expanduser(os.path.expandvars(pth)))
    return expand_pth


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
        description="""Run megacc in an easy way"""
    )
    parser.add_argument(
        "--target",
        required=True,
        type=str,
        action=FullPaths,
        help="""The path to the target file (2bit/fasta)"""
    )
    #parser.add_argument(
        #"--query",
        #required=True,
        #type=str,
        #action=FullPaths,
        #help="""The path to the query file (2bit/fasta)"""
    #)
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

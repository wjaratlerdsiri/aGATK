#!/usr/bin/env python
# encoding: utf-8

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

class SimpleAlign():
    '''for lastz'''
    def __init__(self, target, query, out=False):
        if not out:
            fd, self.output = tempfile.mkstemp(suffix='.lastz')
            os.close(fd)
        else:
            self.output = out
        self.cli = '{3} {0}[multiple,nameparse=full] {1}[nameparse=full] \
                --output={2} --ambiguous=iupac \
                --format=general-:score,name1,strand1,zstart1,end1,length1,name2,strand2,zstart2,end2,length2,diff,cigar,identity,continuity'.format(target, query, self.output, os.path.abspath("lastz"))

    def run(self):
        lastz_stdout, lastz_stderr = subprocess.Popen(
            self.cli,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE).communicate(None)
        return lastz_stdout, lastz_stderr

class Align():
    '''for lastz'''
    def __init__(self, target, query, coverage, identity, out=False, min_match=None):
        if not out:
            fd, self.output = tempfile.mkstemp(suffix='.lastz')
            os.close(fd)
        else:
            self.output = out
        if identity and not min_match:
            self.cli = '{5} {0}[multiple,nameparse=full] {1}[nameparse=full] \
                --strand=both \
                --seed=12of19 \
                --transition \
                --nogfextend \
                --nochain \
                --gap=400,30 \
                --xdrop=910 \
                --ydrop=8370 \
                --hspthresh=3000 \
                --gappedthresh=3000 \
                --noentropy \
                --coverage={2} \
                --identity={3} \
                --output={4} --ambiguous=iupac \
                --format=general-:score,name1,strand1,zstart1,end1,length1,name2,strand2,zstart2,end2,length2,diff,cigar,identity,continuity'.format(target, query, coverage, identity, self.output, os.path.abspath("lastz"))
        elif min_match:
            self.cli = '{5} {0}[multiple,nameparse=full] {1}[nameparse=full] \
                --strand=both \
                --seed=12of19 \
                --transition \
                --nogfextend \
                --nochain \
                --gap=400,30 \
                --xdrop=910 \
                --ydrop=8370 \
                --hspthresh=3000 \
                --gappedthresh=3000 \
                --noentropy \
                --matchcount={2} \
                --identity={3} \
                --output={4} --ambiguous=iupac \
                --format=general-:score,name1,strand1,zstart1,end1,length1,name2,strand2,zstart2,end2,length2,diff,cigar,identity,continuity'.format(target, query, min_match, identity, self.output, os.path.abspath("lastz"))

    def run(self):
        lastz_stdout, lastz_stderr = subprocess.Popen(
            self.cli,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE).communicate(None)
        return lastz_stdout, lastz_stderr


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

    def next(self):
        """read next fasta sequence"""
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


class FullPaths(argparse.Action):
    """Expand paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

def get_args():
    """Get arguments"""
    parser = argparse.ArgumentParser(
        description="""Run lastz in an easy way"""
    )
    parser.add_argument(
        "--target",
        required=True,
        type=str,
        action=FullPaths,
        help="""The path to the target file (fasta)"""
    )
    parser.add_argument(
        "--query",
        required=True,
        type=str,
        action=FullPaths,
        help="""The path to the query file (fasta)"""
    )
    parser.add_argument(
        "--output",
        required=True,
        type=str,
        action=FullPaths,
        help="""The path to the output file"""
    )
    parser.add_argument(
        "--identity",
        type=float,
        default=92.5,
        help="""The minimum percent identity to require for a match"""
    )
    cov_or_match = parser.add_mutually_exclusive_group(required=False)
    cov_or_match.add_argument(
        "--coverage",
        type=float,
        default=83.0,
        help="""The minimum coverage (%%) required for a match"""
    )
    cov_or_match.add_argument(
        "--min_match",
        type=int,
        default=None,
        help="""The minimum number of base pairs required for a match"""
    )
    return parser.parse_args()


def main():
    start_time = time.time()
    print 'Started: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(start_time))
    args = get_args()
    alignment = Align(args.target, args.query, args.coverage, args.identity, args.output, args.min_match)
    lzstdout, lztstderr = alignment.run()
    if lztstderr:
        raise IOError(lztstderr)
    end_time = time.time()
    print 'Ended: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(end_time))
    print 'Time for execution: ', (end_time - start_time)/60, 'minutes'

if __name__ == '__main__':
    main()



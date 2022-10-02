#! /usr/bin/env python3
'''
Convert a CRAM/BAM/SAM file to a consensus sequence
'''

# imports
from datetime import datetime
from os.path import isdir, isfile
from sys import argv, stderr
import argparse
import pysam

# constants
DEFAULT_BUFSIZE = 1048576 # 1 MB

# print log
def print_log(s='', end='\n'):
    tmp = "[%s] %s" % (datetime.now().strftime("%Y-%m-%d %H:%M:%S"), s)
    print(tmp, file=stderr); stderr.flush()

# error message
def error(s=None):
    if s is None:
        print_log("ERROR")
    else:
        print_log("ERROR: %s" % s)
    exit(1)

# parse user args
def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=False, type=str, default='stdin', help="Input File (CRAM/BAM/SAM)")
    parser.add_argument('-r', '--reference', required=True, type=str, help="Reference Genome (FASTA)")
    parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Output File (BSON)")
    args = parser.parse_args()
    if args.input.lower() != 'stdin' and not isfile(args.input):
        error("Input file not found: %s" % args.input)
    if args.output.lower() != 'stdout' and (isfile(args.output) or isdir(args.output)):
        error("Output file already exits: %s" % args.output)
    return args

# compute length of reference genome
def compute_ref_length(fn, bufsize=DEFAULT_BUFSIZE):
    return sum(len(l.strip()) for l in open(fn, buffering=bufsize) if not l.startswith('>'))

# open input CRAM file
def open_aln(fn):
    if fn.lower() == 'stdin':
        fn = '-'
    tmp = pysam.set_verbosity(0) # disable htslib verbosity to avoid "no index file" warning
    aln = pysam.AlignmentFile(fn, mode='r', reference_filename=args.reference)
    pysam.set_verbosity(tmp) # re-enable htslib verbosity and finish up
    return aln

# main content
if __name__ == "__main__":
    # prep user input
    args = parse_args()
    ref_len = compute_ref_length(args.reference)
    aln = open_aln(args.input)
    count_mapped = 0 # TODO DELETE
    count_unmapped = 0 # TODO DELETE MAYBE
    for read_num, read in enumerate(aln.fetch(until_eof=True)):
        if read.is_unmapped:
            count_unmapped += 1 # TODO DELETE MAYBE
        else:
            count_mapped += 1 # TODO DELETE
    print("Mapped:   %d" % count_mapped) # TODO DELETE
    print("Unmapped: %d" % count_unmapped) # TODO DELETE MAYBE

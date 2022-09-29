#! /usr/bin/env python3
'''
Convert a BAM file to reference-guided BSON
'''

# imports
from datetime import datetime
from os.path import isdir, isfile
from sys import argv, stderr
import argparse
import pysam

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
    parser.add_argument('-i', '--input', required=False, type=str, default='stdin', help="Input File (SAM/BAM)")
    parser.add_argument('-r', '--reference', required=True, type=str, help="Reference Genome (FASTA)")
    parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Output File (BSON)")
    parser.add_argument('--force_bam', action="store_true", help="Force BAM Input (otherwise infer from filename)")
    args = parser.parse_args()
    if args.output.lower() != 'stdout' and (isfile(args.output) or isdir(args.output)):
        error("Output already exits: %s" % args.output)
    return args

# open input SAM/BAM file
def open_sam(fn, force_bam=False):
    tmp = pysam.set_verbosity(0) # disable htslib verbosity to avoid "no index file" warning
    if fn.lower() == 'stdin':
        if force_bam:
            aln = pysam.AlignmentFile('-', 'rb')
        else:
            aln = pysam.AlignmentFile('-', 'r') # standard input default --> SAM
    elif not isfile(fn):
        error("File not found: %s" % fn)
    elif force_bam or fn.lower().endswith('.bam'):
        aln = pysam.AlignmentFile(fn, 'rb')
    elif fn.lower().endswith('.sam'):
        aln = pysam.AlignmentFile(fn, 'r')
    else:
        error("Invalid input alignment file extension: %s" % fn)
    pysam.set_verbosity(tmp) # re-enable htslib verbosity and finish up
    return aln

# main content
if __name__ == "__main__":
    # prep user input
    args = parse_args()
    aln = open_sam(args.input, args.force_bam)
    count_mapped = 0 # TODO DELETE
    count_unmapped = 0 # TODO DELETE MAYBE
    for read_num, read in enumerate(aln.fetch(until_eof=True)):
        if read.is_unmapped:
            count_unmapped += 1 # TODO DELETE MAYBE
        else:
            count_mapped += 1 # TODO DELETE
    print("Mapped:   %d" % count_mapped) # TODO DELETE
    print("Unmapped: %d" % count_unmapped) # TODO DELETE MAYBE

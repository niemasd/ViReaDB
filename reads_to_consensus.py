#! /usr/bin/env python3
'''
Convert a CRAM/BAM/SAM file to a consensus sequence
'''

# imports
from datetime import datetime
from numpy import uintc, zeros
from os.path import isdir, isfile
from sys import argv, stderr
import argparse
import pysam

# constants
DEFAULT_BUFSIZE = 1048576 # 1 MB
DEFAULT_MIN_QUAL = 20
DEFAULT_MIN_DEPTH = 10
DEFAULT_MIN_FREQ = 0.5
DEFAULT_AMBIG = 'N'
BASE_TO_NUM = {'A':0, 'C':1, 'G':2, 'T':3, None:4}

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
    parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Output File (FASTA)")
    parser.add_argument('-q', '--min_qual', required=False, type=int, default=20, help="Minimum Quality to Count Base")
    parser.add_argument('-m', '--min_depth', required=False, type=int, default=10, help="Minimum Depth to Call Consensus Base")
    parser.add_argument('-t', '--min_freq', required=False, type=float, default=0.5, help="Minimum Frequency to Call Consensus Base")
    parser.add_argument('-n', '--ambig', required=False, type=str, default='N', help="Ambiguous Symbol")
    args = parser.parse_args()
    if args.input.lower() != 'stdin' and not isfile(args.input):
        error("Input file not found: %s" % args.input)
    if args.output.lower() != 'stdout' and (isfile(args.output) or isdir(args.output)):
        error("Output file already exits: %s" % args.output)
    if args.min_qual < 0:
        error("Minimum Quality must be non-negative: %s" % args.min_qual)
    if args.min_depth < 1:
        error("Minimum Depth must be positive: %s" % args.min_depth)
    if args.min_freq < 0 or args.min_freq > 1:
        error("Minimum Frequency must be in range [0,1]: %s" % args.min_freq)
    if len(args.ambig) != 1:
        error("Ambiguous Symbol must be 1 character in length: %s" % args.ambig)
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

# compute base counts
def compute_base_counts(aln, ref_len, min_qual=DEFAULT_MIN_QUAL, min_depth=DEFAULT_MIN_DEPTH, min_freq=DEFAULT_MIN_FREQ, ambig=DEFAULT_AMBIG):
    pos_counts = zeros((ref_len,len(BASE_TO_NUM)), dtype=uintc)
    ins_counts = dict()
    for read_num, read in enumerate(aln.fetch(until_eof=True)):
        # skip unmapped reads
        if read.is_unmapped:
            continue

        # prepare some helper variables
        seq = read.query_sequence; quals = read.query_qualities
        aln_start = read.query_alignment_start; aln_end = read.query_alignment_end
        aligned_pairs = tuple(read.get_aligned_pairs())
        deletion_start_inds = list(); deletion_end_inds = list()

        # increment base counts and keep track of deletions
        for pair_ind, pair in enumerate(aligned_pairs):
            read_pos, ref_pos = pair
            if read_pos is None:
                if len(deletion_end_inds) != 0 and deletion_end_inds[-1] == pair_ind-1:
                    deletion_end_inds[-1] = pair_ind
                else:
                    deletion_start_inds.append(pair_ind); deletion_end_inds.append(pair_ind)
            elif read_pos < aln_start or read_pos >= aln_end or quals[read_pos] < args.min_qual:
                continue # skip soft-clipped or low-quality bases
            elif ref_pos is None: # insertion
                pass # TODO HANDLE INSERTIONS
            else:
                pos_counts[ref_pos][BASE_TO_NUM[seq[read_pos]]] += 1

        # handle deletions
        if len(deletion_end_inds) != 0:
            for deletion_start_ind, deletion_end_ind in zip(deletion_start_inds, deletion_end_inds):
                if deletion_start_ind == 0 or deletion_end_ind == len(aligned_pairs)-1:
                    continue # this deletion is at the very beginning or end of the read, so skip
                prev_read_pos = aligned_pairs[deletion_start_ind-1][0]
                if prev_read_pos < aln_start:
                    continue # this deletion is right after a soft-clipped start, so skip
                next_read_pos = aligned_pairs[deletion_end_ind+1][0]
                if next_read_pos >= aln_end:
                    continue # this deletion is right before a soft-clipped end, so skip
                if quals[prev_read_pos] < args.min_qual or quals[next_read_pos] < args.min_qual:
                    continue # base just before or just after this deletion is low quality, so skip
                for deletion_ind in range(deletion_start_ind, deletion_end_ind+1):
                    pos_counts[aligned_pairs[deletion_ind][1]][4] += 1 # '-' = 4
    return pos_counts, ins_counts

# main content
if __name__ == "__main__":
    args = parse_args()
    ref_len = compute_ref_length(args.reference)
    aln = open_aln(args.input)
    pos_counts, ins_counts = compute_base_counts(aln, ref_len, args.min_qual, args.min_depth, args.min_freq, args.ambig)

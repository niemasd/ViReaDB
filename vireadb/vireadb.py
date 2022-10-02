#! /usr/bin/env python3
'''
vireadb: Viral Read Database
'''

# imports
from .common import *
from numpy import uintc, zeros
from os import remove
from os.path import isdir, isfile
from sys import argv
import argparse
import pysam
import sqlite3

# constants
META_TABLE_COLS = ('key', 'val')
READ_TABLE_COLS = ('ID', 'CRAM', 'NUC_COUNTS', 'INS_COUNTS', 'CONSENSUS')
DEFAULT_BUFSIZE = 1048576 # 1 MB
DEFAULT_MIN_QUAL = 20
DEFAULT_MIN_DEPTH = 10
DEFAULT_MIN_FREQ = 0.5
DEFAULT_AMBIG = 'N'
BASE_TO_NUM = {'A':0, 'C':1, 'G':2, 'T':3, None:4}
NUM_TO_BASE = 'ACGT-'

class ViReaDB:
    '''``ViReaDB`` database class'''
    def __init__(self, db_fn):
        '''``ViReaDB`` constructor

        Args:
            ``db_fn`` (``str``): The filename of the SQLite3 database file representing this database

        Returns:
            ``ViReaDB`` object
        '''
        self.con = sqlite3.connect(db_fn)
        self.cur = self.con.cursor()

def create_db(db_fn, overwrite=False):
    '''Create a new ViReaDB database

    Args:
        ``db_fn`` (``str``): The filename of the SQLite3 database file representing this database

        ``overwrite`` (``bool``): Overwrite ``db_fn`` if it already exists

    Returns:
        ``ViReaDB`` object
    '''
    # check valid database filename
    if isdir(db_fn):
        error("db_fn exists as a directory: %s" % db_fn)
    if isfile(db_fn):
        if overwrite:
            remove(db_fn)
        else:
            error("db_fn exists: %s" % db_fn)

    # create SQLite3 database and populate with `meta` table
    db = ViReaDB(db_fn)
    db.cur.execute("CREATE TABLE meta(%s)" % ', '.join(META_TABLE_COLS))
    db.cur.execute("""
        INSERT INTO meta VALUES
            ('VERSION', '%s')
    """ % VERSION)
    db.con.commit()
    return db

def load_db(db_fn):
    '''Load a ViReaDB database from file

    Args:
        ``db_fn`` (``str``): The filename of the SQLite3 database file representing this database

    Returns:
        ``ViReaDB`` object
    '''
    if not isfile(db_fn):
        error("db_fn not found: %s" % db_fn)
    return ViReaDB(db_fn)

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
def compute_base_counts(aln, ref_len, min_qual=DEFAULT_MIN_QUAL):
    # prepare some global helper variables
    pos_counts = zeros((ref_len,len(BASE_TO_NUM)), dtype=uintc)
    ins_counts = dict() # ins_counts[ref_pos] stores the counts of all insertions that occurred *before* position ref_pos
    insertion_start_inds = list(); insertion_end_inds = list()
    deletion_start_inds = list(); deletion_end_inds = list()

    # iterate over all reads in CRAM/BAM/SAM
    for read_num, read in enumerate(aln.fetch(until_eof=True)):
        # skip unmapped reads
        if read.is_unmapped:
            continue

        # prepare some within-loop helper variables
        seq = read.query_sequence; quals = read.query_qualities
        aln_start = read.query_alignment_start; aln_end = read.query_alignment_end
        aligned_pairs = tuple(read.get_aligned_pairs())

        # increment base counts and keep track of deletions
        for pair_ind, pair in enumerate(aligned_pairs):
            read_pos, ref_pos = pair
            if read_pos is None: # deletion
                if len(deletion_end_inds) != 0 and deletion_end_inds[-1] == pair_ind-1:
                    deletion_end_inds[-1] = pair_ind
                else:
                    deletion_start_inds.append(pair_ind); deletion_end_inds.append(pair_ind)
            elif read_pos < aln_start or read_pos >= aln_end or quals[read_pos] < args.min_qual:
                continue # skip soft-clipped or low-quality bases
            elif ref_pos is None: # insertion
                if len(insertion_end_inds) != 0 and insertion_end_inds[-1] == pair_ind-1:
                    insertion_end_inds[-1] = pair_ind
                else:
                    insertion_start_inds.append(pair_ind); insertion_end_inds.append(pair_ind)
            else: # (mis)matches
                pos_counts[ref_pos][BASE_TO_NUM[seq[read_pos]]] += 1

        # handle insertions (if any)
        if len(insertion_end_inds) != 0:
            for insertion_start_ind, insertion_end_ind in zip(insertion_start_inds, insertion_end_inds):
                if insertion_start_ind == 0 or insertion_end_ind == len(aligned_pairs)-1:
                    continue # this insertion is at the very beginning or end of the read, so skip
                prev_read_pos = aligned_pairs[insertion_start_ind-1][0]
                if prev_read_pos < aln_start:
                    continue # this insertion is right after a soft-clipped start, so skip
                next_read_pos, next_ref_pos = aligned_pairs[insertion_end_ind+1]
                if next_read_pos >= aln_end:
                    continue # this insertion is right before a soft-clipped end, so skip
                if quals[prev_read_pos] < args.min_qual or quals[next_read_pos] < args.min_qual:
                    continue # base just before or just after this insertion is low quality, so skip
                ins_seq = seq[aligned_pairs[insertion_start_ind][0] : aligned_pairs[insertion_end_ind][0]+1]
                if next_ref_pos not in ins_counts:
                    ins_counts[next_ref_pos] = dict()
                if ins_seq not in ins_counts[next_ref_pos]:
                    ins_counts[next_ref_pos][ins_seq] = 0
                ins_counts[next_ref_pos][ins_seq] += 1
            insertion_start_inds.clear(); insertion_end_inds.clear()


        # handle deletions (if any)
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
            deletion_start_inds.clear(); deletion_end_inds.clear()
    return pos_counts, ins_counts

# generate consensus sequence
def generate_consensus(pos_counts, ins_counts, min_depth=DEFAULT_MIN_DEPTH, min_freq=DEFAULT_MIN_FREQ, ambig=DEFAULT_AMBIG):
    parts = ['']*(len(pos_counts)+len(ins_counts)); ind = 0
    pos_count_tots = [float(sum(row)) for row in pos_counts]
    for ref_pos in range(len(pos_counts)+1):
        # handle insertions before ref_pos
        if ref_pos in ins_counts:
            curr_ins_counts = ins_counts[ref_pos]
            ins_seqs = sorted(((curr_ins_counts[s], s) for s in curr_ins_counts), reverse=True)
            best_c, best_s = ins_seqs[0]; tot = sum(c for c,s in ins_seqs)
            if ref_pos != 0:
                tot += pos_count_tots[ref_pos-1]
            if ref_pos != len(pos_counts):
                tot += pos_count_tots[ref_pos+1]
            if tot >= min_depth:
                freq = best_c / tot
                if freq >= min_freq:
                    parts[ind] = best_s; ind += 1

        # handle ref_pos
        if ref_pos != len(pos_counts):
            curr_counts = sorted(((c,i) for i,c in enumerate(pos_counts[ref_pos])), reverse=True)
            best_c, best_i = curr_counts[0]
            tot = pos_count_tots[ref_pos]
            if tot < min_depth:
                parts[ind] = ambig; ind += 1; continue
            freq = best_c / tot
            if freq < min_freq:
                parts[ind] = ambig; ind += 1; continue
            parts[ind] = NUM_TO_BASE[best_i]; ind += 1
    return ''.join(parts)

# main content
if __name__ == "__main__":
    args = parse_args()
    ref_len = compute_ref_length(args.reference)
    aln = open_aln(args.input)
    pos_counts, ins_counts = compute_base_counts(aln, ref_len, args.min_qual)
    cons_seq = generate_consensus(pos_counts, ins_counts, args.min_depth, args.min_freq, args.ambig)
    if args.output.lower() == 'stdout':
        from sys import stdout as out_f
    else:
        out_f = open(args.output, 'w')
    out_f.write(">%s\n%s\n" % (' '.join(argv), cons_seq))
    out_f.close()

#! /usr/bin/env python3
'''
Functions dealing with FASTA files
'''

# import statements
from .common import *
from os.path import isfile

def load_ref(ref_fn, bufsize=DEFAULT_BUFSIZE):
    '''Load the reference genome FASTA

    Args:
        ``ref_fn`` (``str``): The reference genome FASTA file

    Returns:
        A ``str`` containing the contents of the one-line FASTA file of the reference genome (i.e., sequence whitespace removed)
    '''
    if not isfile(ref_fn):
        raise ValueError("File not found: %s" % ref_fn)
    f = open(ref_fn, buffering=bufsize); lines = f.read().strip().splitlines(); f.close()
    header = None; seq_parts = list()
    for line_num, line in enumerate(lines):
        if line.startswith('>'):
            if line_num != 0:
                raise ValueError("Reference FASTA must have a single sequence")
            header = line.strip()
        else:
            seq_parts.append(line.strip())
    return header, ''.join(seq_parts)

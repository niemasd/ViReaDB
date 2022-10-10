#! /usr/bin/env python3
'''
vireadb: Viral Read Database
'''

# imports
from .common import *
from .cram import *
from .fasta import *
from os import remove
from os.path import isdir, isfile
from subprocess import check_output
from sys import argv
from tempfile import NamedTemporaryFile
from warnings import warn
import argparse
import sqlite3

# constants
META_TABLE_COLS = ('key', 'val')
SEQS_TABLE_COLS = ('ID', 'CRAM', 'POS_COUNTS', 'INS_COUNTS', 'CONSENSUS')

class ViReaDB:
    '''``ViReaDB`` database class'''
    def __init__(self, db_fn, bufsize=DEFAULT_BUFSIZE):
        '''``ViReaDB`` constructor

        Args:
            ``db_fn`` (``str``): The filename of the SQLite3 database file representing this database

        Returns:
            ``ViReaDB`` object
        '''
        self.con = sqlite3.connect(db_fn)
        self.cur = self.con.cursor()
        self.version = self.cur.execute("SELECT val FROM meta WHERE key='VERSION'").fetchone()[0]
        self.ref_name = self.cur.execute("SELECT val FROM meta WHERE key='REF_NAME'").fetchone()[0]
        self.ref_seq = self.cur.execute("SELECT val FROM meta WHERE key='REF_SEQ'").fetchone()[0]
        self.ref_len = len(self.ref_seq)
        self.ref_f = NamedTemporaryFile('w', prefix='vireadb', suffix='.fas', buffering=bufsize)
        self.ref_f.write('%s\n%s\n' % (self.ref_name, self.ref_seq)); self.ref_f.flush()

    def __del__(self):
        '''``ViReaDB`` destructor'''
        self.con.close()

    def commit(self):
        '''Commit the SQLite3 database'''
        self.con.commit()

    def add_entry(self, ID, reads_fn, filetype=None, include_unmapped=False, bufsize=DEFAULT_BUFSIZE, threads=DEFAULT_THREADS, commit=True):
        '''Add a CRAM/BAM/SAM entry to this database

        Args:
            ``ID`` (``str``): The unique ID of the entry to add

            ``reads_fn`` (``str``): The input reads file

            ``filetype`` (``str``): The format of the input reads file (CRAM, BAM, or SAM), or None to infer from ``reads_fn``

            ``bufsize`` (``int``): Buffer size for reading from file
            
            ``threads`` (``int``): Number of threads to use for compression

            ``commit`` (``bool``): Commit database after adding this entry
        '''
        # check for validity
        if self.cur.execute("SELECT ID FROM seqs WHERE ID='%s'" % ID).fetchone() is not None:
            raise ValueError("ID already exists in database: %s" % ID)
        if not isfile(reads_fn):
            raise ValueError("File not found: %s" % reads_fn)
        if filetype is None:
            filetype = reads_fn.split('.')[-1].upper()
        if not isinstance(filetype, str):
            raise TypeError("Invalid filetype: %s (must be CRAM, BAM, or SAM)" % filetype)
        filetype = filetype.strip().upper()
        if not isinstance(threads, int) or threads < 1:
            raise ValueError("Invalid number of threads: %s" % threads)

        # handle CRAM (just read all data)
        if filetype == 'CRAM':
            f = open(reads_fn, 'rb', buffering=bufsize); cram_data = f.read(); f.close()

        # handle BAM/SAM (convert to CRAM)
        elif filetype == 'BAM' or filetype == 'SAM':
            command = [
                'samtools', 'view',
                '--output-fmt-option', 'version=3.0', # TODO update to 3.1 when stable
                '--output-fmt-option', 'use_lzma=1',
                '--output-fmt-option', 'archive=1',
                '--output-fmt-option', 'level=9',
                '--output-fmt-option', 'lossy_names=1',
                '-@', str(threads),
                '-T', self.ref_f.name,
                '-C', # CRAM output
            ]
            if not include_unmapped:
                command += ['-F', '4'] # only include mapped reads
            command += [reads_fn]
            cram_data = check_output(command)

        # invalid filetype
        else:
            raise TypeError("Invalid filetype: %s (must be CRAM, BAM, or SAM)" % filetype)

        # add this CRAM to the database
        curr_row = (ID, cram_data, None, None, None)
        self.cur.execute("INSERT INTO seqs VALUES(?, ?, ?, ?, ?)", curr_row)
        if commit:
            self.commit()

    def del_entry(self, ID, commit=True):
        '''Remove an entry to this database

        Args:
            ``ID`` (``str``): The unique ID of the entry to remove

            ``commit`` (``bool``): Commit database after removing this entry
        '''
        if self.cur.execute("SELECT ID FROM seqs WHERE ID='%s'" % ID).fetchone() is None:
            raise ValueError("ID doesn't exist in database: %s" % ID)
        self.cur.execute("DELETE FROM seqs WHERE ID='%s'" % ID)
        if commit:
            self.commit()

    def compute_counts(self, ID, overwrite=False, commit=True):
        '''Compute position and insertion counts for a given entry

        Args:
            ``ID`` (``str``): The unique ID of the entry whose counts to compute

            ``overwrite`` (``bool``): ``True`` to recompute (and overwrite) counts if they already exist

            ``commit`` (``bool``): Commit database after updating this entry
        '''
        tmp = self.cur.execute("SELECT POS_COUNTS, INS_COUNTS FROM seqs WHERE ID='%s'" % ID).fetchone()
        if tmp is None:
            raise ValueError("ID doesn't exist in database: %s" % ID)
        pos_counts, ins_counts = tmp
        raise RuntimeError("TODO HERE")


def create_db(db_fn, ref_fn, overwrite=False):
    '''Create a new ViReaDB database

    Args:
        ``db_fn`` (``str``): The filename of the SQLite3 database file representing this database

        ``ref_fn`` (``str``): The filename of the viral reference genome to use for this database

        ``overwrite`` (``bool``): Overwrite ``db_fn`` if it already exists

    Returns:
        ``ViReaDB`` object
    '''
    # check valid inputs
    if not isfile(ref_fn):
        raise ValueError("File not found: %s" % ref_fn)
    if isdir(db_fn):
        raise ValueError("db_fn exists as a directory: %s" % db_fn)
    if isfile(db_fn):
        if overwrite:
            remove(db_fn)
        else:
            raise ValueError("db_fn exists: %s" % db_fn)

    # load reference genome
    ref_name, ref_seq = load_ref(ref_fn)

    # create SQLite3 database and populate with `meta` table
    con = sqlite3.connect(db_fn); cur = con.cursor()
    cur.execute("CREATE TABLE meta(%s)" % ', '.join(META_TABLE_COLS))
    cur.execute("INSERT INTO meta VALUES(?, ?)", ('VERSION', VERSION))
    cur.execute("INSERT INTO meta VALUES(?, ?)", ('REF_NAME', ref_name))
    cur.execute("INSERT INTO meta VALUES(?, ?)", ('REF_SEQ', ref_seq))
    cur.execute("CREATE TABLE seqs(%s)" % ', '.join(SEQS_TABLE_COLS))
    con.commit(); con.close()
    return ViReaDB(db_fn)

def load_db(db_fn):
    '''Load a ViReaDB database from file

    Args:
        ``db_fn`` (``str``): The filename of the SQLite3 database file representing this database

    Returns:
        ``ViReaDB`` object
    '''
    if not isfile(db_fn):
        raise ValueError("db_fn not found: %s" % db_fn)
    return ViReaDB(db_fn)

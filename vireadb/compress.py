#! /usr/bin/env python3
'''
Functions dealing with file compression
'''

# import statements
from .common import *
from io import BytesIO
from lzma import LZMACompressor, LZMADecompressor, PRESET_EXTREME
import json

def compress_str(s):
    '''Compress a string

    Args:
        ``s`` (``str``): The string to compress

    Returns:
        ``bytes`` object containing the LZMA-compressed string
    '''
    if s is None:
        return None
    lzma_comp = LZMACompressor(preset=PRESET_EXTREME)
    return lzma_comp.compress(s.encode("ascii")) + lzma_comp.flush()

def decompress_str(s_xz):
    '''Decompress a compressed string

    Args:
        ``s_xz`` (``bytes``): The LZMA-compressed string to decompress

    Returns:
        ``str`` object containing the decompressed string
    '''
    lzma_decomp = LZMADecompressor()
    return lzma_decomp.decompress(s_xz).decode()

def compress_pos_counts(pos_counts):
    '''Compress position counts

    Args:
        ``pos_counts`` (``numpy.array``): Position counts

    Returns:
        ``bytes`` object containing the LZMA-compressed position counts
    '''
    if pos_counts is None:
        return None
    lzma_comp = LZMACompressor(preset=PRESET_EXTREME)
    pos_counts_vf = BytesIO(); numpy.save(pos_counts_vf, pos_counts, allow_pickle=False); pos_counts_vf.seek(0)
    return lzma_comp.compress(pos_counts_vf.read()) + lzma_comp.flush()

def decompress_pos_counts(pos_counts_xz):
    '''Decompress position counts

    Args:
        ``pos_counts_xz`` (``bytes``): The LZMA-compressed position counts to decompress

    Returns:
        ``numpy.array`` object containing position counts
    '''
    if pos_counts_xz is None:
        return None
    lzma_decomp = LZMADecompressor()
    return numpy.load(BytesIO(lzma_decomp.decompress(pos_counts_xz)), allow_pickle=False)

def compress_ins_counts(ins_counts):
    '''Compress insertion counts

    Args:
        ``ins_counts`` (``dict``): Insertion counts

    Returns:
        ``bytes`` object containing the LZMA-compressed insertion counts
    '''
    if ins_counts is None:
        return None
    lzma_comp = LZMACompressor(preset=PRESET_EXTREME)
    return lzma_comp.compress(json.dumps(ins_counts).encode("ascii")) + lzma_comp.flush()

def decompress_ins_counts(ins_counts_xz):
    '''Decompress insertion counts

    Args:
        ``ins_counts_xz`` (``bytes``): The LZMA-compressed insertion counts to decompress

    Returns:
        ``dict`` object containing insertion counts
    '''
    if ins_counts_xz is None:
        return None
    lzma_decomp = LZMADecompressor()
    return json.loads(lzma_decomp.decompress(ins_counts_xz).decode("ascii"))

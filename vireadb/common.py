#! /usr/bin/env python3
'''
Common Functions
'''

# imports
from datetime import datetime
from sys import stderr

# constants
VERSION = '0.0.1'
DEFAULT_BUFSIZE = 1048576 # 1 MB

def print_log(s='', end='\n'):
    '''Print log message to standard error

    Args:
        ``s`` (``str``): The log message to print

        ``end`` (``str``): The symbol to print after printing ``s``
    '''
    tmp = "[%s] %s" % (datetime.now().strftime("%Y-%m-%d %H:%M:%S"), s)
    print(tmp, file=stderr); stderr.flush()

def error(s=None, exit_status=None):
    '''Throw an error

    Args:
        ``s`` (``str``): The error message

        ``exit_status`` (``int``): The exit status with which to kill the program, otherwise ``None`` to not kill
    '''
    if s is None:
        message = "ERROR"
    else:
        message = "ERROR: %s" % s
    if exit_status is None:
        raise RuntimeError(message)
    else:
        print_log(message); exit(exit_status)

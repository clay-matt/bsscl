################################

# Matt Clay
# version 130725

################################

from scl import *
import argparse # parse command line arguments

# define command-line parser
parser = argparse.ArgumentParser(prog='bsscl.py',
                                description='Compute lower bound for scl(g) in BS(m,l)',
                                epilog='Created by: Matt Clay, email:mattclay@uark.edu')
parser.add_argument('word', metavar='g')
parser.add_argument('m', type=int)
parser.add_argument('l', type=int)
parser.add_argument('-v','--verbose', action='count', help='increase output verbosity')

# process command-line arguments
args = parser.parse_args()

# set variables from parser
g = args.word
m = args.m
l = args.l
verbose = args.verbose

# compute scl
scl_g = scl(g,m,l,verbose)

# if verbose == True value is already displayed to screen
if not verbose: print scl_g

################################

# Matt Clay
# version 130507

################################

import argparse # parse command line arguments
from scl import *

# define command-line parser
parser = argparse.ArgumentParser(prog='bsscl.sage',
                                 description='Compute lower and upper bound for scl(g) in BS(m,l)',
                                 epilog='Created by: Matt Clay, email: mattclay@uark.edu')
parser.add_argument('word', metavar='g')
parser.add_argument('m', type=int)
parser.add_argument('l', type=int)
parser.add_argument('-b','--bound', type=int, choices=[0,1], metavar='B', default=0,
                    help='bound toggle: 0 = lower bound (default), 1 = upper bound')
parser.add_argument('-v','--verbose', action='store_true', help='increase output verbosity')

# process command-line arguments
args = parser.parse_args()

# set variables from parser
g = args.word
m = args.m
l = args.l
bound = args.bound
verbose = args.verbose

# compute scl
scl_g = scl(g,m,l,bound,verbose)

# if verbose == True value is already displayed to screen
if not verbose: print scl_g

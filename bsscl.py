#*****************************************************************************
#       Copyright (C) 2013 Matt Clay <mattclay@uark.edu>
# 
#  Distributed under the terms of the GNU General Public License (GPL) 
#                  http://www.gnu.org/licenses/ 
#***************************************************************************** 

from utils import *
from scl import *
import argparse # parse command line arguments

# define command-line parser
parser = argparse.ArgumentParser(prog='bsscl.py',
                                description='Compute lower bound for scl(g) in BS(m,l)',
                                epilog='Created by: Matt Clay, email:mattclay@uark.edu')
parser.add_argument('g', metavar='g')
parser.add_argument('m', type=int)
parser.add_argument('l', type=int)
parser.add_argument('-v','--verbose', action='count', help='increase output verbosity')

# process command-line arguments
args = parser.parse_args()

# set variables from parser
g = args.g
m = args.m
l = args.l
verbose = args.verbose

# compute scl
scl_g = scl(g,m,l,verbose)

# if verbose > 0 value is already displayed to screen
if verbose == 0: print scl_g

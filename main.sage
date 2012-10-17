################################

# Matt Clay
# version 121017

################################

import sys # get command line arguments
# 1 = word, 2 = m, 3 = n, 4 = verbose

from scl import *

# parse input
if len(sys.argv) < 4:
    print 'Invalid input: not enough arguments'
    quit
    
g = sys.argv[1]
m = int(sys.argv[2])
n = int(sys.argv[3])
verbose = True if len(sys.argv)==5 else False

# compute scl
scl_g = scl(g,m,n,verbose)

if verbose == False: print scl_g

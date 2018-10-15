#!usr/bin/env python

"""
Copy a cube file, remove all but specified atoms from the atom list

cube_remove_atoms.py filename.cub outfilename.cub 1 2 3 4
cube_remove_atoms.py filename.cub outfilename.cub {1,2,3,4}

specify atoms to NOT be removed

"""

import sys
import numpy as np


f_in = open(sys.argv[1],"r")
f_out = open(sys.argv[2],"w")

if sys.argv[3][0] == '{':
    keep_atoms = [int(x) for x in sys.argv[3][1:-2].split(',')]
    print(keep_atoms)
else:
    keep_atoms = [int(x) for x in sys.argv[3:]]
#print(keep_atoms)

f_out.write(f_in.readline())  # copy first two lines (comments)
f_out.write(f_in.readline())


# Next line: replace number of atoms, keep xyz origins
l = f_in.readline().split()
n_atoms_in = -int(l[0])
n_atoms_out = len(keep_atoms)
f_out.write("{}\t{}\t{}\t{}\t{}\n".format(-n_atoms_out, l[1], l[2], l[3], l[4]))

# Write xyz increment lines
f_out.write(f_in.readline())
f_out.write(f_in.readline())
f_out.write(f_in.readline())


for i in range(1, n_atoms_in+1):
    #print("{}\t".format(i))
    l = f_in.readline()
    if i in keep_atoms:
        #print(i)
        #print(l)
        f_out.write(l)

f_out.write(f_in.readline())

for line in f_in:
    f_out.write(line)

f_in.close()
f_out.close()

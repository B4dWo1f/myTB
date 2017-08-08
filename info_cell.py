#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import numpy as np
import sys


try: fname = sys.argv[1]
except IndexError: 
   print('XYZ file not specified')
   sys.exit(1)


pos = np.loadtxt(fname, skiprows=2,usecols=(1,2,3))

## Center
C = np.mean(pos,axis=0)
pos -= C

X = pos[:,0]
Y = pos[:,1]
Z = pos[:,2]


## Requires hexagon with 1 side parallel to X or Y
lx = np.max(X) - np.min(X)   # diameter
ly = np.max(Y) - np.min(Y)
hex_side = min((lx/2.,ly/2.))

l = hex_side*np.sqrt(3)/2. # RADIO circulo inscrito


print('Information of the cell described in file:',fname)
print('Limits of the cell:')
print('  X:',np.min(X),np.max(X))
print('  Y:',np.min(Y),np.max(Y))
print('  Z:',np.min(Z),np.max(Z))
print('')
print('Hexagon parameters:',' a:',hex_side, ';  a\':',l)

print('Distances between vacancies:')
print('max allowed:',2*l)  # dist be2een parallel sides
print('max recommended:',2*l/3.)

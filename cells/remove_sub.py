#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import os
import sys
try: fnames = sys.argv[1:]
except IndexError:
   print('File not specified')
   exit()
import numpy as np


for fname in fnames:
   print(fname)
   tot = open(fname).read().splitlines()
   Natoms = tot[0]
   vec = tot[1]
   del tot
   ats = np.loadtxt(fname,skiprows=2,usecols=(0,),dtype=bytes)
   ats = [str(x,'utf-8') for x in ats]
   X,Y,Z = np.loadtxt(fname,skiprows=2,usecols=(1,2,3),unpack=True)
   print(X.shape,Y.shape,Z.shape)
   print(ats.shape)
   with open('/tmp/TEST.xyz','w') as f:
      s = '   '
      f.write(Natoms+'\n')
      f.write(str(vec)+'\n')
      for i in range(ats.shape[0]):
         f.write(str(ats[i])+s+str(X[i])+s+str(Y[i])+s+str(Z[i])+'\n')
   f.close()
   os.system('mv /tmp/TEST.xyz %s'%(fname))

#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import numpy as np
import os
import matplotlib.pyplot as plt

com = 'find ../../Documents/data/artificial_lattices/graphene/OUTS/1orb/simple/n45_l2/nv3_na0/d126/alpha0.0/e*/dfct.bands'
files = os.popen(com).read().splitlines()
files = sorted(files,key= lambda x: float(x.split('/')[-2].replace('e','')))

for f in files:
   print(f)
   elec = float(f.split('/')[-2].replace('e',''))
   X0,Y0,Z0 = np.loadtxt(f.replace('dfct','pris'),unpack=True)
   X,Y,Z = np.loadtxt(f,unpack=True)
   fig, ax = plt.subplots()
   ax.scatter(X0,Y0,alpha=0.2)
   ax.scatter(X,Y)
   ax.set_title('$E=%s$ $(eV)$'%(elec),fontsize=20)
   fname = '/tmp/kagome/e%s.png'%(elec)
   fig.savefig(fname)
   with open('/tmp/kagome/files.txt','a') as fi:
      fi.write(fname+'\n')
   fi.close()
   print('  ',fname)
   plt.close(fig)

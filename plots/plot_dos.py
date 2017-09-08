#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import sys

try: fnames = sys.argv[1:]
except IndexError: exit() #fname = ['mag.dat']

Xdos,Ydos = [],[]
for fname in fnames:
   dos = fname
   M = np.loadtxt(dos,unpack=True)
   Xdos.append( M[0] )
   Ydos.append( np.sum(M[1:],axis=0) )


import matplotlib.pyplot as plt
fig, ax = plt.subplots()
for x,y in zip(Xdos,Ydos):
   ax.plot(x,y)
#ax.axhline(0.01)
ax.fill_between(x, 0, 0.01, facecolor='grey', alpha=0.5)
ax.set_xlim([-3.5,3.5])
ax.set_ylim([0,0.61])

plt.show()

#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import sys
import numpy as np
import matplotlib.pyplot as plt


try: fname = sys.argv[1]
except IndexError:
   print('File not specified')
   exit()


fig, ax = plt.subplots()

try: 
   X0,Y0,Z0 = np.loadtxt(fname + 'pris.bands',unpack=True)
   X,Y,Z    = np.loadtxt(fname + 'dfct.bands',unpack=True)
except OSError:
   M0 = np.load(fname + 'pris.bands.npy')
   M = np.load(fname + 'dfct.bands.npy')
   X0 = M0[:,0]
   Y0 = M0[:,1]
   Z0 = M0[:,2]
   X = M[:,0]
   Y = M[:,1]
   Z = M[:,2]

ax.scatter(X0,Y0,alpha=0.2,zorder=10)
ax.scatter(X,Y)

plt.show()

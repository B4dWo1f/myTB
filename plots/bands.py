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

X0,Y0,Z0 = np.loadtxt(fname + 'pris.bands',unpack=True)
X,Y,Z    = np.loadtxt(fname + 'dfct.bands',unpack=True)

ax.scatter(X0,Y0,alpha=0.2)
ax.scatter(X,Y)

plt.show()

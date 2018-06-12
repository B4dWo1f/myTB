#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata

cmap = plt.cm.jet
cmap = plt.cm.viridis

#X,Y,E,A = np.loadtxt('mono_2d.dat',unpack=True)
#X,Y,E,A = np.loadtxt('bi_2d.dat',unpack=True)
#X,Y,E,A = np.loadtxt('bi_2d_e-0.2.dat.1',unpack=True)
X,Y,E,A = np.loadtxt('datos.dat',unpack=True)


def grid_data(X,Y,Z):
   # define grid.
   xi = np.linspace(min(X),max(Y),100)
   yi = np.linspace(min(Y),max(Y),100)
   # grid the data.
   try: zi = griddata(X,Y,Z,xi,yi)
   except:
      print('Problem with natgrid')
      zi = griddata(X,Y,Z,xi,yi,interp='linear')
   return xi,yi,zi

xe,ye,e = grid_data(X,Y,E)
xa,ya,a = grid_data(X,Y,A)


# contour the gridded data, plotting dots at the randomly spaced data points.
fig, ax = plt.subplots()
#CS = ax.contour(xi,yi,zi,15,linewidths=0.5,colors='k')
CS = ax.contourf(xe,ye,e,100,cmap=cmap,vmax=0.0)
fig.colorbar(CS) # draw colorbar
ax.set_xlim([-20,0])
ax.set_ylim([0,20])
ax.set_title('$E_0$ $(eV)$')
fig.tight_layout()

fig, ax = plt.subplots()
#CS = ax.contour(xi,yi,zi,15,linewidths=0.5,colors='k')
CS = ax.contourf(xa,ya,a,100,cmap=cmap)
fig.colorbar(CS) # draw colorbar
#plt.colorbar() # draw colorbar
ax.set_xlim([-20,0])
ax.set_ylim([0,20])
ax.set_title('$\mathcal{A}$ $(MHz)$')
fig.tight_layout()

plt.show()

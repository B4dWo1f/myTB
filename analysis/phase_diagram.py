#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import mwe_exchange as ex

fname = 'datos.dat'

d,e,JF,D,tRL,tLR,UR,UL,E1,E2 = np.loadtxt(fname,unpack=True)

X,Y,Z = [],[],[]
for i in range(len(d)):
   X.append( e[i] )
   Y.append( d[i] )
   H = ex.blue(JF[i],D[i],tRL[i],tLR[i],UR[i],UL[i],E1[i],E2[i])
   es = np.linalg.eigvalsh(H)
   Ea = np.std(es[0:3])
   Eb = np.std(es[1:4])
   #Z.append(JF[i])
   Z.append(Ea-Eb)

X = np.array(X)
Y = np.array(Y)
Z = np.array(Z)

def re_size(x):
   inds_p = np.where(x>=0)[0]
   inds_m = np.where(x<0)[0]

   xp = x[inds_p]
   xm = x[inds_m]

   xp /= np.max(xp)
   xm /= np.max(np.abs(xm))

   ip,im = 0,0
   X = []
   for i in range(len(x)):
      if i in inds_p:
         X.append(xp[ip])
         ip += 1
      else:
         X.append(xm[im])
         im += 1
   return np.array(X)

#Z = re_size(Z)
#Z = np.where(Z>=0,1,-1)



l = np.max(np.abs(Z))


xi = np.linspace(np.min(X),np.max(X),100)
yi = np.linspace(np.min(Y),np.max(Y),100)
# grid the data.
from matplotlib.mlab import griddata
zi = griddata(X,Y,Z,xi,yi,interp='linear')
# contour the gridded data, plotting dots at the randomly spaced data points.


fig, ax = plt.subplots()
cmap = 'tab20b'
cmap = 'coolwarm'
ax.contourf(xi,yi,zi,15,cmap=cmap,vmin=np.min(Z),vmax=np.max(Z),alpha=0.5,zorder=0) #,colors='k')
C = ax.scatter(X,Y,c=Z,cmap=cmap,vmin=np.min(Z),vmax=np.max(Z),zorder=10)
fig.colorbar(C) #, ticks=[-1, 0, 1])
ax.set_xlabel('$\mathcal{E}$ $(eV)$')
ax.set_ylabel('$d$ $(\AA)$')
import datetime as dt
now = dt.datetime.now()
name = now.strftime('%Y%m%d_%H%M%S.png')
#fig.savefig(name)

zflat = np.sort(Z.flatten())

fig, ax = plt.subplots()
ax.plot(zflat)
plt.show()

#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import os
import numpy as np
from exchange import Spectrum
import matplotlib.pyplot as plt

R = '../../Documents/data/Vpps/'
#R = '../../Documents/data/Vsss/'
tail = 'OUTs/4orb/ac/n35_l1/nv0_na1/d0.0/alpha0.0/e0.0/'
fols = os.popen('ls %s'%(R)).read().splitlines()
fols = sorted(map(float,fols))
X,Y,val,con = [],[],[],[]
A = []
for fo in fols:
   f = R + str(fo) + '/' + tail
   try: S = Spectrum(f)
   except FileNotFoundError: continue
   X.append(fo)
   Y.append(S.E[4])
   val.append(S.E[3])
   con.append(S.E[5])
   v = S.V[4]
   v = v*np.conj(v)
   hyper = v[-1]*1420
   A.append(hyper)
   print(fo,'',S.E[4],'',hyper)

fig, ax = plt.subplots()
ax.axhline(0.0,c='k',ls='--')
ax.plot(X,val,'b--')
ax.plot(X,con,'b--')
ax.plot(X,Y,'r-')
ax.scatter(X,Y,c='r',edgecolors='none')
ax.set_xlim([-100,100])
ax.set_ylim([-0.06,0.06])
ax.set_ylabel('In-gap state')
ax.set_xlabel(R.split('/')[-2])
ax.grid()


import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.plot(X,hyper)
plt.show()

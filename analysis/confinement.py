#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import os
import numpy as np
from exchange import Spectrum

R = '../../Documents/data/confinement/OUTs/4orb/ac/'
tail = '/nv0_na1/d0.0/alpha0.0/e0.0/'

fols = os.popen('ls %s'%(R)).read().splitlines()
N,L,G,A = [],[],[],[]
for f in fols:
   n,l = f.split('_')
   n = int(n.replace('n',''))
   l = int(l.replace('l',''))
   fname = R+f+tail
   try: S = Spectrum(fname)
   except: continue
   print(f)
   N.append(n)
   L.append(l)
   G.append(S.gapP)
   v = S.V[-1]
   v = v*np.conj(v)
   hyper = v[-1]*1024
   A.append(hyper)

N = np.array(N)
L = np.array(L)
G = np.array(G)
A = np.array(A)

N1 = N[L==1]
G1 = G[L==1]
A1 = A[L==1]
inds = np.argsort(N1)
N1 = N1[inds]
G1 = G1[inds]
A1 = A1[inds]

N2 = N[L==2]
G2 = G[L==2]
A2 = A[L==2]
inds = np.argsort(N2)
N2 = N2[inds]
G2 = G2[inds]
A2 = A2[inds]

print('l=1')
for n,g,a in zip(N1,G1,A1):
   print(n,'',g,'',a)
print('\nl=2')
for n,g,a in zip(N2,G2,A2):
   print(n,'',g,'',a)

import matplotlib.pyplot as plt
fig, ax = plt.subplots()
#print(N1)
print(N2)
ax.plot(N1,G1,'o-',label='monolayer')
ax.plot(N2,G2,'o-',label='bilayer')
ax.legend()
ax.grid()
plt.show()

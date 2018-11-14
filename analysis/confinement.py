#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import os
import numpy as np
from mwe_exchange import Spectrum
#from exchange import Spectrum

R = '../../../Documents/data/artificial_lattices/confinement/OUTS/4orb/ac/'
tail = '/nv0_na1/dNone/alpha0.0/e0.0/'
tail = '/nv0_na1/dNone/alpha0.0/e-0.2/'

fols = os.popen('ls %s'%(R)).read().splitlines()
print('%s fols to analyze'%(len(fols)))
N,L,G,A = [],[],[],[]
for f in fols:
   n,l = f.split('_')
   n = int(n.replace('n',''))
   l = int(l.replace('l',''))
   fname = R+f+tail
   try: S = Spectrum(fname)
   except:
      print('Error reading',fname)
      continue
   N.append(n)
   L.append(l)
   G.append(S.gapP)
   v = S.V[-1]
   v = v*np.conj(v)
   v = np.bincount(S.Ba.n, weights=v)
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
ax.set_ylabel('$\delta$ $(eV)$')
ax.set_xlabel('N')

fig, ax = plt.subplots()
ax.plot(N1,A1,'o-',label='monolayer')
ax.plot(N2,A2,'o-',label='bilayer')
ax.legend()
ax.set_ylabel('$\mathcal{A}$ $(MHz)$')
ax.set_xlabel('N')
plt.show()

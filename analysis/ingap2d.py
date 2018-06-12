#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import os
import numpy as np
from exchange import Spectrum

#R = '../../Documents/data/space_parameters/bi/'
#tail = 'OUTs/4orb/ac/n25_l2/nv0_na1/d0.0/alpha0.0/e0.0/'
R = '../../../Documents/data/space_parameters/bi/'
R = '../../../Documents/data/space_parameters/e0.2/'
tail = 'OUTs/4orb/ac/n25_l2/nv0_na1/d0.0/alpha0.0/e0.2/'

fout = 'datos.dat'
## Reset file
f = open(fout,'w')
f.close()

fols = os.popen('ls %s'%(R)).read().splitlines()
Vsss,Vsps,E0,A = [],[],[],[]
for f in fols:
   s,p = map(float,f.split('_'))
   if True: #p>0.0: # and s>0.0:
      #print(R + f)
      try: S = Spectrum(R + f + '/' + tail   )
      except FileNotFoundError: continue
      except: continue
      Vsss.append(s)
      Vsps.append(p)
      e = S.E[4]
      E0.append( e) #-S.cond )
      v = S.V[4]
      v = v*np.conj(v)
      hyper = v[-1] * 1420
      A.append(hyper)
      print(s,'',p,'',e,'',hyper)
      with open(fout,'a') as f:
         f.write(str(s)+'  '+str(p)+'  '+str(e)+'  '+str(hyper)+'\n')

import matplotlib.pyplot as plt
fig, ax = plt.subplots()
C = ax.scatter(Vsss,Vsps,c=E0,s=100,edgecolors='none')
ax.axhline(7.81,c='k',ls='--')
ax.axvline(-6.84,c='k',ls='--')
ax.set_xlabel('$Vss\sigma$',fontsize=15)
ax.set_ylabel('$Vsp\sigma$',fontsize=15)
L=20
ax.set_xlim([-L,0])
ax.set_ylim([0,L])
fig.colorbar(C)
#plt.colorbar(C)

fig, ax = plt.subplots()
H = ax.scatter(Vsss,Vsps,c=A,s=100,edgecolors='none')
ax.axhline(7.81,c='k',ls='--')
ax.axvline(-6.84,c='k',ls='--')
ax.set_xlabel('$Vss\sigma$',fontsize=15)
ax.set_ylabel('$Vsp\sigma$',fontsize=15)
ax.set_xlim([-L,0])
ax.set_ylim([0,L])
fig.colorbar(H)
plt.show()

exit()
fols = sorted(map(float,fols))
X,Y,val,con = [],[],[],[]
for fo in fols:
   f = R + str(fo) + '/' + tail
   try: S = Spectrum(f)
   except FileNotFoundError: continue
   print(fo,'',S.E[4])
   X.append(fo)
   Y.append(S.E[4])
   val.append(S.E[3])
   con.append(S.E[5])

import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.axhline(0.0,c='k',ls='--')
ax.plot(X,val,'b--')
ax.plot(X,con,'b--')
ax.plot(X,Y,'r-')
ax.scatter(X,Y,c='r',edgecolors='none')
ax.set_ylabel('In-gap state')
ax.set_xlabel(R.split('/')[-2])
ax.grid()
plt.show()

#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import sys

try: fname = sys.argv[1]
except IndexError: exit() #fname = ['mag.dat']

try: tag = sys.argv[2]
except IndexError: tag = 'dfct'

#fname = '/tmp/OUTs/simple1_l2/nv0_d0.0_alpha0.0/e0.0/'
ban = fname + tag + '.bands'
dos = fname + tag + '.dos'
bas = fname + tag + '.basis'
ind,nat = np.loadtxt(bas,usecols=(0,1),dtype=int,unpack=True)
ats,orbs,sub = np.loadtxt(bas,usecols=(2,3,4),dtype=bytes,unpack=True)
ats = np.array([str(x,'utf-8') for x in ats])
orbs = np.array([str(x,'utf-8') for x in orbs])
sub = np.array([str(x,'utf-8') for x in sub])


#Xdos,Ydos = np.loadtxt(dos,unpack=True)
print(dos)
M = np.loadtxt(dos,unpack=True)
Mm = M[1:][orbs=='pz']
Xdos = M[0]
Ydos = np.sum(M[1:],axis=0)
Yesp = np.sum(Mm,axis=0)
Xban,Yban,Zban = np.loadtxt(ban,unpack=True)

print('-- Bands ----')
for x,y,z in zip(Xban,Yban,Zban):
   print(x,y,z)
print('\n-- DOS ----')
for x,y,ys in zip(Xdos,Ydos,Yesp):
   print(x,y,ys)

from matplotlib import gridspec

fig = plt.figure()
gs = gridspec.GridSpec(1, 2)
fig.subplots_adjust(wspace=0.,hspace=0.1)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[0,1],sharey=ax1)
## Bands
ax1.scatter(Xban,Yban,c='b',edgecolor='none')
#ax1.scatter(Xban,Yban,c=Zban,edgecolor='none')
ax1.set_xlim([min(Xban),max(Xban)])
ax1.set_ylim([-10,10])
ax1.set_xticks([])
ax1.set_xlabel('K-Path')
ax1.set_ylabel('$E$ $(eV)$',fontsize=17)
## DOS
ax2.plot(Ydos,Xdos,'b-',lw=2,label='$DOS$')
ax2.plot(Yesp,Xdos,'r-',lw=2,label='$DOS_s$')
ax2.axvspan(0,0.01, alpha=0.5, color='k')
ax2.set_xlim([0,1.75])
ax2.set_xlabel('$DOS$ $(a.u.)$',fontsize=17)

ax2.grid()
ax2.legend(loc=4)

plt.show()

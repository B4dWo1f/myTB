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
bas = fname + tag + '.basis'
ind,nat = np.loadtxt(bas,usecols=(0,1),dtype=int,unpack=True)
ats,orbs,sub = np.loadtxt(bas,usecols=(2,3,4),dtype=bytes,unpack=True)
ats = np.array([str(x,'utf-8') for x in ats])
orbs = np.array([str(x,'utf-8') for x in orbs])
sub = np.array([str(x,'utf-8') for x in sub])


#Xdos,Ydos = np.loadtxt(dos,unpack=True)
try: Xban,Yban,Zban = np.loadtxt(ban,unpack=True)
except:
   M = np.load(ban+'.npy')
   Xban = M[:,0]
   Yban = M[:,1]
   Zban = M[:,2]

print('-- Bands ----')
for x,y,z in zip(Xban,Yban,Zban):
   print(x,y,z)


from matplotlib.colors import LinearSegmentedColormap
cdict={'red'  : ((0.,   0,   0),(0.6,0.0,0.0),(1, 1.0, 1.0)), 'green': ((0., 0.0, 0.0),(0.4,1.0,1.0),(0.6,1.0,1.0),(1, 0.0, 0.0)), 'blue' : ((0., 1.0, 1.0),(0.4,0.0,0.0),(1, 0.0, 0.0))}
my_cmap = LinearSegmentedColormap('my_colormap',cdict,256)



import matplotlib.pyplot as plt
fig, ax = plt.subplots()
## Bands
ax.scatter(Xban,Yban,c=Zban)#,cmap=my_cmap,edgecolor='none')
#ax1.scatter(Xban,Yban,c=Zban,edgecolor='none')
ax.set_xlim([min(Xban),max(Xban)])
ax.set_xticks([])
ax.set_xlabel('K-Path')
ax.set_ylabel('$E$ $(eV)$',fontsize=17)
ax.grid()

plt.show()

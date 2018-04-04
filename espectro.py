#!/usr/bin/python3
# -*- coding: UTF-8 -*-


def state_space(v,X,Y,Z,S=100):
   """
     Plot real space distribution of a given wave-function
      v: vector to plot, with shape: (N,) 
      X,Y,Z: x,y,z coordinates of each of the atoms
   """
   #y = [(v[i]*np.conj(v[i])).real for i in range(v.shape[0])]
   xx = np.linspace(-5,5,50)
   yy = np.sin(xx)
   fig, ax = plt.subplots()
   #ax.scatter(X,Y,c='k',alpha=0.5)
   #ax.scatter(X,Y,c='r',s=S*y,alpha=0.5)
   ax.plot(xx,yy)
   ax.set_aspect('equal')
   plt.show()


def state_base(v,base):
   """
     Plot each of the components of a wave-function in a provided basis
      v: vector to plot, with shape: (N,)
      base: names of basis states, list of length N
   """
   y = [(v[i]*np.conj(v[i])).real for i in range(v.shape[0])]
   fig, ax = plt.subplots() #figsize=(17,3))
   ind = np.array(range(v.shape[0]))
   width = 0.7
   ax.bar(ind,y)
   ax.set_xticks(ind + width/1.9)
   ax.set_xticklabels(base)
   ax.set_ylim([0,1])
   ax.grid()
   fig.tight_layout()
   plt.show()


def picker_wrapper(maxd=0.075):
   def picker(line, Mevent):
      """
      find the points within a certain distance from the mouseclick in
      data coords and attach some extra attributes, pickx and picky
      which are the data points that were picked
      """
      if Mevent.xdata is None:
         return False, dict()
      x0 = line.get_xdata()
      y0 = line.get_ydata()
      d = np.sqrt((x0 - Mevent.xdata)**2. + (y0 - Mevent.ydata)**2.)
      ind = np.nonzero(np.less_equal(d, maxd))
      if len(ind):
         pickx = np.take(x0, ind)
         picky = np.take(y0, ind)
         props = dict(ind=ind, pickx=pickx, picky=picky)
         return True, props
      else:
         return False, dict()
   return picker


def onpick_wrapper(vecs,base,plot):
   def onpick(event):
      try:
         ind = event.pickx[0,0]
         v = vecs[ind,:]
         print('Selected:',ind)
         plot(v,base)
      except IndexError: pass   #print('no data selected')
   return onpick


import numpy as np
import exchange as ex
import matplotlib.pyplot as plt

fname = '../../Documents/data/OUTs/1orb/ac/n20_l2/nv1_na0/d0.0/alpha0.0/e1.0/dfct_spectrum.npy'
M = np.load(fname)
E = M[:,0]
V = M[:,1:]
names = ['' for _ in E]

exit()

my_onpick = onpick_wrapper(V,names,state_space)
my_picker = picker_wrapper()

fig, ax = plt.subplots()
line, = ax.plot(range(len(E)),E, 'o', picker=my_picker)
ax.set_xlim([-1,len(E)+1])
ax.grid()

fig.canvas.mpl_connect('pick_event', my_onpick)
plt.show()

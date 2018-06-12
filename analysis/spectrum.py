#!/usr/bin/python3
# -*- coding: UTF-8 -*-


"""
 this will plot the  spectrum for different electric fields
"""

import numpy as np
import exchange as ex
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import sys
here = os.path.dirname(os.path.realpath(__file__)) + '/'
mpl.rcParams['font.size'] = 15

## Plot hyperfine & gap ########################################################
def plot_exchange(e):
   print('Plotting spectrum for E=',e)
   print('... Please wait a few moments')
   ffol = fol+'e%s/'%(e)
   com = 'python3 exchange.py %s &'%(ffol)
   os.system(com)

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
      #ind = np.nonzero(np.less_equal(d, maxd))
      ind = (np.argmin(d),)
      if len(ind):
         pickx = np.take(x0, ind)
         picky = np.take(y0, ind)
         props = dict(ind=ind, pickx=pickx, picky=picky)
         return True, props
      else:
         return False, dict()
   return picker

def onpick_wrapper(plot):
   def onpick(event):
      try:
         elec = event.pickx[0]
         plot(elec)
      except IndexError: pass   #print('no data selected')
   return onpick

################################################################################
try: fol = sys.argv[1]
except IndexError: 
   print('No folder provided')
   exit()

folders = []
for a in os.walk(fol):
   folders.append( a[0]+'/' )
folders = folders[1:]
folders = sorted(folders,key=lambda x: float(x.split('/')[-2][1:]))

print('Analyzing %s folders'%(len(folders)))

X,hyper,gap,gapP = [],[],[],[]
Xplt,Yplt,YPplt = [],[],[]
SP,SL = [], []
#gapP,gapD,lc,E0 = [],[],[],[]
for f in folders:
   #A = ex.Spectrum(f,nv=1)
   try: A = ex.Spectrum(f,nv=1)
   except:
      print('ERROR reading:',f)
      continue
   v = A.V_ingap[0,:]
   vv = np.conj(v) * v
   SL.append(A.SL_eig[0])
   SP.append(A.SP_eig[0])
   X.append(A.elec)
   hyper.append( vv[-1]*1420)
   print(A.elec, vv[-1]*1420)
   gap.append(A.gap)
   gapP.append(A.gapP)
   for e,ep in zip(A.E,A.Ep):
      Xplt.append(A.elec)
      Yplt.append(e)
      YPplt.append(ep)

## Plot spectrums ##############################################################

my_onpick = onpick_wrapper(plot_exchange)
my_picker = picker_wrapper()

#fig, (ax,ax1,ax2) = plt.subplots(3,1,sharex=True)
import matplotlib.gridspec as gridspec
fig = plt.figure()
gs = gridspec.GridSpec(4, 1)
ax = plt.subplot(gs[0:2, 0])
ax1=plt.subplot(gs[2, 0])
ax2=plt.subplot(gs[3, 0])

ax.scatter(Xplt,Yplt,c='r',s=50,edgecolors='none')
ax.scatter(Xplt,YPplt,c='b',s=20,edgecolors='none',alpha=0.7)
#ax.set_xlabel('$\lambda_E$ $(eV)$',fontsize=15)
ax.set_ylabel('$E$ $(eV)$',fontsize=15)

line, = ax1.plot(X,hyper,'o-',picker=my_picker)
ax1.set_ylim(ymin=0) #[-1,65])
#ax1.set_xlabel('$\lambda_E$ $(eV)$',fontsize=15)
ax1.set_ylabel('$\mathcal{A}$ $(MHz)$',fontsize=15)

ax2.plot(X,gap,'ro-',picker=my_picker)
ax2.plot(X,gapP,'k--')
ax2.set_ylim(ymin=0)
ax2.set_xlabel('$\lambda_E$ $(eV)$',fontsize=15)
ax2.set_ylabel('$\Delta$ $(eV)$',fontsize=15)

ax.set_xlim([-0.2,0.2])
ax1.set_xlim([-0.2,0.2])
ax2.set_xlim([-0.2,0.2])
ax.grid()
ax1.grid()
ax2.grid()

#fig.subplots_adjust(hspace=0)
fig.canvas.mpl_connect('pick_event', my_onpick)
fig.tight_layout()


import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.plot(X,SL,label='layer')
ax.plot(X,SP,label='sublatt')
ax.set_ylim([-1,1])
ax.grid()
ax.legend()
plt.show()
plt.show()

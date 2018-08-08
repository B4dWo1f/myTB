#!/usr/bin/python3
# -*- coding: UTF-8 -*-

"""
  This script shows the evolution with the electric field of the spectrum and
  some other properties.
  Default: only spectrum (faster)
  Other quantities:
     - [P] layer/sublattice polarization
     - [G] gap
     - [E0] in-gap energy
     X- [A] Hyperfine   # not yet
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from tqdm import tqdm
import mwe_exchange as ex

## Plot hyperfine & gap ########################################################
def plot_exchange(e):
   print('Plotting spectrum for E=',e)
   print('... Please wait a few moments')
   ffol = fol+'e%s/'%(e)
   com = 'python3 mwe_exchange.py %s &'%(ffol)
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

my_onpick = onpick_wrapper(plot_exchange)
my_picker = picker_wrapper()
################################################################################


import sys
try: fols = sys.argv[1:]
except IndexError:
   print('Folder not specified')
   exit()



for fol in fols:
   folders = []
   for a in os.walk(fol):
      folders.append( a[0]+'/' )
   folders = folders[1:]
   folders = sorted(folders,key=lambda x: float(x.split('/')[-2][1:]))
   folders = folders


   X,P,L,G,E0 = [],[],[],[],[]
   Xplt,Yplt,YPplt = [],[],[]
   for f in tqdm(folders):
      A = ex.Spectrum(f)
      #try: A = ex.Spectrum(f)
      #except:
      #   print('ERROR reading:',f)
      #   continue
      for e,ep in zip(A.E,A.Ep):
         Xplt.append(A.elec)
         Yplt.append(e)
         YPplt.append(ep)
      X.append(A.elec)
      P.append(A.SP)
      L.append(A.SL)
      G.append(A.gap)
      E0.append(A.E_ingap)
   mx,Mx = np.min(X), np.max(X)

   fig = plt.figure(figsize=(9,10))
   gs = gridspec.GridSpec(5, 1)
   #fig, ax = plt.subplots()
   ax = plt.subplot(gs[0:2, 0])
   ax_P  = plt.subplot(gs[2, 0], sharex=ax)
   ax_G  = plt.subplot(gs[3, 0], sharex=ax)
   ax_E0 = plt.subplot(gs[4, 0], sharex=ax)

   # Spectrum
   ax.plot(Xplt,YPplt,'o',markersize=10,alpha=0.4)
   line, = ax.plot(Xplt,Yplt,'o',picker=my_picker)
   ax.set_ylabel('$E$ $(eV)$')
   ax.set_xlim([mx,Mx])

   # Polarization
   P = np.array(P)
   L = np.array(L)
   for i in range(P.shape[1]):
      iP = P[:,i]
      iL = L[:,i]
      l, = ax_P.plot(X,iP,lw=2)
      ax_P.plot(X,iL,color=l.get_color(),ls='--',lw=2)
   ax_P.set_ylabel('Polarization')
   ax_P.set_ylim([-1,1])
   ax_P.set_xlim([mx,Mx])

   # Gap
   ax_G.plot(X,G,lw=2)
   ax_G.set_ylabel('$\Delta$ $(eV)$')
   ax_G.set_xlim([mx,Mx])

   # In-gap Energies
   ax_E0.plot(X,E0,lw=2)
   ax_E0.set_ylabel('$E_0$ $(eV)$')
   ax_E0.set_xlim([mx,Mx])

   fig.canvas.mpl_connect('pick_event', my_onpick)
   fig.tight_layout()

plt.show()

#!/usr/bin/python3
# -*- coding: UTF-8 -*-


"""
 this will plot the  spectrum for different electric fields
"""


import numpy as np
import exchange as ex
import os
import sys
here = os.path.dirname(os.path.realpath(__file__)) + '/'


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

#from scipy.constants import physical_constants
#import pint
#
#U = pint.UnitRegistry()
#
#
### Electron Charge
#e = physical_constants['atomic unit of charge']
#e = e[0] * U.parse_expression(e[1])
#print('e =',e)
#
### Atomic stuff
#rb = physical_constants['Bohr radius']
#rb = rb[0] * U.parse_expression(rb[1])
#z0 = 3*rb*0.620*U.angstrom/(0.529*U.angstrom)
#print('z0 =',z0)
#
### Slater-Koster parameter
#Vsps = 5.58*U.eV
#print('Vsps =',Vsps)
#
### atomic SOC
#xi = 6 * U.meV
#print('xi =',xi)

#xxx = (0.01*U.eV*3*Vsps)/(e*z0*xi)
#print(xxx.to('V/nm'))
#exit()

X,hyper,gap,gapP = [],[],[],[]
Xplt,Yplt,YPplt = [],[],[]
#gapP,gapD,lc,E0 = [],[],[],[]
for f in folders:
   #A = ex.Spectrum(f,nv=1)
   try: A = ex.Spectrum(f,nv=1)
   except:
      print('ERROR reading:',f)
      continue
   v = A.V_ingap[0,:]
   vv = np.conj(v) * v
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
import matplotlib.pyplot as plt
fig1, ax1 = plt.subplots()
ax1.scatter(Xplt,Yplt,c='r',edgecolors='none')
ax1.scatter(Xplt,YPplt,c='b',edgecolors='none',alpha=0.7)
ax1.set_xlabel('$\lambda_E$ $(eV)$',fontsize=15)
ax1.set_ylabel('$E$ $(eV)$',fontsize=15)


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


my_onpick = onpick_wrapper(plot_exchange)
my_picker = picker_wrapper()

fig, (ax,ax1) = plt.subplots(2,1,sharex=True)
line, = ax.plot(X,hyper,'o-',picker=my_picker)
ax.set_ylim([-1,50])
ax.set_xlabel('$\lambda_E$ $(eV)$',fontsize=15)
ax.set_ylabel('$\mathcal{A}$ $(MHz)$',fontsize=15)
ax.grid()

#fig1, ax1 = plt.subplots()
ax1.plot(X,gap,'ro-',picker=my_picker)
ax1.plot(X,gapP,'k--')
ax1.set_ylim(ymin=0)
ax1.set_xlabel('$\lambda_E$ $(eV)$',fontsize=15)
ax1.set_ylabel('$\Delta$ $(eV)$',fontsize=15)
ax1.grid()

#fig.subplots_adjust(hspace=0)
fig.canvas.mpl_connect('pick_event', my_onpick)

plt.show()

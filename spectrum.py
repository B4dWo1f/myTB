#!/usr/bin/python3
# -*- coding: UTF-8 -*-


"""
 this will plot the  spectrum for different electric fields
"""


import numpy as np
import exchange as ex
import os
import sys


try: fol = sys.argv[1]
except IndexError: 
   print('No folder provided')
   exit()

folders = []
for a in os.walk(fol):
   folders.append( a[0]+'/' )
folders = folders[1:]
folders = sorted(folders,key=lambda x: float(x.split('/')[-2].replace('e','')))

print('Analyzing %s folders'%(len(folders)))

X,hyper = [],[]
Xplt,Yplt,YPplt = [],[],[]
for f in folders:
   try: A = ex.Spectrum(f)
   except: continue
   X.append( A.elec )
   Xplt.append( [A.elec for _ in A.E] )
   Yplt.append( A.E )
   YPplt.append( A.Ep )
   v = A.V_ingap[0]
   vv = np.conj(v) * v
   hyper.append(vv[-1]*1420)

import matplotlib.pyplot as plt
fig, ax = plt.subplots()
mx,Mx = 1000,-1000
for x,y,yp in zip(Xplt,Yplt,YPplt):
   ax.scatter(x,y,c='b',s=50,edgecolor='none',alpha=0.5,zorder=0)
   ax.scatter(x,yp,c='k',s=20,edgecolor='none',alpha=0.5,zorder=1)
   mx = min([min(x),mx])
   Mx = max([max(x),Mx])

ax.set_xlim([mx,Mx])
ax.grid()


import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.plot(X,hyper)
ax.grid()


plt.show()



#import matplotlib.pyplot as plt
#import numpy as np
#import os
#import sys


#def find_nearest(array,value): return (np.abs(array-value)).argmin()
#
### Decide sign
#def get_sign(v,e):
#   ind_p = find_nearest(v, e)
#   ind_m = find_nearest(v,-e)
#   if abs(v[ind_p]-e) < abs(v[ind_m]+e): return e
#   else: return -e
#
#def mean_val(v,M): return np.sum(np.conj(v)*M*v)
#
#
#
#try: fol = sys.argv[1]
#except IndexError: fol = 'OUTs/ac30_l2/nv2_d120.0_alpha0.0/'
#
#Ueff = 1.5 * 2.7
#
#
#
#folders = []
#for a in os.walk(fol):
#   folders.append( a[0]+'/' )
#folders = folders[1:]
#folders = sorted(folders,key=lambda x: float(x.split('/')[-2].replace('e','')))
#
#print('Analyzing %s folders'%(len(folders)))
#
#
### Parse information from folder structure
#l = fol.split('/')
#Nav = l[-4]
#dist = l[-3]
#alpha = l[-2]
#
#Nv = int(Nav.split('_')[0].replace('nv',''))
#Na = int(Nav.split('_')[1].replace('na',''))
#Nv = Na   # XXX
#dist = float(dist.replace('d',''))
#alpha = float(alpha.replace('alpha',''))
#
##patt = r'nv(\S+)_d(\S+)_alpha(\S+)'
##import re
##match = re.search(patt,l,re.UNICODE)
##Nv,dist,alpha = match.groups()
##try:
##   Nv = int(Nv)
##   dist = float(dist)
##   alpha = float(alpha)
##except: 
##   print('ERROR reading Nv, dist, or alpha')
##   exit()
#
#
#
#elecs = []
## Variables to store sublattice and layer polarization
#Xplt,Ysubs,Ylays = [],[[] for _ in range(Nv)],[[] for _ in range(Nv)]
## Variables to store spectrum
#Xs,Ys,Ysp = [],[],[]
## Variables to store exchange couplings
#J_f,J_af = [],[]
#for fol in folders:
#   ## Files
#   pris = fol + 'pris_spectrum.npy'
#   dfct = fol + 'dfct_spectrum.npy'
#   pris_pos = fol + 'base_pris.xyz'
#   dfct_pos = fol + 'base_dfct.xyz'
#   log = fol + 'main.log'
#
#   ## Get Vacancies indices
#   vacs = os.popen('grep " Changing onsite of atom:" %s'%(log)).readlines()
#   vacs = [v.split(':')[-1] for v in vacs]
#   vacs = [int(v.lstrip().rstrip()) for v in vacs]
#
#   ## Get Electric Field
#   elec = float(fol.split('/')[-2].replace('e',''))
#
#   print('Plotting states from:',dfct)
#   ## Read the Eigen-Energies and Eigen-Functions
#   # Pristine
#   M = np.load(pris)
#   Ep = M[:,0]
#   if np.max(np.abs(Ep.imag)) > 0.00000001: print('COMPLEX ENERGIES')
#   Ep = Ep.real
#   # Defected
#   M = np.load(dfct)
#   E = M[:,0]
#   if np.max(np.abs(E.imag)) > 0.00000001: print('COMPLEX ENERGIES')
#   E = E.real
#   V = M[:,1:]
#
#   pos = np.loadtxt(dfct_pos,skiprows=2,usecols=(1,2,3))
#   Rx = pos[:,0]   # Operator position
#   Z = pos[:,2]
#
#
#   # Store data for spectrum
#   x,xp = [],[]
#   for _ in E:
#      x.append(elec)
#      xp.append(elec)
#   elecs.append(elec)
#   Xs.append( np.array(x) )
#   Ys.append( E )
#   Ysp.append( Ep )
#
#
#   # Treat data for Polarizations
#   ## Sublattice operator
#   sub = np.loadtxt(dfct_pos,skiprows=2,usecols=(4,),dtype=str)
#   sub = np.array([a.replace('b\'','') for a in sub])
#   sub = np.array([a.replace('\'','') for a in sub])
#   subdic = {'A':1,'B':-1}     # sublatice to num
#   sub = np.array([subdic[s] for s in sub])
#   ## Layer operator
#   layer = np.where(Z>0,1,-1)
#
#   ## Get conduction, valence, and in-gap states and indices
#   absE = np.abs(E)
#   ind_E = np.argsort(absE)
#   N = Nv   #XXX this assumes 2 in-gap states
#   ingap = [get_sign(E,e) for e in absE[ind_E[0:N]]]
#   EE = E
#   for ein in ingap:
#      EE = EE[EE!=ein]
#   cond = np.max(EE[EE<0])
#   vale = np.min(EE[EE>0])
#   ## indices
#   inds = [find_nearest(E,ie) for ie in ingap]
#   Xplt.append(elec)
#   for ii in range(Nv):
#      s = mean_val(V[inds[ii]],sub)
#      l = mean_val(V[inds[ii]],layer)
#      print('E:',E[inds[ii]])
#      print('  Sublattice Polarization:',s)
#      print('  Layer Polarization:',l)
#      Ysubs[ii].append(s)
#      Ylays[ii].append(l)
#   ## Define and sort left and right states
#   if Nv != 2: continue
#   v1 = (V[inds[0]] + V[inds[1]])/2   #XXX np.sqrt(2)?
#   v2 = (V[inds[0]] - V[inds[1]])/2   #XXX np.sqrt(2)?
#   Rx = pos[:,0]   # Operator position
#   A = mean_val(v1,Rx)
#   B = mean_val(v2,Rx)
#   if A > B: vs = [v2,v1]      # Left, Right
#   else: vs = [v1,v2]          # Left, Right
#   vL,vR = vs
#   F = 0.
#   for iv in range(len(vL)):
#      F += ( vL[iv]*np.conj(vL[iv]) ) * ( vR[iv]*np.conj(vR[iv]) )
#   J_f.append(F*2*Ueff)
#   J_af.append( (np.abs(ingap[0]-ingap[1])**2)/Ueff )
#
#
#
##### Plots #####################################################################
#import matplotlib.pyplot as plt
#fig, ax = plt.subplots()
##from matplotlib import gridspec
##fig = plt.figure()
##gs = gridspec.GridSpec(3, 1)
##fig.subplots_adjust(wspace=0.,hspace=0.)
##ax1 = plt.subplot(gs[0,0])   # Spectrum
##ax2 = plt.subplot(gs[1,0], sharex=ax1)   # Layer
##ax3 = plt.subplot(gs[2,0], sharex=ax1)   # Sublattice
#
#print('     ------ Spectrum -----')
#d = 0.0005
#for x,y,yp in zip(Xs,Ys,Ysp):
#   for ix,iy,iyp in zip(x,y,yp):
#      print(ix,iy,iyp)
#   I = np.linspace(1,1,len(x))
#   ax.scatter(x-d*I,yp,c='b',s=50,alpha=0.3,zorder=0)   # pris
#   ax.scatter(x+d*I,y,c='k',s=20,alpha=0.5,zorder=1)    # dfct
#ax.set_xlim(np.min(Xs)-3*d,np.max(Xs)+3*d)
#ax.set_ylim(np.min(Ys)-3*d,np.max(Ys)+3*d)
#ax.set_xlabel('$\lambda_E$ $(eV)$',fontsize=17)
#ax.set_ylabel('$E$ $(eV)$',fontsize=17)
#
#
#ax.grid()
#
#
#plt.show()

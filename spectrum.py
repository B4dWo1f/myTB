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
folders = sorted(folders,key=lambda x: float(x.split('/')[-2][1:]))

print('Analyzing %s folders'%(len(folders)))

from scipy.constants import physical_constants
import pint

U = pint.UnitRegistry()


## Electron Charge
e = physical_constants['atomic unit of charge']
e = e[0] * U.parse_expression(e[1])
print('e =',e)

## Atomic stuff
rb = physical_constants['Bohr radius']
rb = rb[0] * U.parse_expression(rb[1])
z0 = 3*rb*0.620*U.angstrom/(0.529*U.angstrom)
print('z0 =',z0)

## Slater-Koster parameter
Vsps = 5.58*U.eV
print('Vsps =',Vsps)

## atomic SOC
xi = 6 * U.meV
print('xi =',xi)


X,hyper = [],[]
Xplt,Yplt,YPplt = [],[],[]
for f in folders:
   try: A = ex.Spectrum(f,nv=2)
   except: continue
   elec = (A.elec*U.eV*3*Vsps)/(e*z0*xi)
   elec = elec.to('V/nm').magnitude
   X.append( elec)  #A.elec )
   Xplt.append( [A.elec for _ in A.E] )
   Yplt.append( A.E )
   YPplt.append( A.Ep )
   #v = A.V_ingap[0]
   #vv = np.conj(v) * v
   #hyper.append(vv[-1]*1420)
   hyper.append( abs(A.E_ingap[0]-A.E_ingap[1]) )
X = np.array(X)
hyper = np.array(hyper)
h0 = hyper[X==np.min(np.abs(X))]

import matplotlib.pyplot as plt
fig, ax = plt.subplots()
mx,Mx = 1000,-1000
for x,y,yp in zip(Xplt,Yplt,YPplt):
   ax.scatter(x,yp,c='b',s=50,edgecolor='none',alpha=0.5,zorder=0)
   ax.scatter(x,y,c='k',s=50,edgecolor='none',alpha=0.5,zorder=1)
   mx = min([min(x),mx])
   Mx = max([max(x),Mx])

#xx = np.linspace(mx,MX)
ax.set_xlim([mx,Mx])
ax.grid()


import matplotlib.pyplot as plt
print('#V/m   meV')
for x,h in zip(X,(hyper-h0)*1000):
   print(x,' ',h)
fig, ax = plt.subplots()
ax.plot(X,(hyper-h0)*1000,'o-')
ax.set_xlim([min(X),max(X)])
ax.set_xlabel('$E$ $(V/nm)$',fontsize=20)
ax.set_ylabel('$\Delta-\Delta_0$ $(meV)$',fontsize=20)
ax.grid()
fig.tight_layout()

plt.show()

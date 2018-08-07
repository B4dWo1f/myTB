#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import mwe_exchange as ex
from tqdm import tqdm
import numpy as np
import os


def get_d_e(f):
   if f[-1] != '/': f += '/'
   ff = f.split('/')
   d = float(ff[-4].replace('d',''))
   e = float(ff[-2].replace('e',''))
   return d,e

fs = ['../../../Documents/data/artificial_lattices/dimer/OUTS/1orb/ac/n50_l2/nv2_na0', '../../../Documents/data/artificial_lattices/dimer_random/OUTS/1orb/ac/n50_l2/nv2_na0']
fs = ['../../../Documents/data/artificial_lattices/dimer/OUTS/1orb/ac/n50_l2/nv2_na0']
fout = 'datos_grid_30.dat'

fols = []
for f in fs:
   for root, dirs, files in os.walk(f):
      if len(files) > 1:
         if "alpha30.0" in root: fols.append(root)


X,Y = [],[]
JF,D,tRL,tLR,UR,UL,E1,E2 = [],[],[],[],[],[],[],[]
f_out = open(fout,'w')
f_out.write('#d   e   JF   D   tRL   tLR   UR   UL   E1   E2\n')
for fol in tqdm(fols):
   try: S = ex.Spectrum(fol)
   except FileNotFoundError: continue
   S.analyze_ingap()
   if hasattr(S, 'warning'):
      print(fol)
      continue
   d,e = get_d_e(fol)
   X.append(e)
   Y.append(d)
   f_out.write(str(d)+'   '+str(e)+'   ')
   jf,d,trl,tlr,ur,ul,e1,e2 = S.get_blue_parameters()
   JF.append( jf )
   D.append(d)
   tRL.append( trl )
   tLR.append( tlr )
   UR.append( ur )
   UL.append( ul )
   E1.append( e1 )
   E2.append( e2 )
   f_out.write(str(jf)+'   '+str(d)+'   '+str(trl)+'   '+str(tlr)+'   ')
   f_out.write(str(ur)+'   '+str(ul)+'   '+str(e1)+'   '+str(e2)+'\n')
   f_out.flush()
f_out.close()

import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.scatter(X,Y,c=JF)
ax.set_ylim([2,100])
ax.set_xlim([-0.2,0.2])
plt.show()

#H = ex.blue(JF,D,tRL,tLR,UR,UL,e1,e2)
#es,v = np.linalg.eigh(H)
#v = v.transpose()
#
#for i in range(len(es)):
#   print(i,es[i],np.round(v[i,:],4))
#import matplotlib.pyplot as plt
#plt.close('all')
#fig, ax = plt.subplots()
#ax.plot(es,'o')
#ax.grid()
#plt.show()

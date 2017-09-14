#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import mygreen_tools as green
import numpy as np
import os


def get_DOS(Emin, Emax, vintra, h, path_slf='/tmp', nE=101, use_all=True,
                                               fol='./', delta=0.01,add0=True):
   de = 0.1* (Emax - Emin)
   E = np.linspace(int(Emin-de),int(Emax+de),nE)
   if add0:
      nn = int(2/delta)+1
      E = np.sort(np.append(E,np.linspace(-1,1,nn)))
   if use_all:
      files = os.popen('ls -1 %s*.npy'%(path_slf)).read().splitlines()
      es = [x.split('/')[-1].replace('.npy','') for x in files]
      es = [x.split('/')[-1].replace('self','') for x in es]
      es = np.array(sorted([float(x) for x in es]))
      es = es[(es>Emin) & (es<Emax)]
      E = np.sort(np.append(E,es))
   DOSd,DOSp = [],[]
   f = open(fol+'dfct.dos','w')
   g = open(fol+'pris.dos','w')
   for e in E:
      Gd = green.green_function(e,vintra,h,path_selfes=path_slf,delta=delta)
      Gp = green.green_function(e,h.intra,h,path_selfes=path_slf,delta=delta)
      dd = -Gd.trace()[0,0].imag/np.pi
      dp = -Gp.trace()[0,0].imag/np.pi
      f.write(str(e)+'   ')
      g.write(str(e)+'   ')
      vd = -Gd.diagonal().imag/np.pi
      for ix in range(max(vd.shape)):
         f.write(str(vd[0,ix])+'   ')
      f.write('\n')
      vp = -Gp.diagonal().imag/np.pi
      for ix in range(max(vp.shape)):
         g.write(str(vp[0,ix])+'   ')
      g.write('\n')
      DOSd.append(dd)
      DOSp.append(dp)
   f.close()
   g.close()
   return E,DOSp,DOSd


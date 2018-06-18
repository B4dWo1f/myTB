#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import numpy as np
from mayavi import mlab
from read import xyz
import xcrysden

fname = '../../cells/ac_n1_l1.xyz'

ats,pos,latt,sub = xyz(fname)

def plot_crystal(ats,pos,latt):
   sizefig=(800,540)
   mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=sizefig)
   mlab.clf()
   difats = np.unique(ats)
   for at in difats:
      col = xcrysden.atom_color[at]
      siz = xcrysden.atom_size[at]
      Ps = pos[ats==at]
      X = Ps[:,0]
      Y = Ps[:,1]
      Z = Ps[:,2]
      mlab.points3d(X, Y, Z,\
              scale_factor=siz,\
              resolution=20,\
              opacity=0.99,\
              color=col,\
              scale_mode='none')
   for r in latt:
      for at in difats:
         col = xcrysden.atom_color[at]
         siz = xcrysden.atom_size[at]
         Ps = pos[ats==at] + r
         X = Ps[:,0]
         Y = Ps[:,1]
         Z = Ps[:,2]
         mlab.points3d(X, Y, Z,\
                 scale_factor=siz,\
                 resolution=20,\
                 opacity=0.25,\
                 color=col,\
                 scale_mode='none')
   #pos = (-0.7,0,0)
   #com = (1.4,0,0)
   #mlab.quiver3d(*pos,*com,color=(1,1,1),mode='cylinder')
   mlab.show()

plot_crystal(ats,pos,latt)

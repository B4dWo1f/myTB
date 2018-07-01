#!/usr/bin/python3
# -*- coding: UTF-8 -*-

# Standard libraries
# Third party libraries
import numpy as np
from mayavi import mlab
from tqdm import tqdm
# User defined libraries
import numeric as num


def lDOS(centers,heights,c=0.3,norm=True,show=False,fortran=True):
   """
     This function returns a XYZ grid with a gaussian with a given height in
     each center
   """
   def calc_lDOS(x,y,centers,heights,c=0.3):
      def gaussian3d(x,y,r0,a=1,c=0.3):
         r = np.sqrt((x-r0[0])**2 + (y-r0[1])**2)
         return a*np.exp(-(r*r/(2*c*c)))
      Z = 0.0*x
      #for i in tqdm(range(len(centers))):
      for i in range(len(centers)):
         cent,high = centers[i], heights[i]
         Z += gaussian3d(x,y,cent,high)
      return Z
   
   X = centers[:,0]
   Y = centers[:,1]
   d = 3*c
   xv, yv = np.mgrid[np.min(X)-d:np.max(X)+d:0.1, np.min(Y)-d:np.max(Y)+d:0.1]
   
   if fortran: zv = num.ldos(xv,yv,centers,heights,0.3)
   else: zv = calc_lDOS(xv,yv,centers,heights)
   if norm: zv = zv/np.max(zv)
   
   
   if show:
      mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0)) #, size=(1600,1080))
      mlab.clf()
      f = 10 #np.mean( [np.max(xv)-np.min(xv), np.max(yv)-np.min(yv)] )/10
      mlab.surf(xv,yv,f*zv)
      #mlab.axes(z_axis_visibility=False)
      # ==========
      #    Ejes
      # ==========
      def ejes(n=3):
         mlab.plot3d([0,n],[0,0],[0,0])
         mlab.plot3d([0,0],[0,n],[0,0])
         mlab.plot3d([0,0],[0,0],[0,n])
         mlab.text3d(n,-0.3,0,'X',scale=0.2)
         mlab.text3d(0,n,0,'Y',scale=0.2)
         mlab.text3d(0,0,n-0.15,'Z',scale=0.2)
      ejes()
      mlab.show()

   return xv,yv,zv

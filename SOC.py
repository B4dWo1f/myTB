#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import numpy as np

def ylm2xyz_l2():
   """Return the matrix that converts the cartesian into spherical harmonics"""
   m = np.matrix([[0.0j for i in range(5)] for j in range(5)])
   r2 = np.sqrt(2.)
   m[2,0] =  1. # dz2
   m[1,1] =  1./r2 # dxz
   m[3,1] = -1./r2 # dxz
   m[1,2] =  1j/r2 # dyz
   m[3,2] =  1j/r2 # dyz
   m[0,3] =  1j/r2 # dxy
   m[4,3] = -1j/r2 # dxy
   m[0,4] =  1./r2 # dx2y2
   m[4,4] =  1./r2 # dx2y2
#   m = m.H # no fucking idea if there is a transpose...
#   m = m.T # inverse
   #from increase_hilbert import spinful
   m = spinful(m) # with spin degree of freedom
   return m # return change of bassi matrix


def ylm2xyz_l1():
   """Return the matrix that converts the cartesian into spherical harmonics"""
   m = np.matrix([[0.0j for i in range(3)] for j in range(3)])
   s2 = np.sqrt(2.)
   m[1,2] = 1. # pz
   m[0,0] = 1./s2 # dxz
   m[2,0] = -1./s2 # dxz
   m[0,1] = 1j/s2 # dyz
   m[2,1] = 1j/s2 # dyz
   #from increase_hilbert import spinful
   m = spinful(m) # with spin degree of freedom
   return m

def m2spin(matin,matin2=[]):
  n=len(matin)
  from numpy import matrix
  if len(matin2)==0:
    matin2=matrix([[0.0j for i in range(n)] for j in range(n)])
  matout=matrix([[0.0j for i in range(2*n)] for j in range(2*n)])
  for i in range(n):
    for j in range(n):
      matout[2*i,2*j]=matin[i,j].copy()
      matout[2*i+1,2*j+1]=matin2[i,j].copy()
  return matout


def spinful(m):
  """ Return a spinful hamiltonian"""
  return m2spin(m,matin2=m)


def soc_l(l):
   """
     Calculate the spin orbit coupling in a basis of spherical harmonics
     and up,down, up,down, .....
   """
   nm = 2*l + 1 # number of components
   zero = np.matrix([[0.0j for i in range(2*nm)] for j in range(2*nm)])
   # initialize matrices
   lz = zero.copy()
   lx = zero.copy()
   ly = zero.copy()
   lm = zero.copy()
   lp = zero.copy()
   # create l+ and l- and lz
   for m in range(-l,l): # loop over m components
      val = np.sqrt((l-m)*(l+m+1)) # value of the cupling
      im = m + l
      lp[2*(im+1),2*im] = val # up channel
      lp[2*(im+1)+1,2*im+1] = val # down channel
   for m in range(-l,l+1): # loop over m components
      im = m + l
      lz[2*im,2*im] = m # value of lz, up channel
      lz[2*im+1,2*im+1] = m # value of lz, down channel
   lm = lp.H # adjoint
   lx = (lp + lm) /2.
   ly = -1j*(lp - lm) /2.
   # create spin matrices
   sz = zero.copy()
   sx = zero.copy()
   sy = zero.copy()
   for m in range(-l,l+1): # loop over m components
      im = m + l
      sx[2*im,2*im+1] = 1.0
      sx[2*im+1,2*im] = 1.0
      sy[2*im,2*im+1] = -1j
      sy[2*im+1,2*im] = 1j
      sz[2*im,2*im] = 1.
      sz[2*im+1,2*im+1] = -1.
   # check that the matrix is fine
   sx = sx/2.
   sy = sy/2.
   sz = sz/2.
   #if True:
   #  comm_zero(sx,lx)
   #  comm_zero(sx,ly)
   #  comm_zero(sx,lz)
   #  comm_zero(sy,ly)
   #  comm_zero(sz,lz)
   #  comm_angular(sx,sy,sz)
   #  comm_angular(sy,sz,sx)
   #  comm_angular(sz,sx,sy)
   #  comm_angular(lx,ly,lz)
   #  comm_angular(ly,lz,lx)
   #  comm_angular(lz,lx,ly)
   ls = lx*sx + ly*sy + lz*sz  # SOC matrix
   import scipy.linalg as lg
   from scipy.sparse import csc_matrix as csc
   ## rotate to cartesian orbitals
   # Choose rotation matrix
   if l==0: R = np.matrix([[0,0],[0,0]],dtype=complex)
   elif l==1: R = ylm2xyz_l1()
   elif l==2: R = ylm2xyz_l2()
   # Rotation
   ls = R.H * ls * R
   return ls

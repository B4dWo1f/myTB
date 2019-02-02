#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import numpy as np
from scipy.sparse import coo_matrix, csc_matrix
import logging
import log_help
LG = logging.getLogger(__name__)   # Logger for this module


def m2spin(M,Delt=0.0):
   """
     Duplicate a matrix like this:
                a+Delt   0    b+Delt   0
     a b  --->    0    a-Delt   0    b-Delt
     c d        c+Delt   0    d+Delt   0
                  0    c-Delt   0    d-Delt
     only for square matrices by now
   """
   try: n,m = M.shape
   except: n,m = len(M),0
   if isinstance(M,coo_matrix):
      Rn,Cn,Dn = [],[],[]
      for r,c,d in zip(M.row,M.col,M.data):
         ## Rows
         Rn.append(2*r)
         Rn.append(2*r+1)
         ## Cols
         Cn.append(2*c)
         Cn.append(2*c+1)
         ## Data
         try:
            Dn.append(d-Delt)
            Dn.append(d+Delt)
         except TypeError:    # In case of arrays of other type (str,...)
            Dn.append(d)      #
            Dn.append(d)      #
      return csc_matrix((Dn, (Rn, Cn)), shape=(2*n, 2*m),dtype=M.dtype)
   elif isinstance(M,np.matrix):
      LG.warning('Using dense matrix in m2spin (mat: %s,%s)'%(n,m))
      keep_type = 0.0*M[0,0]  # Preserve type (int, float, complex...)
      matout = np.matrix([[keep_type for i in range(2*n)] for j in range(2*n)])
      matout = np.matrix(np.zeros(M.shape),dtype=M.dtype)
      for i in range(n):
         for j in range(n):
            try:
               matout[2*i+1,2*j+1] = M[i,j]-Delt
               matout[2*i,2*j] = M[i,j]+Delt
            except TypeError:
               matout[2*i,2*j] = M[i,j]
               matout[2*i+1,2*j+1] = M[i,j]
      return csc_matrix(matout)
   elif isinstance(M,list):
      return np.repeat(M,2)
   else:
      if isinstance(M,np.ndarray): return np.repeat(M,2)  #XXX Warning
         #if M.dtype == '<U2': return np.repeat(M,2)
      LG.debug('Unkown matrix type in m2spin, trying to convert to coo_matrix')
      return m2spin(coo_matrix(M),Delt=Delt)


def isinb2een(x,a,b):
   try: return np.array([isinb2een(ix,a,b) for ix in x])
   except:
      if (min(a.real,b.real) <= x.real <= max(a.real,b.real)) and \
         (min(a.imag,b.imag) <= x.imag <= max(a.imag,b.imag)):
         c = a-b
         m0 = abs(c.imag/c.real)
         c = x-a
         m1 = abs(c.imag/c.real)
         if m0 == m1: return True
      else: return False

def dens2band(H):
   M = []
   for j in range(H.shape[1]):
      aux = []
      for im in range(H.shape[0]-j):
         aux.append(H[im,im+j])
      for _ in range(H.shape[0]-j,H.shape[1]):
         aux = [0] + aux
      M.append(aux)
   M.reverse()
   M = np.array(M)
   return M


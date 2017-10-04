#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import numpy as np
from scipy.sparse import coo_matrix



def m2spin(M,Delt=0.0): #,matin2=None):
   """
     Duplicate a matrix like this:
                a+Delt   0    b+Delt   0
     a b  --->    0    a-Delt   0    b-Delt
     c d        c+Delt   0    d+Delt   0
                  0    c-Delt   0    d-Delt
     only for square matrices by now
   """
   n,m = M.shape
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
         Dn.append(d-Delt)
         Dn.append(d+Delt)
      return csr_matrix((Dn, (Rn, Cn)), shape=(2*n, 2*m))
   elif isinstance(M,np.matrix):
      print('WARNING!!! Dense!!!')
      keep_type = 0.0*M[0,0]  # Preserve type (int, float, complex...)
      matout = np.matrix([[keep_type for i in range(2*n)] for j in range(2*n)])
      matout = np.matrix(np.zeros(M.shape),dtype=M.dtype)
      for i in range(n):
         for j in range(n):
            matout[2*i,2*j] = M[i,j]+Delt
            matout[2*i+1,2*j+1] = M[i,j]-Delt
      return csr_matrix(matout)
   else: m2spin(coo_matrix(M))


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

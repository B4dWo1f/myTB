#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import numpy as np


def m2spin(matin,Delt=0.0): #,matin2=None):
   """
     Duplicate a matrix like this:
                a+Delt   0    b+Delt   0
     a b  --->    0    a-Delt   0    b-Delt
     c d        c+Delt   0    d+Delt   0
                  0    c-Delt   0    d-Delt
   """
   n = matin.shape[0]
   keep_type = 0*matin[0,0]  # Preserve type (int, float, complex...)
   matout = np.matrix([[keep_type for i in range(2*n)] for j in range(2*n)])
   for i in range(n):
      for j in range(n):
         matout[2*i,2*j] = matin[i,j]+Delt
         matout[2*i+1,2*j+1] = matin[i,j]-Delt
   return matout

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

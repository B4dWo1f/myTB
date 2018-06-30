#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import numpy as np
import algebra as alg
from scipy.sparse import coo_matrix,bmat,csc_matrix
import logging
LG = logging.getLogger(__name__) # Logger for this module


def spin(base,dire=(0,0,1)):
   N = 0                   #  shouldn't this be the ham dimension?
   for E in base:          #
      N += len(E.onsite)   #
   c = [2*i+1 for i in range(N)]
   r = [2*i for i in range(N)]
   row = np.array(r+c)
   col = np.array(c+r)
   ## Sx
   data = np.array([1 for _ in row])
   Sx = coo_matrix((data, (row, col)), shape=(2*N, 2*N))
   ## Sy
   data = np.array([1j*(-1)**(i+1) for i in row])
   Sy = coo_matrix((data, (row, col)), shape=(2*N, 2*N))
   ## Sz
   c = [i for i in range(2*N)]
   r = [i for i in range(2*N)]
   row = np.array(r)
   col = np.array(c)
   data = np.array([(-1)**(i) for i in row])
   Sz = coo_matrix((data, (row, col)), shape=(2*N, 2*N))
   return dire[0]*Sx + dire[1]*Sy + dire[2]*Sz


def sublattice(base):   #,Subs=['A']):
   """
     This function returns a matrix with 1's in the sublattice "sublattice"
     positions
   """
   n = base.ndim
   row,dat = [],[]
   for i in range(len(base.SUBS)):
      row.append(i)
      dat.append(base.SUBS[i])
   return csc_matrix((dat,(row,row)),shape=(n,n))


def orbital(base,Orbs=['pz']):
   """
     This function returns a matrix with 1's in the place of the provided
     orbitals
   """
   if isinstance(Orbs,str): Orbs = [Orbs]
   elif isinstance(Orbs,list): pass
   else: LG.error('Wrong input in orbital operator')
   diag_name = []
   for E in base:
      for o in E.orbitals:
         diag_name.append(o)
   diag = [0 for _ in diag_name]
   for i in range(len(diag_name)):
      if diag_name[i] in Orbs: diag[i] = 1
   if base.DOspin: return alg.m2spin(coo_matrix(np.diag(diag)))
   else: return coo_matrix(np.diag(diag))


def atom(base,Ats=[1]):
   """
     This function returns a matrix with identity matrices in the places of
     the provided atoms (provided by order in the basis)
   """
   if isinstance(Ats,int): Ats = [Ats]
   elif isinstance(Ats,list): pass
   else: LG.error('Wrong input in atom operator')
   aux = [[None for _ in base] for _ in base]
   for E in base:
      dim = len(E.onsite)
      if E.place in Ats: aux[E.place][E.place] = np.identity(dim,dtype=int)
      else: aux[E.place][E.place] = np.zeros((dim,dim),dtype=int)
   if base.DOspin: return alg.m2spin(bmat(aux))
   else: return bmat(aux)


def layer(base):
   """
     This function returns a matrix with the layer of each element in the
     diagonal
   """
   aux = [[None for _ in base] for _ in base]
   for E in base:
      N = len(E.onsite)
      aux[E.place][E.place] = E.layer * np.identity(N)
   if base.DOspin: return alg.m2spin(bmat(aux))
   else: return bmat(aux)


def dist(base,r0=np.array([0.,0.,0.])):
   """
     Returns a matrix with the distance of each atom to a given point in the
     diagonal
   """
   aux = [[None for _ in base] for _ in base]
   for E in base:
      r = np.linalg.norm( E.position - r0 )
      N = len(E.onsite)
      aux[E.place][E.place] = r * np.identity(N)
   if base.DOspin: return alg.m2spin(bmat(aux))
   else: return bmat(aux)


def position(base,coor=2):
   """
     Returns a matrix with the given coordinate of each atom in the diagonal
   """
   aux = [[None for _ in base] for _ in base]
   for E in base:
      r = E.position
      N = len(E.onsite)
      aux[E.place][E.place] = r[coor] * np.identity(N)
   if base.DOspin: return alg.m2spin(bmat(aux))
   else: return bmat(aux)

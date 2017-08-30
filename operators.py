#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import numpy as np
from scipy.sparse import coo_matrix,bmat
import logging

LG = logging.getLogger(__name__) # Logger for this module
### File Handler
#LG = logging.getLogger('operators') # Logger for this module
#fh = logging.FileHandler('operators.log',mode='w')
#fmt = logging.Formatter('%(asctime)s %(name)s:%(levelname)s - %(message)s')
#fh.setFormatter(fmt)
#fh.setLevel(logging.DEBUG)
#LG.addHandler(fh)
#LG.setLevel(logging.DEBUG)
### Screen Handler
##sh = logging.StreamHandler()
##fmt = logging.Formatter('%(name)s -%(levelname)s- %(message)s')
##sh.setFormatter(fmt)
##sh.setLevel(logging.WARNING)
##LG.addHandler(sh)


def spin(base,dire=(0,0,1)):
   def pauli_matrix(n):
      """
        Builds 3 nxn matrices with the each of the Pauli Matrices in the
        diagonal
      """
      Sx = np.matrix([[0,1],[1,0]],dtype=complex)
      Sy = np.matrix([[0,-1j],[1j,0]],dtype=complex)
      Sz = np.matrix([[1,0],[0,-1]],dtype=complex)
      auxX = [[None for _ in range(n)] for _ in range(n)]
      auxY = [[None for _ in range(n)] for _ in range(n)]
      auxZ = [[None for _ in range(n)] for _ in range(n)]
      for i in range(n):
         auxX[i][i] = Sx
         auxY[i][i] = Sy
         auxZ[i][i] = Sz
      return bmat(auxX), bmat(auxY), bmat(auxZ)
   N = 0
   for E in base:
      N += len(E.onsite)
   sig = pauli_matrix(N)
   return dire[0]*sig[0] + dire[1]*sig[1] + dire[2]*sig[2]


def sublattice(base,Subs=['A']):
   """
     This function returns a matrix with 1's in the sublattice "sublattice"
   positions
   """
   if isinstance(Subs,str): Zorb = [Subs]
   elif isinstance(Subs,list): pass
   else: LG.error('Wrong input in sublattice operator')
   aux = [[None for _ in base] for _ in base]
   for E in base:
      i = E.place
      dim = len(E.onsite)
      if E.sublattice in Subs: aux[i][i] = np.identity(dim,dtype=int)
      else: aux[i][i] = np.zeros((dim,dim),dtype=int)
   return bmat(aux)


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
   return coo_matrix(np.diag(diag))


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
   return bmat(aux)


def layer(base):
   """
     This function returns a matrix with the layer of each element in the
     diagonal
   """
   aux = [[None for _ in base] for _ in base]
   for E in base:
      N = len(E.onsite)
      aux[E.place][E.place] = E.layer * np.identity(N)
   return bmat(aux)


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
   return bmat(aux)


def position(base,coor=2):
   """
     Returns a matrix with the given coordinate of each atom in the diagonal
   """
   aux = [[None for _ in base] for _ in base]
   for E in base:
      r = E.position
      N = len(E.onsite)
      aux[E.place][E.place] = r[coor] * np.identity(N)
   return bmat(aux)
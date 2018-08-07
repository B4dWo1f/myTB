#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import matplotlib.pyplot as plt
import numpy as np

def H0(e0,t1):
   """ Intra-cell """
   return np.matrix([[e0,t1,t1],\
                     [t1,e0,t1],\
                     [t1,t1,e0]])

def H1(t1,t2,t3,t4=0,t5=0):
   """ Hopping with cell at a1 """
   return np.matrix([[t3,t2,t1],\
                     [t4,t3,t2],\
                     [t5,t4,t3]])

def H2(t1,t2,t3,t4=0,t5=0):
   """ Hopping with cell at a2 """
   return np.matrix([[t3,t4,t2],\
                     [t2,t3,t1],\
                     [t4,t5,t3]])

def H1m2(t1,t2,t3,t4=0,t5=0):
   """ Hopping with cell at a1-a2 """
   return np.matrix([[t3,t1,t2],\
                     [t5,t3,t4],\
                     [t4,t2,t3]])


def hamil(k,a1,a2,e0,t1,t2,t3,t4=0,t5=0):
   """ k-dependent hamiltonian """
   h0 = H0(e0,t1)
   h1 = H1(t1,t2,t3,t4,t5)
   h2 = H2(t1,t2,t3,t4,t5)
   h1m2 = H1m2(t1,t2,t3,t4,t5)
   # ** There is no hopping to cell a1+a2
   return h0 +\
          np.exp(-1j*np.dot(k,a1))*h1 + np.exp(1j*np.dot(k,a1))*h1.H+\
          np.exp(-1j*np.dot(k,a2))*h2 + np.exp(1j*np.dot(k,a2))*h2.H+\
          np.exp(-1j*np.dot(k,a1-a2))*h1m2 + np.exp(1j*np.dot(k,a1-a2))*h1m2.H


def reciprocal(latt_vec):
   """
    Returns the reciprocal vectors given the lattice vectors
    TO_-DO: Check 3D case
   """
   if len(latt_vec) == 0: # Islands, no Translation symmetry
      return []
   else:
      if len(latt_vec) == 1:
         a1 = latt_vec[0]
         # Define more direct vectors to use the same formula
         a2 = np.cross( a1,np.array([rand(),rand(),rand()]) )
         a2 = a2*np.linalg.norm(a1)/np.linalg.norm(a2)
         a3 = np.cross(a1,a2)
         a3 = a3*np.linalg.norm(a1)/np.linalg.norm(a3)
      elif len(latt_vec) == 2:
         a1 = latt_vec[0]
         a2 = latt_vec[1]
         a3 = np.cross(a1,a2)
         a3 = a3*np.linalg.norm(a1)/np.linalg.norm(a3)
      elif len(latt_vec) == 3:
         a1 = latt_vec[0]
         a2 = latt_vec[1]
         a3 = latt_vec[2]
      b1 = 2*np.pi*(np.cross(a2,a3)/(np.dot(a1,np.cross(a2,a3))))
      b2 = 2*np.pi*(np.cross(a3,a1)/(np.dot(a1,np.cross(a2,a3))))
      b3 = 2*np.pi*(np.cross(a1,a2)/(np.dot(a1,np.cross(a2,a3))))
      vecs = [b1,b2,b3]
      recip_vec = []
      for i in range(len(latt_vec)):   # returns only the necessary
         recip_vec.append(vecs[i])     # reciprocal vectors
      return recip_vec

def recorrido(points,nk):
   """
     Returns nk*len(points) points along the path defined by points in a
     N-dimensional space
    points: list of points between which to interpolate
    nk: may be an integer or a sequence of len(points)-1 integers
   """
   ## clean nk
   try: itr = iter(nk)
   except TypeError: nk = [nk for _ in range(len(points)-1)]
   else:
      if len(nk) < len(points)-1:
         msg = 'WARNING: incorrect number of points. '
         msg += 'Using %s for every interval'%(nk[0])
         LG.warning(msg)
         nk = [nk[0] for _ in range(len(points)-1)]
      elif len(nk) > len(points)-1: pass    # Report to log?
      else: pass    # Report to log?

   ## interpolate between points
   RECORRIDO = []
   ret = (1,False)   # don't skip last point
   lim = len(points)-1
   for ipunto in range(lim):
      N = nk[ipunto]
      if ipunto == lim-1: ret = (0,True)   # skip last point
      P1 = points[ipunto]
      P2 = points[ipunto+1]
      coors = []
      for idim in range(len(P1)): # loop over dimensionality
         coors.append(np.linspace(P1[idim],P2[idim],N+ret[0],endpoint=ret[1]))
      for p in zip(*coors):
         RECORRIDO.append(np.array(p))
   return RECORRIDO



def cell(N,a=1.4,buck=0.0,cent=True,show=False):
   """
   Returns the list of atoms, atomic positions, lattice vectors and sublattice
   of a Kagome lattice
   The Kagome lattice is not bipartite so no sublattice will be returned
   """
   if N != 1:
      print('WARNING: N!=1 not implemented yet. Using N=1 instead')
      N = 1
   ap = np.sqrt(3)/2.   # mathematical constant
   r3 = np.sqrt(3)
   brick = [np.array([-a/2,0,0]),
            np.array([ a/2,0,0]),
            np.array([0,r3*a/2,0])]
   e=1.
   vectors = [e*np.array([2*a,0,0]),
              e*np.array([a,r3*a,0])]
   pos = []
   for i in range(N):
      for j in range(N):
         for ir in range(len(brick)):
            r = brick[ir]
            p = r + i*vectors[0] + j*vectors[1]
            pos.append(p)
   ats = np.array(['C' for _ in pos])
   latt = vectors
   return ats,np.array(pos),latt,None   # No sublattice


################################################################################
if __name__ == '__main__':
   ## Atomic positions and lattice vectors
   _,pos,latt,_ = cell(1)
   if True:   # Plot lattice
      X0 = pos[:,0]
      Y0 = pos[:,1]
      X,Y = [],[]
      for v in latt:
         for p in pos:
            r = v+p
            X.append(r[0])
            Y.append(r[1])
            r = -v+p
            X.append(r[0])
            Y.append(r[1])
      fig, ax = plt.subplots()
      ax.scatter(X0,Y0)
      ax.scatter(X,Y)
      #plt.show()


   ## High-symmetry points to define the K-path
   # Reciprocal vecctors
   recip = reciprocal(latt)
   G = 0*recip[0] + 0*recip[1]
   K = (2*recip[0]+recip[1])/3.
   Kp = (recip[0]+2*recip[1])/3.
   M = (recip[0]+recip[1])/2.

   ## K-path
   points = [G,K,Kp,G]
   rec = recorrido(points,50)

   ## Bands for first-neighbors only
   e0,t1,t2,t3 = 0.,-1.,0.,0.
   X0,Y0 = [],[]
   for i in range(len(rec)):
      k = rec[i]
      H = hamil(k,latt[0],latt[1],e0,t1,t2,t3)
      for e in np.linalg.eigvalsh(H):
         X0.append(i)
         Y0.append(e)

   ## Bands for first and second neighbors
   e0,t1,t2,t3 = 0.,-1.,-0.3,0.1
   X1,Y1 = [],[]
   for i in range(len(rec)):
      k = rec[i]
      H = hamil(k,latt[0],latt[1],e0,t1,t2,t3)
      for e in np.linalg.eigvalsh(H):
         X1.append(i)
         Y1.append(e)


   fig, ax = plt.subplots()
   ax.scatter(X0,Y0)
   ax.scatter(X1,Y1)
   plt.show()

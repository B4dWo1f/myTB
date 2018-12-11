#!/usr/bin/python3
# -*- coding: UTF-8 -*-

"""
  All the cell-creating functions should return 4 items:
             return ats, pos, latt, subs
  ats: np.array, dtype=str. Contains the atomic elements of the unit cell
  pos: np.array, dtype=float. Contains the atomic positions of the unit cell
  latt: list containing N np.arrays(float). Contains the lattice vectors
  subs: np.array, dtype=int. Contains the corresponding sublattice.
  Extended XYZ format:
    Natoms
    [],[],[]  Lattice vectors
    E   x   y   z sublatt
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from itertools import product
import geometry as geo
## LOG
import logging
LG = logging.getLogger(__name__)

class UnitCell(object):
   def __init__(self,ats,pos,latt,subs=[],pasivate=False):
      self.ats = np.array(ats,str)
      self.pos = np.array(pos)
      self.latt = np.array(latt)
      self.sub = np.array(subs)
      #self.sub = np.array(subs,int)   #TODO check type
      if self.check() and len(self.ats) > 0: self.get_geo_info()
      else: pass
      if pasivate: self.pasivate()
   def get_geo_info(self):
      """ Calculate geometric attributes """
      if not self.check():
         LG.critical('Different number of atoms, positions and sublatties')
      self.center = np.mean(self.pos,axis=0)
      lmin = np.min(self.pos,axis=0)
      lmax = np.max(self.pos,axis=0)
      self.lims = [lmin,lmax]
   def pasivate(self):
      """
        Add Hidrogen atoms on the C3 missing positions
      """
      latt = []
      hs,subh = pasivate(self.pos,sub=self.sub)
      self.ats = np.append(self.ats,['H' for _ in hs])
      aux = []
      for i in range(pos.shape[0]):
         aux.append(pos[i,:])
      for i in range(len(hs)):
         aux.append(hs[i])
      self.pos = np.array(aux) #np.append(self.pos,np.array(hs))
      self.sub = np.append(self.sub,subh)
   def multilayer(self,lN):
      self.ats,self.pos,self.sub = multilayer(self.pos,self.ats,self.sub,N=lN)
      #self.center = np.mean(np.array(self.pos),axis=0)
      self.get_geo_info()
   def check(self):
      if len(self.ats) == len(self.pos) == len(self.sub): return True
      else: return False
   def to_xyz(self,fname=None):
      """
        Save the info about the UCell to a file in extended xyz format
      """
      #TODO use IO.write.xyz
      print(len(self.ats))
      print('')
      for i in range(len(ats)):
      #for a,r,s in zip(self.ats,self.pos,self.sub):
         a = self.ats[i]
         r = self.pos[i]
         s = self.sub[i]
         print('%s   %s   %s   %s   %s'%(a,r[0],r[1],r[2],s))
      if fname != None:
         f = open(fname,'w')
         f.write(str(len(self.ats))+'\n')
         ## Write lattice vectors
         for v in latt:
            f.write('[%s,%s,%s]'%(v[0],v[1],v[2]))
         f.write('\n')
         for a,r,s in zip(self.ats,self.pos,self.sub):
            f.write('%s   %s   %s   %s   %s\n'%(a,r[0],r[1],r[2],s))
         f.close()
   def from_xyz(self,fname,pasivate=False):
      """ Read Unit Cell information from extended xyz file """
      from IO.read import xyz
      self.ats,self.pos,self.latt,self.sub = xyz(fname)
      #self.center = np.mean(np.array(self.pos),axis=0)
      self.get_geo_info()
   def plot(self,fname=None):
      plot_cell(self.pos,self.latt,fname=fname)
   def lot_sublattice(self,fname=None):
      X = self.pos[:,0]
      Y = self.pos[:,1]
      S = self.sub
      fig, ax = plt.subplots()
      ax.scatter(X,Y,c=S,s=100,edgecolors='none')
      plt.show()
   def __str__(self):
      """ Overload of the string method for pretty printing the class """
      C = self.center
      L = self.lims
      names = ['X','Y','Z']
      msg = 'Unit cell with %s atoms\n'%(len(self.ats))
      msg += 'centered at: (%.2f,%.2f,%.2f) and with lims:\n'%(C[0],C[1],C[2])
      for i in range(len(names)):
         msg += names[i]+'   %.3f -- %.3f\n'%(L[0][i],L[1][i])
      return msg



def ribbon_armc(Nx,Ny,a=1.4,buck=0.0,cent=True,show=False):
   """
      The function returns 2 lists, one containing the positions of the atoms
     and the other containing the lattice vectors for an ARMCHAIR island.
     The parameters are the following:
       N: [int] number of repetitions of the brick
       a: [float] atomic distance
       buck: [folat] buckling of the atoms (introduced by sublattice)
       show: show a 2D-plot of the unit cell and lattice vectors
   """
   ap = np.sqrt(3)/2.   # mathematical constant
   ## Positions of a benzene
   b = buck/2.
   brick = [a*np.array([1.,0.,b]),
            a*np.array([1/2.,ap,-b]),
            a*np.array([-1/2.,ap,b]),
            a*np.array([-1.,0.,-b])]
   sublatt = [1,-1,1,-1]  # XXX check order
   
   vectors = [a*np.array([3.,0.,0.]),
              a*np.array([0.,2.*ap,0])]  # to expand the unit cell
              #a*np.array([0,3.*ap,0])]  # to expand the unit cell
   latt = [Nx*vectors[0]]

   cell_x,sub_aux = [],[]
   for i in range(Nx):
      for p,s in zip(brick,sublatt):
         cell_x.append(i*vectors[0]+p)
         sub_aux.append(s)
   pos,subs = [],[]
   for i in range(Ny):
      for j in range(len(cell_x)):
         p = cell_x[j]
         pos.append(i*vectors[1]+p)
         subs.append(sub_aux[j])
   ## Re-Center the unit cell
   if cent:
      X = [p[0] for p in pos]
      Y = [p[1] for p in pos]
      Z = [p[2] for p in pos]
      C = np.array( [np.mean(X),np.mean(Y),np.mean(Z)] )
      for i in range(len(pos)):
         pos[i] -= C
   ats = np.array(['C' for _ in pos])
   pos = np.array(pos)
   subs = np.array(subs)
   return ats,pos,latt,subs
   #return UnitCell(ats,pos,latt,subs=[])

def armchair(N,a=1.4,buck=0.0,show=False):
   """
      The function returns 2 lists, one containing the positions of the atoms
     and the other containing the lattice vectors for an ARMCHAIR island.
     The parameters are the following:
       N: [int] number of repetitions of the brick
       a: [float] atomic distance
       buck: [folat] buckling of the atoms (introduced by sublattice)
       show: show a 2D-plot of the unit cell and lattice vectors
   """
   ap = np.sqrt(3)/2.   # mathematical constant
   ## Positions of a benzene
   b = buck/2.
   brick = [a*np.array([1.,0.,b]),
            a*np.array([1/2.,ap,-b]),
            a*np.array([-1/2.,ap,b]),
            a*np.array([-1.,0.,-b]),
            a*np.array([-1/2.,-ap,b]),
            a*np.array([1/2.,-ap,-b])]
   sublatt = [1,-1,1,-1,1,-1]  # XXX check order
   
   vectors = [a*np.array([3/2.,3.*ap,0.]),
              a*np.array([3.,0.,0.])]
   latt = [(N+1)*vectors[0]+N*vectors[1],
           -N*vectors[0]+(2*N+1)*vectors[1]]
   
   ## Start combinations
   lista = range(-N,N+1)
   perms = [p for p in product(lista, repeat=2)]
   lim = N+1
   all_vecs = []  # all combis where to replicate the brick
   for p in perms:
      if abs(np.sum(p)) < lim:
         vec = np.array([0.,0.,0.])
         for i in range(len(p)):
            vec += p[i]*vectors[i]
         all_vecs.append(vec)
   pos,sub = [],[]
   for v in all_vecs:
      #for r in brick:
      for i in range(len(brick)):
         r = brick[i]
         s = sublatt[i]
         w = r+v
         pos.append(w)  # All the atomic positions
         sub.append(s)  # All the atomic positions
   ### Plot
   ats = np.array(['C' for _ in pos])
   pos = np.array(pos)
   sub = np.array(sub)
   return ats,pos,latt,sub
   #return UnitCell(ats,pos,latt,subs=[])


def zigzag(N,a=1.4,buck=0.0,show=False):
   """
      The function returns 2 lists, one containing the positions of the atoms
     and the other containing the lattice vectors for an ZIGZAG island.
     The parameters are the following:
       N: [int] number of repetitions of the brick
       a: [float] atomic distance
       buck: [folat] buckling of the atoms (introduced by sublattice)
       show: show a 2D-plot of the unit cell and lattice vectors
   """
   ap = np.sqrt(3)/2.   # mathematical constant
   ## Positions of a benzene
   b = buck/2.
   brick = [a*np.array([1.,0.,b]),
            a*np.array([1/2.,ap,-b]),
            a*np.array([-1/2.,ap,b]),
            a*np.array([-1.,0.,-b]),
            a*np.array([-1/2.,-ap,b]),
            a*np.array([1/2.,-ap,-b])]
   
   vectors = [a*np.array([3/2.,-ap,0]),
              a*np.array([3/2.,ap,0])]
   latt = [(N+1)*(vectors[0]+vectors[1]),
           (N+1)*(-vectors[0]+2*vectors[1])]
   subs = [1,-1,1,-1,1,-1]

   ## Start combinations
   lista = range(-N,N+1)
   perms = [p for p in product(lista, repeat=2)]
   lim = N+1
   all_vecs = []  # all combis where to replicate the brick
   for p in perms:
      if abs(np.sum(p)) < lim:
         vec = np.array([0.,0.,0.])
         for i in range(len(p)):
            vec += p[i]*vectors[i]
         all_vecs.append(vec)
   pos,sub = [],[]   #TODO sublattice
   for v in all_vecs:
      for ir in range(len(brick)):
         r = brick[ir]
         s = subs[ir]
         w = r+v
         # Avoid repited atoms
         if not vec_in_list(w,pos):
            pos.append(w)
            sub.append(s)
   ### Plot
   if show: plot_cell(pos,latt,tit='ZigZag Cell (%s)'%(N))
   ats = np.array(['C' for _ in pos])
   pos = np.array(pos)
   subs = np.array(subs)
   return ats,pos,latt,sub
   #return UnitCell(ats,pos,latt,subs=[])


def kagome(N,a=1.4,buck=0.0,cent=True,show=False):
   """
   A Kagome lattice is not bipartite so no sublattice will be returned
   """
   #if N != 1:
   #   print('WARNING: N!=1 not implemented yet. Using N=1 instead')
   #   N = 1
   r3 = np.sqrt(3)
   ap = r3/2.   # mathematical constant
   brick = [np.array([-a/2,0,0]),
            np.array([ a/2,0,0]),
            np.array([0,r3*a/2,0])]

   vectors = [a*np.array([2,0,0]),
              a*np.array([1,r3,0])]
   latt = [N*vectors[0],
           N*vectors[1]]
   subs = [-1,0,1]

   pos,sub = [],[]   #TODO check sublattice
   for i in range(N):
      for j in range(N):
         for ir in range(len(brick)):
            r = brick[ir]
            p = r + i*vectors[0] + j*vectors[1]
            pos.append(p)
            sub.append(subs[ir])
   ats = np.array(['C' for _ in pos])
   return ats,pos,latt,sub

def simple(N,a=1.4,buck=0.0,cent=True,show=False):
   """
      The function returns 2 lists, one containing the positions of the atoms
     and the other containing the lattice vectors for the simplest graphene
     super-cell.
     The parameters are the following:
       N: [int] number of repetitions of the brick
       a: [float] atomic distance
       buck: [folat] buckling of the atoms (introduced by sublattice)
       cent: [boolean] center the unit cell at (0,0,0)
       show: show a 2D-plot of the unit cell and lattice vectors
   """
   if N == 0:
      print('WARNING: N=0 is ill-defined. Using N=1 instead')
      N = 1
   ap = np.sqrt(3)/2.   # mathematical constant
   b = buck/2.
   brick = [a*np.array([-1/2.,0.,-b]),
            a*np.array([ 1/2.,0.,b]) ]
   vectors = [a*np.array([3/2.,-ap,0.]),
              a*np.array([3/2., ap,0.])]
   latt = [N*vectors[0],
           N*vectors[1]]
   sublatt = [1,-1]

   pos,sub = [],[]
   for i in range(N):
      for j in range(N):
         for ir in range(len(brick)):
            r = brick[ir]
            p = r + i*vectors[0] + j*vectors[1]
            s = sublatt[ir]
            pos.append(p)
            sub.append(s)
   ## Re-Center the unit cell
   if cent:
      C = np.mean(pos,axis=0)
      for i in range(len(pos)):
         pos[i] -= C
   if show: plot_cell(pos,latt,tit='Simple Cell %sx%s'%(N,N))
   ats = np.array(['C' for _ in pos])
   pos = np.array(pos)
   sub = np.array(sub)
   return ats,pos,latt,sub
   #return UnitCell(ats,pos,latt,subs=[])

def zigzag_triangle(N,a=1.4,buck=0.0,show=False):
   ap = np.sqrt(3)/2.   # mathematical constant
   ## Positions of a benzene
   b = buck/2.
   brick = [a*np.array([1.,0.,b]),
            a*np.array([1/2.0,ap,-b]),
            a*np.array([-1/2.,ap,b]),
            a*np.array([-1.,0.,-b]),
            a*np.array([-1/2.,-ap,b]),
            a*np.array([1/2.,-ap,-b])]
   vectors = [a*np.array([3/2.,-ap,0.]),
              a*np.array([3/2.,ap,0.])]
   subs = [1,-1,1,-1,1,-1]
   ## Start combinations
   lista = range(N+1)
   perms = [p for p in product(lista, repeat=2)]
   lim = N+1
   all_vecs = []  # all combis where to replicate the brick
   for p in perms:
      if abs(np.sum(p)) < lim:
         vec = np.array([0.,0.,0.])
         for i in range(len(p)):
            vec += p[i]*vectors[i]
         all_vecs.append(vec)
   pos,sub = [],[]
   for v in all_vecs:
      for ir in range(len(brick)):
         w = brick[ir]+v
         s = subs[ir]
         # Avoid repited atoms
         if not vec_in_list(w,pos):
            pos.append(w)
            sub.append(s)
   if show: plot_cell(pos,tit='Triangular Island (%s)'%(N))
   ats = np.array(['C' for _ in pos])
   pos = np.array(pos)
   subs = np.array(subs)
   return ats,pos,[],sub
   #return UnitCell(ats,pos,latt,subs=[])

def mullen(Nx,Ny=4,pas=False):
   ##XXX sublattice is wrong!!!!
   ats,pos,_,sub = ribbon_armc(Nx,Ny)
   #print(len(ats))
   #print('')
   #for a,p,s in zip(ats,pos,sub):
   #   print(a,p[0],p[1],p[2],s)
   #exit()
   ## Remove extra atoms
   pos = np.array(pos)
   sub = np.array(sub)
   X = pos[:,0]
   Y = pos[:,1]
   Z = pos[:,2]

   X = X[Y!=np.min(Y)]
   Z = Z[Y!=np.min(Y)]
   sub = list(sub[Y!=np.min(Y)])
   Y = Y[Y!=np.min(Y)]
   ## Sotre valid atoms
   pos = [np.array((x,y,z)) for x,y,z in zip(X,Y,Z)]
   #pos,latt,sub = mullen(Nx,Ny)
   #hs,subh = pasivate(pos,sub=sub)
   #pos += hs
   #ats += ['H' for _ in hs]
   #sub += subh
   ats = np.array(['C' for _ in pos])
   pos = np.array(pos)
   sub = np.array(sub)
   return ats,pos,[],sub
   #return UnitCell(ats,pos,latt,subs=[])




def vec_in_list(v,l,eps=1e-9):
   """ Returns True if vector v is in the list of vectors l (also in util.py)"""
   for x in l:
      if np.linalg.norm(x-v) < eps: return True
   return False

def plot_cell(pos,latt=[],tit=None,fname=None,show=True):
   """
     Plots the unit cell, lattice vectors, and first neighbouring unit cells.
   """
   fig = plt.figure() #figsize=(20,10))
   gs = gridspec.GridSpec(1, 1)
   fig.subplots_adjust(wspace=0.25,hspace=0.0)
   ax = plt.subplot(gs[0])  # Original plot

   ## Plot Unit cell
   X,Y = [],[]
   for r in pos:
      X.append(r[0])
      Y.append(r[1])
   ax.scatter(X,Y,c='k',s=100,edgecolors='none')
   ## Plot neighbouring cells
   if len(latt) > 0:
      cs = ['b','r','g','y','c','m']
      v_norm = np.mean([np.linalg.norm(v) for v in latt])
      ## Empiric size of the arrow head
      hw = v_norm * 0.17/3.46410161514
      hl = hw * 0.3/0.2
      i = 0
      for v in latt:
         vn = v/np.linalg.norm(v) # normalized vector
         vv = (np.linalg.norm(v)-hl)* vn #vector minus the length of the arrow
         X,Y = [],[]
         for r in pos:
            w = r+v
            X.append(w[0])
            Y.append(w[1])
         ax.scatter(X,Y,c=cs[i],s=100,edgecolors='none')
         ax.arrow(0,0,vv[0],vv[1],head_width=hw,head_length=hl,fc='b', ec='b')
         ax.text(v[0],v[1], r'$\vec{a}_{%s}$'%(i+1), fontsize=20,
                             bbox={'facecolor':'white', 'alpha':0.7, 'pad':5})
         i+=1
      ## Extra cells  XXX Error for 1D
      if len(latt) == 1: latt.append(np.array([0.,0.,0.])) # XXX Shame on you!!!
      latt2 = []
      for v in latt:
         latt2.append(-v)
      latt2.append(latt[0]-latt[1])
      latt2.append(-latt[0]+latt[1])
      for v in latt2:
         X,Y = [],[]
         for r in pos:
            w = r+v
            X.append(w[0])
            Y.append(w[1])
         ax.scatter(X,Y,c=cs[i],s=90,edgecolors='none')
         X = [0,v[0]]
         Y = [0,v[1]]
         ax.plot(X,Y,'b--')
         i+=1

   if tit != None: ax.set_title(tit)
   ax.axis('equal')
   ax.grid()
   if fname != None: fig.savefig(fname)
   if show: plt.show()

def multilayer(pos,ats,sub=[],N=2,vec=np.array([1.4,0,1.4])):
   """ Generates the positions for a multilayer ABC... """
   new_ats, new_pos, new_sub = [], [], []
   rs = [i*vec for i in range(N)]
   for r in rs:
      for j in range(len(pos)):
         new_ats.append(ats[j])
         new_pos.append(pos[j]+r)
         new_sub.append(sub[j])
         #try: new_sub.append(sub[j])
         #except: pass
   return new_ats, new_pos, new_sub


import numeric as num
def pasivate(pos,sub=[],nneig=3):
   """ Return the position of the H atoms to pasivate the edges. """
   #TODO include consideration of lattice vectors for ribbons
   ## List all the atoms of a given kind with less than nneig neighbours
   nn = num.count_neig(pos,pos,1.5)
   rows,cols = num.dists(pos,pos,nn,1.5)
   rows -= 1   # because python starts counting at 0
   cols -= 1   #
   needH,aux_sub = [],[]
   for i in range(len(pos)):
      if len(cols[rows==i]) < nneig:
         needH.append( (i,cols[rows==i]) )
         aux_sub.append( sub[i] )
   new_atoms, new_sub = [],[]
   for i in range(len(needH)):
      at,neig = needH[i]
      at2neig = pos[at]
      v1 = pos[neig[0]]
      v2 = pos[neig[1]]
      r1 = v1-at2neig
      r2 = v2-at2neig
      r_orto = np.cross(r1,r2) # this vector determines the plane
      # angle for the new atom
      ang = np.arccos(np.dot(r1,r2)/(np.linalg.norm(r1)*np.linalg.norm(r2)))
      listvecs = [r1,r2]
      #  try to put new atom in any of the possible positions until finding
      # the missing atom
      ivec = 0
      r3 = geo.rotation(listvecs[ivec],ang,r_orto,deg = False)
      while geo.vecinlist(r3,listvecs):
         r3 = geo.rotation(listvecs[ivec],ang,r_orto,deg = False)
         ivec += 1
      # position of the new atom
      v3 = r3 + pos[at]
      new_atoms.append(v3)
      new_sub.append(-1*aux_sub[i])
      #try: new_sub.append(-1*sub[i])
      #except IndexError: pass
   return new_atoms,new_sub


if __name__ == '__main__':
   n = 35
   l = 2
   passivate = False
   for n in range(10,65,5):
      for l in [1,2]:
         ats,pos,latt,subs = armchair(n)
         #ats,pos,latt,subs = simple(n)
         #ats,pos,latt,subs = kagome(n)
         A = UnitCell(ats,pos,latt,subs)
         if passivate: A.pasivate()
         r = 'cells/ac'
         if l>1: A.multilayer(l)
         if passivate: A.to_xyz(r+'_n%s_l%s_H.xyz'%(n,l))
         else: A.to_xyz(r+'_n%s_l%s.xyz'%(n,l))

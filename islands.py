#!/usr/bin/python3
# -*- coding: UTF-8 -*-

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
      self.sub = np.array(subs,int)
      if self.check() and len(self.ats) > 0:
         print('here')
         self.get_geo_info()
      else: pass
      if pasivate: self.pasivate()
   def get_geo_info(self):
      """ Calculate geometric attributes """
      if not self.check(): print('WARNING')   #TODO Log
      self.center = np.mean(self.pos,axis=0)
      lmin = np.min(self.pos,axis=0)
      lmax = np.max(self.pos,axis=0)
      self.lims = [lmin,lmax]
   def pasivate(self):
      latt = []
      hs,subh = pasivate(self.pos,sub=self.sub)
      self.ats = np.append(self.ats,['H' for _ in self.pos])
      self.pos = np.append(self.pos,hs)
      self.sub = np.append(self.sub,subh)
   def multilayer(self,lN):
      self.ats,self.pos,self.sub = multilayer(self.pos,self.sub,N=lN)
      #self.center = np.mean(np.array(self.pos),axis=0)
      self.get_geo_info()
   def check(self):
      if len(self.ats) == len(self.pos) == len(self.sub): return True
      else: return False
   def to_xyz(self,fname=None):
      print(len(self.ats))
      print('')
      for a,r,s in zip(self.ats,self.pos,self.sub):
         print('%s   %s   %s   %s   %s'%(a,r[0],r[1],r[2],s))
      if fname != None:
         f = open(fname,'w')
         f.write(str(len(self.ats))+'\n\n')
         for a,r,s in zip(self.ats,self.pos,self.sub):
            f.write('%s   %s   %s   %s   %s\n'%(a,r[0],r[1],r[2],s))
         f.close()
   def from_xyz(self,fname,pasivate=False):
      from IO import xyz
      self.ats,self.pos,self.latt,self.sub = xyz(fname)
      #self.center = np.mean(np.array(self.pos),axis=0)
      self.get_geo_info()
   def plot(self,fname=None):
      plot_cell(self.pos,self.latt,fname=fname)
   def plot_sublattice(self,fname=None):
      X = self.pos[:,0]
      Y = self.pos[:,1]
      S = self.sub
      fig, ax = plt.subplots()
      ax.scatter(X,Y,c=S,s=100,edgecolors='none')
      plt.show()
   def __str__(self):
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
   ats = ['C' for _ in pos]
   return ats,pos,latt,subs

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
   if show: plot_cell(pos,latt,tit='Armchair Cell (%s)'%(N))
   ats = ['C' for _ in pos]
   return ats,pos,latt,sub


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
   sub = np.array(sub)
   ### Plot
   if show: plot_cell(pos,latt,tit='ZigZag Cell (%s)'%(N))
   ats = ['C' for _ in pos]
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
   pos = np.array(pos)
   sub = np.array(sub)
   ## Re-Center the unit cell
   if cent:
      C = np.mean(pos,axis=0)
      for i in range(len(pos)):
         pos[i] -= C
   if show: plot_cell(pos,latt,tit='Simple Cell %sx%s'%(N,N))
   ats = ['C' for _ in pos]
   return ats,pos,latt,sub

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
   sub = np.array(sub)
   if show: plot_cell(pos,tit='Triangular Island (%s)'%(N))
   ats = ['C' for _ in pos]
   return ats,pos,[],sub

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
   ats = ['C' for _ in pos]
   #hs,subh = pasivate(pos,sub=sub)
   #pos += hs
   #ats += ['H' for _ in hs]
   #sub += subh
   return ats,pos,[],sub




def vec_in_list(v,l,eps=0.000000001):
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
         try: new_sub.append(sub[j])
         except: pass
   return new_ats, new_pos, new_sub

def pos2xyz(pos,latt,at='C',sub=[],fname='lattice.xyz'):
   """ at has to be a string or a list/array of strings """
   LG = logging.getLogger('IO.pos2xyz')
   if isinstance(at,str):
      LG.info('Only one atom provided. Using %s for all the atoms'%(at))
      at = [at for _ in pos]
   with open(fname,'w') as f:
      ## Write number of atoms
      f.write(str(len(pos))+'\n')
      LG.debug('Written the number of atoms')
      ## Write lattice vectors
      for v in latt:
         f.write('[%s,%s,%s]'%(v[0],v[1],v[2]))
      f.write('\n')
      LG.debug('Written the lattice vectors')
      ## Write atoms positions
      for i in range(len(pos)):
         a,r = at[i],pos[i]
         try:
            s = sub[i]
            f.write(a+'   %s   %s   %s   %s\n'%(r[0],r[1],r[2],s))
         except: f.write(a+'   %s   %s   %s\n'%(r[0],r[1],r[2]))
      LG.debug('Written all the atomic positions')


import numeric as num
def pasivate(pos,sub=[],nneig=3):
   """ Return the position of the H atoms to pasivate the edges. """
   #TODO include consideration of lattice vectors for ribbons
   ## List all the atoms of a given kind with less than nneig neighbours
   nn = num.count_neig(pos,pos)
   rows,cols = num.dists(pos,pos,nn)
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
   A = UnitCell([],[],[],[])
   A.from_xyz('cells/ac_n1_l1.xyz') #'test.xyz')
   A.get_geo_info()
   print(A)
   A.plot(fname='test.png')
   A.plot_sublattice()

   exit()
   ats,pos,latt,subs = armchair(2)
   A = UnitCell(ats,pos,latt,subs)
   A.pasivate()
   A.to_xyz('test.xyz')

   #print(len(ats)+len(hs))
   #print('')
   #for a,p,s in zip(ats,pos,subs):
   #   print(a,p[0],p[1],p[2],s)
   #for p,s in zip(hs,subh):
   #   print('H',p[0],p[1],p[2],s)

   exit()
   ## Read island type, size and layers from standard input  (TODO argparse)
   # Usage: python islands.py armchair 20 2  ---> armchair island 20x20 bilayer
   import sys
   pas = False
   try: func = sys.argv[1]
   except:
      print('Type of unit cell not specified.')
      print('Available unit cells:')
      print('   - simple      - armchair')
      print('   - zigzag      - zigzag_triangle')
      print('   - Mullen')
      print('\nExample of usage:')
      print('python islands.py armchair 3 1')
      exit()
   try: N = int(sys.argv[2])
   except IndexError:
      print('WARNING: Supercell index not specified. Assumed 0')
      N = 0
   try: lN = int(sys.argv[3])
   except IndexError:
      print('WARNING: Multilayer index not specified. Assumed 1')
      lN = 1

   ##Setup the proper function
   funcs = [armchair, zigzag, zigzag_triangle,simple,mullen]
   funcs_names = [f.__name__ for f in funcs]
   funcs = dict(zip(funcs_names,funcs))
   acronym = {'armchair':'ac','zigzag':'zz','zigzag_triangle':'zzt', # islands
              'simple':'simple', 'mullen':'mullen'}
   try: func = funcs[func]
   except KeyError:
      print('Requested geometry (%s) not implemented'%(func))
      exit()

   ##Do the calculation
   print('Using funcion',func.__name__,'with index',N)
   ats,pos,latt,sub = func(N)

   ## Pasivation
   if pas:
      hs,subh = pasivate(pos,sub=sub)
      pos += hs
      ats += ['H' for _ in hs]
      sub += subh

   ## Multilayer
   if lN > 1:
      ats,pos,sub = multilayer(pos,sub,N=lN)

   if pas: nam = acronym[func.__name__]+'_n%s_l%s_H.xyz'%(N,lN)
   else: nam = acronym[func.__name__]+'_n%s_l%s.xyz'%(N,lN)

   C = np.mean(pos,axis=0)
   for i in range(len(pos)):
      pos[i] -= C
   pos2xyz(pos,latt,at=ats,sub=sub,fname=nam)

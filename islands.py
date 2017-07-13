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
   brick = [a*np.array([1.,0,b]),
            a*np.array([1/2.0,ap,-b]),
            a*np.array([-1/2.,ap,b]),
            a*np.array([-1,0,-b])]
   sublatt = ['A','B','A','B']  # XXX check order
   
   vectors = [a*np.array([3.,0.,0.]),
              a*np.array([0,2.*ap,0])]  # to expand the unit cell
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
   brick = [a*np.array([1.,0,b]),
            a*np.array([1/2.0,ap,-b]),
            a*np.array([-1/2.,ap,b]),
            a*np.array([-1,0,-b]),
            a*np.array([-1/2.,-ap,b]),
            a*np.array([1/2.,-ap,-b])]
   sublatt = ['A','B','A','B','A','B']  # XXX check order
   
   vectors = [a*np.array([3/2.,3*ap,0]),
              a*np.array([3,0,0])]
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
   brick = [a*np.array([1.,0,b]),
            a*np.array([1/2.0,ap,-b]),
            a*np.array([-1/2.,ap,b]),
            a*np.array([-1,0,-b]),
            a*np.array([-1/2.,-ap,b]),
            a*np.array([1/2.,-ap,-b])]
   
   vectors = [a*np.array([3/2.,-ap,0]),
              a*np.array([3/2.,ap,0])]
   latt = [(N+1)*(vectors[0]+vectors[1]),
           (N+1)*(-vectors[0]+2*vectors[1])]

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
      for r in brick:
         w = r+v
         # Avoid repited atoms
         if not vec_in_list(w,pos): pos.append(w)
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
   sublatt = ['A','B']

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
   brick = [a*np.array([1.,0,b]),
            a*np.array([1/2.0,ap,-b]),
            a*np.array([-1/2.,ap,b]),
            a*np.array([-1,0,-b]),
            a*np.array([-1/2.,-ap,b]),
            a*np.array([1/2.,-ap,-b])]
   vectors = [a*np.array([3/2.,-ap,0]),
              a*np.array([3/2.,ap,0])]
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
   pos = []
   for v in all_vecs:
      for r in brick:
         w = r+v
         # Avoid repited atoms
         if not vec_in_list(w,pos): pos.append(w)
   if show: plot_cell(pos,tit='Triangular Island (%s)'%(N))
   ats = ['C' for _ in pos]
   return ats,pos,[]

def mullen(Nx,Ny=4,pas=False):
   ats,pos,_,sub = ribbon_armc(Nx,Ny)
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
   hs,subh = pasivate(pos,sub=sub)
   pos += hs
   ats += ['H' for _ in hs]
   sub += subh
   return ats,pos,[],sub




def vec_in_list(v,l,eps=0.000000001):
   """ Returns True if vector v is in the list of vectors l (also in util.py)"""
   for x in l:
      if np.linalg.norm(x-v) < eps: return True
   return False

def plot_cell(pos,latt=[],tit=None):
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
      if len(latt) == 1: latt.append(np.array([0,0,0])) # XXX Shame on you!!!
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
   plt.show()

def multilayer(pos,sub=[],N=2):
   """ Generates the positions for a multilayer ABC... """
   new_ats, new_pos, new_sub = [], [], []
   rs = [i*np.array((1.4,0,1.4)) for i in range(N)]
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
   needH = []
   nn = num.count_neig(pos,pos)
   rows,cols = num.dists(pos,pos,nn)
   rows -= 1
   cols -= 1
   aux_sub = []
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
      try: new_sub.append(sub[i])
      except IndexError: pass
   return new_atoms,new_sub




if __name__ == '__main__':
   ## Read island type, size and layers from standard input  (TODO argparse)
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
   pos2xyz(pos,latt,at=ats,sub=sub,fname=nam)

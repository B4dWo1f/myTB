#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import os
import numpy as np
import numeric as num
from random import uniform,choice
from random import random as rand
from scipy.sparse import coo_matrix
from scipy.spatial import KDTree
from itertools import product
import logging
import log_help
LG = logging.getLogger(__name__) # Logger for this module


def analyze(atoms,points,pairs=[],maxneigh=5,fname='cell.info'):
   """
     Find nearest neighbor distances for each kind of hopping.
     pairs: List of hoppings to analyze. If empty all distances are studied.
     maxneigh: Max number of neighbours to be considered
   """
   diff_names = set(atoms)
   if len(pairs) == 0:
      LG.info('No pairs provided, so all combinations are studied')
      pairs = []
      for i in product(diff_names, repeat=2):
         string = '-'.join(sorted(i))
         if string not in pairs: pairs.append(string)
   keys,values = [],[]
   for bond in pairs:
      at1,at2 = bond.split('-')
      N = 5   # number of neighbors to check
      means = []
      for i in range(len(atoms)):
         if atoms[i] != at1: continue
         dists = []
         for j in range(len(atoms)):
            if i == j: continue
            r = points[i]-points[j]
            dists.append((np.linalg.norm(r),j))
         dists = sorted(dists,key=lambda x:x[0])
         dists = dists[0:N]
         aux = []
         for i_dist in range(len(dists)):
            x = dists[i_dist]
            if atoms[x[1]] == at2: aux.append(x[0])
            #if len(aux) >= N: break
         dists = aux

         relevant_mean = []
         for i in range(1,len(dists)):
            aux = dists[0:i]
            m,s = np.nanmean(aux),np.nanstd(aux)
            if s > 0.1: break
            relevant_mean.append(m)
         if len(relevant_mean) > 0: means.append(np.nanmean(relevant_mean))
      LG.info('%s %s atoms have some %s neighbor'%(len(means),at1,at2))
      M,S = np.nanmean(means),np.nanstd(means)
      keys.append(bond)
      values.append((M,S))
      LG.info('Bond %s, has 1st neig distance: %s'%(bond,M))
      if S < 0.1:
         LG.warning('STD of the %s distance is suspiciously low'%(bond))
   return dict(zip(keys,values))


def snap(P,points,retpoint=False):
   """
     Return the value in points closest to P
     retpoint: if true return the index and the point itself
               if false return only the index
   """
   idx = np.nanargmin(((points - P)**2).sum(axis = -1))
   if retpoint: return idx,points[idx]
   else: return idx


def circle(x,d):
   """
     Returns the y coordinate of a point with x coordinate and diameter d
   """
   return np.sqrt((d/2)**2 - x*x)


def regular_polygon(n,l,s=0.,z=0.):
   """
   Returns the vertices of a regular poligon with n sides of length l.
   First vertex is placed at r*(cos(s),sin(s))
   """
   if n>1: r = l/(2*np.sin(np.pi/n))
   else: r = 0.
   points = []
   for i in range(n):
      theta = np.radians(s) + i*2*np.pi/n
      points.append( (r*np.cos(theta),r*np.sin(theta),z) )  # XXX Z component?
   return points


@log_help.log2screen(LG)
#def defects(pos,d=None,alpha=0.,hollow=True,retpoint=False):
def defects(N,pos,sub,lay,bonds=None,d=None,alpha=0.,hollow=True,retpoint=False):
   """
     Returns pais of atoms at an approximate distance of d from the list pos
     retpoint: if true return the index and the point itself
               if false return only the index
     bonds = Matrix of bonds for intra cell
     ** Assumes the points are centered in 0
     * also assumes unit cell like:
                   -----
                  /     \ 
                  \     /
                   -----
   """
   center = np.mean(pos,axis=0)
   def find_neig(M,ind):
      """ returns the index of the neighbouring atoms of atom ind """
      if not hasattr(ind, '__iter__'): ind = [ind]
      ret = []
      for i in ind:
         r,c = M.row,M.col
         ret.append( c[r==i] )
      return ret
   ## Select sublattice   ---   Hollow Vs Connected
   l = np.max(lay)
   LG.info('Vacancies introduced in layer: %s'%(l))
   auxX = np.array(range(len(lay)))  # original index in pos
   ## Index of atoms sublattice A
   sub_atsA = np.where((lay==l)&(sub==1))[0]
   na = sub_atsA.shape[0]
   ## Index of atoms sublattice B
   sub_atsB = np.where((lay==l)&(sub==-1))[0]
   nb = sub_atsB.shape[0]
   
   lena = np.mean([len(find_neig(bonds,sub_atsA[i])) for i in range(na)])
   lenb = np.mean([len(find_neig(bonds,sub_atsB[i])) for i in range(nb)])

   if hollow: ###### indices of hollow atoms
      if lena < lenb: sub_ats = sub_atsA
      else:           sub_ats = sub_atsB
   else: ########### indices of connected atoms
      if lena > lenb: sub_ats = sub_atsA
      else:           sub_ats = sub_atsB
   ideal = regular_polygon(N,d,s=alpha,z=l)
   ideal = [r-center for r in ideal]   # TODO check
   #indices = [snap(p, pos ) for p in ideal]
   # index referred to the subset A/B
   indices = [snap(p, pos[sub_ats] ) for p in ideal]
   # indices referred to the full list of atomic positions
   indices = [auxX[sub_ats][x] for x in indices]
   return indices
   #### Select atoms for defects
   ##if N == 1:   ## 1 defect
   ##   C = np.mean(self.pos[sub_ats],axis=0)
   ##   ind = geo.snap(C,self.pos[sub_ats])
   ##   indices = [sub_ats[ind]]
   ##elif N == 2: ## 2 defects
   ##   #if d == None: d = np.sqrt(3)*(np.max(self.x)-np.min(self.x))/12
   ##   if d == None: d = latt[0]/np.sqrt(3)
   ##   regular_polygon(N,d) #####################################################
   ##   msg  = 'Including two vacancies %s Ansg apart'%(d)
   ##   msg += ' with an angle %s'%(alpha)
   ##   LG.info(msg)
   ##   # returns indices referred to subset
   ##   inds = geo.defects(self.pos[sub_ats],d=d,alpha=alpha)
   ##   indices = [sub_ats[i] for i in inds]
   ##   msg = 'Vacancies placed at'
   ##   for i in indices:
   ##      msg += ' %s'%(i)
   ##   LG.info(msg)
   ##   rd = self.pos[indices[0]] - self.pos[indices[1]]
   ##   ra = np.degrees(np.arctan2(rd[1],rd[0]))
   ##   rd = np.linalg.norm(rd)
   ##   LG.warning('Requested-dist/Real-dist: %s/%s'%(d,rd))
   ##   LG.warning('Requested-angle/Real-angle: %s/%s'%(alpha,ra))
   ##elif N == 3: ## 3 defects
   ##   LG.critical('Implemention not finished')
   ##   ak = np.linalg.norm(self.latt[0])/2  # size of kagome lattice
   ##   r3 = np.sqrt(3)
   ##   ideal = [np.array([ r3*ak/4,  0 ,0]),
   ##            np.array([-r3*ak/4,ak/2,0]),
   ##            np.array([-r3*ak/4,-ak/2,0])]
   ##   indices = [geo.snap(p, self.pos[sub_ats] ) for p in ideal]
   ##elif N > 3:  #XXX Check!!!!
   ##   LG.warning('Experimental implementation of N>3 defects')
   ##   indices = []
   ##   while len(indices) < N:
   ##      indices = indices + list(set(choice(self.INDS)))

   #X = pos[:,0]
   #Y = pos[:,1]
   #Z = pos[:,2]
   #mx = np.max(X)
   #r = d/2
   #LG.debug('Max distance allowed: %s'%(2*mx))
   #LG.debug('Max recommended distance: %s'%(np.sqrt(3) * mx/2))
   #if r > mx:
   #   LG.critical('Distance between defects bigger than island.')
   #   print('max dist:',2*mx)
   #   print('max recommended dist:',np.sqrt(3) * mx/2)
   #   exit()
   #elif r > np.sqrt(3) * mx/2:
   #   LG.warning('Distance to borders bigger than distances between vacancies')
   #x = np.cos(np.radians(alpha)) * r
   #z0 = choice(Z)          # considering C3 symmetry
   ### Ideal points
   #P0 = np.array( (x, circle(x,d), z0) )
   #P1 = np.array( (-P0[0], -P0[1], z0) )
   ### Real points
   #i0 = snap(P0,pos)
   #i1 = snap(P1,pos)
   #if retpoint: return [i0,i1],[pos[i0],pos[i1]]
   #else: return [i0,i1]


#@log_help.log2screen(LG)
#def layer(pos,dist=1.5,lim=100,eps=0.2):
#   Zs = pos[:,2]
#   Z = np.sort(Zs)
#   dZ = np.diff(Z)
#   mean_dz,std_dz = np.mean(dZ),np.std(dZ)
#   bi = len(dZ[dZ>mean_dz+2*std_dz]) +1
#   hist, bin_edges = np.histogram(Zs,bins=bi-1)
#   L = [x-bi//2 for x in range(bi+1)]
#   lays = [ L[np.argmin(np.abs(bin_edges-z))] for z in Zs]
#   print(lays)
#   exit()
#   return lays

@log_help.log2screen(LG)
def layer(pos,dist=1.5,lim=100,eps=0.2):
   """
    Determine the layer distribution of the atoms (so far only for Z direction)
   """
   LG.info('Guessing the layers')
   zs = list(set(pos[:,2]))
   lay = [i for i in range(len(zs))]
   aux = []
   for i in range(len(zs)):
      z = zs[i]
      l = lay[i]
      for r in pos:
         if abs(r[2]-z) < eps: aux.append(l)
   if len(pos) != len(aux):
      LG.critical('Incorrect assignation of layers')
      exit()
   else:
      LG.info('Layers done')
      return aux


#def fsublattice(intra,dist=1.5,lim=100):
#   print('Sublattice from fortran')
#   aux = np.array((intra.row,intra.col)).transpose()
#   subs = num.sublattice(aux)
#   exit()
@log_help.log2screen(LG)
def check_sublattice(intra,subs,dist=1.5,lim=100):
   """
    This function checks that the sublattice is well calculated (every A atom
    has ONLY B neighbors)
    intra has to be a coo_matrix containing the neighbors of the unit cell
   """
   r,c = intra.row, intra.col
   for i in r:
      aux = list(set([subs[j] for j in c[r==i]]))
      if len(aux) > 1 or aux[0] == subs[i]: return False
   return True

@log_help.log2screen(LG)
def sublattice(intra,dist=1.5,lim=100):
   """
    TODO: Fortran this
    intra has to be a coo_matrix containing the neighbors of the unit cell
   """
   subs = [10 for _ in range(intra.shape[0])]
   subs[0] = 1
   r,c = intra.row, intra.col
   for n in range(10000):
      for i in range(intra.shape[0]):
         #print('Atom',i,'with sub: ',subs[i],'has neigs:',c[r==i])
         for j in c[r==i]:   # neighbors of atom i
            if abs(-1*subs[i]) < 2: subs[j] = -1*subs[i]
            else: 
               #print('   atom',i,'still has no sublattice')
               #print('        but its neighbor',j,'is sublatt:',subs[j])
               #print('        so atom',i,'should be',-1*subs[j])
               if abs(-1*subs[j]) < 2: subs[i] = -1*subs[j]
               else: pass
            #try: subs[j] = -1*subs[i]
            #except KeyError:
            #   print('   atom',i,'still has no sublattice')
            #   print('        but its neighbor',j,'is sublatt:',subs[j])
            #   print('        so atom',i,'should be',-1*sub_dict[subs[j]])
            #   try: subs[i] = -1*subs[j]
            #   except KeyError: pass
      #types = set([type(x) for x in subs])
      N_missing = np.where(np.array(subs)>2)[0].shape[0]
      if N_missing == 0: #len(list(types)) == 1:
         LG.info('Sublattice, %s iterations'%(n))
         return subs
      if n > 100:
         if n%1000 == 0:
            LG.warning('Probably there\'s an error in Sublattice')
            exit()
#def sublattice_old(pos,dist=1.5,lim=100):
#   """
#     This function is in beta testing. I think it should work for any
#     acceptable crystal
#   """
#   LG.info('Guessing the sublattice')
#   class dummy(object):
#      def __init__(self,pos,sub=0):
#         self.pos = pos
#         self.sub = sub
#   new_pos = [dummy(p) for p in pos]
#   tree = KDTree(pos)
#   sub_dict = {'A':1,'B':-1}
#   tcid_bus = {1:'A',-1:'B'}
#   new_pos[0].sub = 'A' ## Random start
#   areall,cont = False, 0
#   while not areall and cont < lim:
#      LG.debug('subalttice, run %s'%(cont))
#      for i in range(len(new_pos)):
#         p = new_pos[i]
#         if isinstance(p.sub,str):  # I know its sublattice
#            d,ind = tree.query(p.pos, k=5)
#            d,ind = d[1:],ind[1:]
#            mask = np.where(d<=dist,True,False)
#            d = d[mask]
#            ind = ind[mask]
#            for j in ind:
#               if not isinstance(new_pos[j].sub,str):
#                  new_pos[j].sub = tcid_bus[-1*sub_dict[p.sub]]
#      are_all = [isinstance(p.sub,str) for p in new_pos]
#      areall = all(are_all)
#      cont += 1
#   LG.debug('Sublattice assignation in %s iterations.'%(cont))
#   if cont == lim: LG.warning('The sublattice assignation may not have worked')
#   subs = []
#   for p in new_pos:
#      subs.append(p.sub)
#   if len(pos) == len(subs): LG.debug('Sublattice assignation done')
#   else: LG.warning('The sublattice assignation didn\'t work')
#   LG.info('Sublattice done')
#   return subs

def reciprocal(latt_vec):
   """
     Should be valid for 0D, 1D and 2D. CAREFUL!!! NOT TESTED!!!
     Uses formula 5.3 in the Ashcroft Book
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


def vecfromcoef(coef,vecs):
   """
     Return the real vector out of the coeficients(cn) of the vectors of a
     basis(an):   v = c1*a1 + c2*a2 + ... + cn*an
   """
   if len(coef) != len(vecs): LG.error('Unable to form vector')
   aux = np.array([0.,0.,0.])
   for i in range(len(vecs)):
      aux += coef[i] * vecs[i]
   return aux


def vecinlist(vec,lista,eps=0.00001):
   """ Checks if a vector vec is in a given list """
   for v in lista:
      eq = vec - v
      if np.linalg.norm(eq) < eps: return True
      else: pass
   return False


#def tree_neig(pos,latt,dist=1.5,nvec=5):
#   """
#     Find nearest neighbours at distance < dist. Returns 2 lists.
#     cosa: contains tuples of the form (cell_index, i_atom, j_atom) where
#          the i_atom and j_atom are the index of the atoms in the unit cell
#     ** Not efficient in parallel
#   """
#   LG.debug('%s atoms'%(len(pos)))
#   aux = [p for p in product([0,1,-1], repeat=len(latt)) ]
#   aux = sorted(aux,key=np.linalg.norm)  #XXX Check correct order
#   perms = []
#   VIL = vecinlist   # Local rename
#   for i in range(len(aux)):
#      v = np.array(aux[i])
#      if not VIL(v,perms) and not VIL(-v,perms): perms.append(aux[i])
#   perms = sorted(perms,key=np.linalg.norm)
#   # all_vecs contains the vectors: 0, a1, a2, a1+a2, a1-a2
#   all_vecs = []
#   for p in perms:
#      all_vecs.append( vecfromcoef(p,latt) )
#   # indices: List containing the cell index for each of the atoms stored in
#   #                                                           the list all_pos
#   # all_pos: List with the position of all the atoms in all the
#   #                                 cells (cetered at the vectors in all_vecs)
#   all_pos,indices = [],[]
#   for i in range(len(all_vecs)):
#      v = all_vecs[i]
#      for r in pos:
#         all_pos.append( v+r )
#         indices.append(i)
#
#   if len(all_pos) != len(indices):   #XXX Raise an exception/error?
#      LG.critical('Wrong assignment of indices and cells')
#
#   tree = KDTree(all_pos)
#   cosa = []
#   for ip in range(len(pos)):
#      LG.debug('Atom %s in the unit cell is connected to:'%(ip))
#      r = pos[ip]
#      d,ind = tree.query(r, k=nvec)
#      d,ind = d[1:],ind[1:]
#      mask = np.where(d<=dist,True,False)
#      # d are the distances between points
#      # ind are the positions in the all_pos list
#      for i in ind[mask]:
#         r1 = all_pos[i]-all_vecs[indices[i]]
#         ind = snap(r1,pos)
#         cosa.append((indices[i],ip,ind))
#         LG.debug('  atom %s in cell %s'%(ind,indices[i]))
#   del tree,indices,mask,all_pos
#   return cosa, all_vecs
#
#def brute_force(pos,latt,dist=1.5,nvec=5,ncpus=4):
#   """ Brute force search for neighbours """
#   aux = [p for p in product([0,1,-1], repeat=len(latt)) ]
#   aux = sorted(aux,key=np.linalg.norm)  #XXX Check correct order
#   perms = []
#   VIL = vecinlist   # Local rename
#   for i in range(len(aux)):
#      v = np.array(aux[i])
#      if not VIL(v,perms) and not VIL(-v,perms): perms.append(aux[i])
#   perms = sorted(perms,key=np.linalg.norm)
#   # all_vecs contains the vectors: 0, a1, a2, a1+a2, a1-a2
#   all_vecs = [ vecfromcoef(p,latt) for p in perms]
#   names = ['intra','x','y','xy','xmy']
#   cosa = []
#   for ir1 in range(len(pos)):
#      LG.debug('Atom %s in the unit cell is connected to:'%(ir1))
#      r1 = pos[ir1]
#      for i in range(len(all_vecs)):
#         v = all_vecs[i]
#         for ir2 in range(len(pos)):
#            r2 = pos[ir2] + v
#            if 0.0 < np.linalg.norm(r1-r2) < dist:
#               LG.debug('  atom %s in cell %s'%(ir2,i))
#               ## cell index;  atom index;  atom index
#               cosa.append( (i,ir1,ir2) )
#   return cosa,all_vecs,names


@log_help.log2screen(LG)
def fneig(pos,latt,fol='./',dist=1.5,nvec=5,ncpus=4,force=False):
   """
    Fortran implementation of the neighbor finding algorithm
   """
   #aux = [p for p in product([0,1,-1], repeat=len(latt)) ]
   #aux = sorted(aux,key=np.linalg.norm)  #XXX Check correct order
   #perms = []
   #VIL = vecinlist   # Local rename
   #for i in range(len(aux)):
   #   v = np.array(aux[i])
   #   if not VIL(v,perms) and not VIL(-v,perms): perms.append(aux[i])
   #perms = sorted(perms,key=lambda x: np.linalg.norm(x))
   # TODO this should be automatic for 3D  extension
   if len(latt) == 2:
      perms = [(0,0),(1,0),(0,1),(1,1),(1,-1)]
      names = ['intra','x','y','xy','xmy']
   elif len(latt) == 1:
      perms = [(0,0),(1,0)]
      names = ['intra','x']
   elif len(latt) == 0:
      perms = [(0,0)]
      names = ['intra']
   else: LG.critical('Dimensionality not implemented')
   # all_vecs contains the vectors: 0, a1, a2, a1+a2, a1-a2
   all_vecs = [ vecfromcoef(p,latt) for p in perms]
   neigs = []
   for i in range(len(all_vecs)):
      LG.info('Looking for neighbors in cell %s'%(names[i]))
      try:
         if force: raise
         if os.stat(fol+names[i]+'.bond').st_size==0:  # empty file
            LG.info('No neighbouts for %s'%(names[i]))
            rows = []
            cols = []
         else:
            M = np.matrix(np.loadtxt(fol+names[i]+'.bond',dtype=int))
            r = np.array(M[:,0])
            c = np.array(M[:,1])
            rows = np.reshape(r,(max((r.shape)),))
            cols = np.reshape(c,(max((c.shape)),))
         LG.info('Read from file: %s'%(fol+names[i]+'.bond'))
      except:
         LG.info('Failed. Calculating with fortran')
         ## Necessary because of f2py problem with allocatable
         nn = num.count_neig(pos,pos+all_vecs[i],dist)
         if nn == 0:
            LG.info('No neighbours in cell %s'%(names[i]))
            if fol != '': np.savetxt(fol+names[i]+'.bond',[],fmt='%d')
            continue
         rows,cols = num.dists(pos,pos+all_vecs[i],nn,dist)
         rows -= 1   # because python counts from 0
         cols -= 1   # because python counts from 0
         aux = np.column_stack((rows,cols))
         if fol != '': np.savetxt(fol+names[i]+'.bond',aux,fmt='%d')
         #if fol != '': np.savetxt(fol+names[i]+'.H',(rows,cols),fmt='%d')
      neigs.append( ((rows,cols),all_vecs[i],names[i]) )
   return neigs

def get_points(recip,N=3):
   """ Get special points in the reciprocal space """
   lista = np.linspace(0,1,N+1)
   perms = [p for p in product(lista, repeat=len(recip)) ]
   points = []
   ii = 0
   for p in perms:
      aux = np.array([0.,0.,0.])
      for i in range(len(p)):
         aux += p[i]*recip[i]
      points.append(aux)
      ii+=1
   return points

def recorrido(points,nk):
   """
     Returns a list of points (np.array) with nk points between each pair of
     points in the provided list of points.
     This function is dimension independent.
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
      for idim in range(len(P1)):
         coors.append(np.linspace(P1[idim],P2[idim],N+ret[0],endpoint=ret[1]))
      for p in zip(*coors):
         RECORRIDO.append(np.array(p))
   return RECORRIDO


from numpy import cos,sin
def rotation(v,Q,u=np.array([0,0,1]),deg=True):
   """
     Rotates a vector v by an angle Q around the axis u
     Q in degrees by default
   """
   u = u/np.linalg.norm(u)
   ux = u[0]
   uy = u[1]
   uz = u[2]
   if deg : q = np.radians(Q)
   else: q = Q
   cos1 = 1.-cos(q)
   R = np.array(
   [[cos(q)+(ux**2)*cos1 , ux*uy*cos1-uz*sin(q), ux*uz*cos1+uy*sin(q)],
    [uy*ux*cos1+uz*sin(q), cos(q)+(uy**2)*cos1 , uy*uz*cos1-ux*sin(q)],
    [uz*ux*cos1-uy*sin(q), uz*uy*cos1+ux*sin(q), cos(q)+(uz**2)*cos1]])
   vec = np.dot(R, v)
   return vec


if __name__ == '__main__':
   import IO
   import sys
   try: fname = sys.argv[1]
   except IndexError:
      print('File not specified')
      exit()
   
   atoms,points,_,_ = IO.read.xyz(fname)
   A = analyze(atoms,points,pairs=[],maxneigh=5,fname='cell.info')
   print('Atomic distances:')
   for k,v in A.items():
      print('  *%s: %.4f +-%.4f'%(k,v[0],v[1]))

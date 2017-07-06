#!/usr/bin/python3
# -*- coding: UTF-8 -*-


import os
import IO
import numpy as np
import geometry as geo
from scipy.sparse import coo_matrix
from copy import deepcopy
from os import listdir
from os.path import isfile, join
## LOG
import logging
import log_help
LG = logging.getLogger(__name__)



class Base_Element(object):
   def __init__(self,n=0,elem='',onsites={},pos=np.array(())):
      self.place = n
      self.element = elem
      # s, px, py, pz, dxy, dyz, dzx, dx2y2, d3z2r2
      order = {'s':0, 'px':1, 'py':2, 'pz':3,
               'dxy':4, 'dyz':5, 'dzx':6, 'dx2y2':7, 'd3z2r2':8}
      self.orbitals = sorted(list(onsites.keys()),key=lambda x:order[x])
      self.onsite = onsites
      self.position = pos
   def update(self,entries): self.__dict__.update(entries)
   def copy(self): return deepcopy(self)
   def __str__(self,w=80):
      msg = '-'*w
      msg += '\n    Atom: %s,  type: %s\n'%(self.place,self.element)
      msg += '   on-site:'
      for o,oe in self.onsite.items():
         msg += '  %s: %s,'%(o,oe)
      msg = msg[0:-1]+'\n'  # To delete last comma
      msg += '   postion: ('
      for r in self.position:
         msg += str(r)+','
      msg = msg[0:-1]+')\n'  # To delete last comma
      try: msg += 'sublattice: %s\n'%(self.sublattice)
      except AttributeError: pass
      try: msg += '     layer: %s'%(self.layer)
      except AttributeError: pass
      msg += '\n'
      return msg


class Base(object):
   def __init__(self, elements=[],latt=[],cent=True):
      self.elements = elements  # List of Base_Elements
      if cent: self.center()    # center unit cell in (0,0,0)
      ## Positions of the atoms
      self.pos = np.array([ e.position for e in self.elements ])
      self.x = self.pos[:,0]
      self.y = self.pos[:,1]
      self.z = self.pos[:,2]
      ## Lattice vectors
      self.latt = latt          # Lattice vectors
      self.recip = geo.reciprocal(latt)
      ## list of index per atom
      self.get_indices()
      self.update_basis()
      ## auxiliary list of indices
      self.inds = np.array([i for i in range(len(self.elements))])
   def copy(self): return deepcopy(self)
   def update_basis(self):
      self.basis = []
      for E in self.elements:
         for o in E.orbitals:
            try: ID = (E.place,E.element,o,E.sublattice)
            except: ID = (E.place,E.element,o)
            self.basis.append(ID)
   def center(self): geo.center_cell(self)
   def __getitem__(self,index):  return self.elements[index]
   def __iter__(self): return (x for x in self.elements)
   def __str__(self):
      msg = 'Lattice vectors: '
      for v in self.latt:
         msg += '(%5.3f, %5.3f, %5.3f)\n                 '%(v[0],v[1],v[2])
      msg += '\n%s Atoms in the unit cell:\n'%(len(self.elements))
      for E in self.elements:
         r = E.position
         msg += ' %s   (%5.3f, %5.3f, %5.3f)\n'%(E.element,r[0],r[1],r[2])
      return msg
   def get_indices(self):
      """
        This function adds a "indices" attribute in each base element
        containing the matrix indices corresponding to that element
      """
      i = 0
      for E in self.elements:
         aux = []
         for e in E.onsite:
            aux.append(i)
            i += 1
         E.indices = aux
   def get_sublattice(self,subs=[]):
      """
        This function adds a Sublattice attribute to the base elements
      """
      LG.info('Calculating sublattice for each element')
      #subs = geo.fsublattice(self.bonds[0][0])
      if len(subs) == 0: subs = geo.sublattice(self.bonds[0][0])
      else:
         if len(subs) != len(self.pos):
            LG.critical('Different number of atoms and sublattice')
      for i in range(len(subs)):
         self.elements[i].sublattice = subs[i]
      self.sublattices = np.array(subs)
      for i in range(len(self.basis)):
         self.basis[i] += (self.sublattices[self.basis[i][0]],)
   def get_layer(self):
      """
        This function adds a sublattice attribute to the base elements
      """
      LG.info('Calculating Layer for each element')
      lays = geo.layer(self.pos)
      lays = 2*np.array(lays) - 1   #TODO general map to interval [-1,1]
      for i in range(len(lays)):
         self.elements[i].layer = lays[i]
      self.layers = lays
   def get_neig(self,nvec=5,fol='./'):
      LG.info('Looking for neighbors')
      self.bonds = []
      for A in geo.fneig(self.pos,self.latt,fol=fol):
         LG.info('Neighbor: %s'%(A[2]))
         v0,v1 = A[0][0], A[0][1]
         data = [1 for _ in range(len(v0))]
         a = coo_matrix( (data,(v0,v1)), shape=(len(self.pos),len(self.pos)) )
         self.bonds.append((a,A[1]))
      # TODO check symmetic intra
      #M = self.bonds[0][0].todense()
      #if not np.allclose(M.transpose(), M):
      #   msg = 'On site neighbouring matrix is not symmetric.'
      #   msg += ' Try base.get_neig(nvec) with nvec bigger than %s'%(nvec)
      #   LG.critical(msg)
      #   exit('Aborted @ get_neig method')
      LG.info('Neighbors done')
      return self.bonds
   @log_help.log2screen(LG)
   def vacancy(self,N=1,d=None,alpha=0.,ind=None,inf=1000000.):
      """
        This function chooses N atoms in the center of the cell and adds an
        infinite on-site energy to them killing the hoppings
        N: number of vacancies to add
        d: distance between vacancies
        alpha: angle (respect to X axis) of the vector between vacancies
        ind: Not implemented
      """
      ## Get atoms in layer 1 and not connected in layer 0
      sub_atsA =np.where((self.layers==1)&(self.sublattices=='A'),self.inds,-1)
      sub_atsA = sub_atsA[sub_atsA>0]
      lena = len(self.find_neig(sub_atsA[0]))
      sub_atsB =np.where((self.layers==1)&(self.sublattices=='B'),self.inds,-1)
      sub_atsB = sub_atsB[sub_atsB>0]
      lenb = len(self.find_neig(sub_atsB[0]))

      ## sub_ats contains the hollow atoms in layer 1
      if lena < lenb: sub_ats = sub_atsA  # indices of hollow atoms
      else: sub_ats = sub_atsB  # indices of hollow atoms

      ## Select atoms for defects
      if N == 1: ## 1 defect
         C = np.mean(self.pos[sub_ats],axis=0)
         ind = geo.snap(C,self.pos[sub_ats])
         indices = [sub_ats[ind]]
         LG.warning('Requested-dist/Real-dist: 0.0/0.0')
         LG.warning('Requested-angle/Real-angle: 0.0/0.0')
      elif N == 2: ## 2 defects
         if d == None: d = np.sqrt(3)*(np.max(self.x)-np.min(self.x))/12
         msg  = 'Including two vacancies %s Ansg apart'%(d)
         msg += ' with an angle %s'%(alpha)
         LG.info(msg)
         # returns indices referred to subset
         inds = geo.defects(self.pos[sub_ats],d=d,alpha=alpha)
         indices = [sub_ats[i] for i in inds]
         msg = 'Vacancies placed at'
         for i in indices:
            msg += ' %s'%(i)
         LG.info(msg)
         rd = self.pos[indices[0]] - self.pos[indices[1]]
         ra = np.degrees(np.arctan2(rd[1],rd[0]))
         rd = np.linalg.norm(rd)
         LG.warning('Requested-dist/Real-dist: %s/%s'%(d,rd))
         LG.warning('Requested-angle/Real-angle: %s/%s'%(alpha,ra))
      for i in indices:
         LG.info('Changing onsite of atom: %s'%(i))
         aux = deepcopy(self.elements[i].onsite)
         for k,v in aux.items():
            aux[k] = inf
         self.elements[i].onsite = aux
      self.basis = []
      for E in self.elements:
         for o in E.orbitals:
            ID = (E.place,E.element,o)
            self.basis.append(ID)
      return indices
   def find_neig(self,ind):
      """ returns the index of the neighbouring atoms of atom ind """
      try: self.bonds
      except AttributeError: self.get_neig()
      if not hasattr(ind,'__iter__'): ind = [ind]
      elif not hasattr(ind,'__getitem__'): ind = [ind]
      ret = []
      for i in ind:
         M = self.bonds[0][0]
         r,c = M.row,M.col
         ret.append( c[r==i] )
      if len(ret) == 1: return ret[0]
      else: return ret
   def save(self,f_ind='base.basis',f_xyz='base.xyz'):
      self.update_basis()
      f = open(f_ind,'w')
      f.write('#ind   n_atom   atom   orb   sublatt\n')
      for i in range(len(self.basis)):
         #e = self.basis[i]
         f.write(str(i))
         for e in self.basis[i]:
            f.write('   '+str(e))
         f.write('\n')
#+'   '+str(e[0])+'   '+str(e[1])+'   '+str(e[2])+'\n')
      f.close()
      self.save_xyz(f_xyz)
   def save_xyz(self,fname='base.xyz'):
      """ Save the base positions and lattice vector to a xyz file """
      pos = [E.position for E in self.elements]
      ats = [E.element for E in self.elements]
      IO.pos2xyz(pos,self.latt,at=ats,sub=self.sublattices,fname=fname)


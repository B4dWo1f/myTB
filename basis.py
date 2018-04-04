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
import algebra as alg
## LOG
import logging
import log_help
LG = logging.getLogger(__name__)



class Base_Element(object):
   #def __init__(self,n=0,elem='',onsites={},pos=np.array(())):
   def __init__(self,n=0,elem='',atoms={},pos=np.array(())):
      self.place = n
      self.element = elem
      onsites = atoms[elem]
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
   def __init__(self, elements=[],latt=[],atoms={},cent=True,dospin=False):
      """
        Capital letter attirbutes must have the dimension of the hamiltonian
        Lower letter attirbutes must have the dimension of the number of atoms
      """
      self.DOspin = dospin
      self.elements = elements     # List of Base_Elements
      self.atoms = atoms           # Dictionary of on-site energies
      ## Lattice vectors
      self.latt = np.array(latt)
      self.recip = np.array(geo.reciprocal(latt))
      ## Fix geometry
      if cent: sc = self.center()    # center unit cell in (0,0,0)
      self.evaluate()
   def evaluate(self):
      self.Natoms = len(self.elements)  # Number of atoms in the system
      ats,orbs,pos = [],[],[]   #TODO adapt to include the sublattice
      for e in self.elements:
         ats.append(e.element)
         orbs.append(e.orbitals)
         pos.append(e.position)
      inds,ATS = [],[]
      for i in range(len(orbs)):
         for _ in orbs[i]:
            inds.append(i)
            ATS.append(ats[i])
      aux_orbs = []
      for o in orbs:
         for io in o:
            aux_orbs.append(io)
      ## Everything to numpy array
      pos = np.array(pos)
      ats = np.array(ats)
      orbs = np.array(aux_orbs)
      ATS = np.array(ATS)
      inds = np.array(inds)
      aux_inds = np.array(range(len(inds)))
      ## Useful attributes
      # Number of atoms
      self.pos = pos
      self.ats = ats
      # Number of orbitals
      self.ndim = len(ATS)  # Hamiltonian dimension
      self.ATS = np.array(ATS)
      self.ORBS = np.array(orbs)
      self.INDS = np.array(inds)          #[0,0,0,1,1,1,.... Nats, Nats]
      self.AUX_INDS = np.array(aux_inds)  #[0,1,2,...Norbs]
      self.SPIN = np.array([1 for _ in self.ORBS])
      ## Positions of the atoms
      self.x = self.pos[:,0]
      self.y = self.pos[:,1]
      self.z = self.pos[:,2]
      del pos,ats,orbs,ATS,inds,aux_inds
   def copy(self): return deepcopy(self)
   def center(self):
      LG.debug('Centering cell')
      pos = np.array([e.position for e in self.elements])
      v = np.mean(pos,axis=0)
      LG.info('Center cell by vector: (%.3f,%.3f,%.3f)'%(v[0],v[1],v[2]))
      for i in range(len(self.elements)):
         self.elements[i].position -= v
      self.pos = np.array([e.position for e in self.elements])
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
      if len(subs) == 0:
         subs = geo.sublattice(self.bonds[0][0])
      else:
         if len(subs) != self.Natoms:
            LG.critical('Different number of atoms and sublattice')
            LG.info('Trying to calculate again the sublattice')
            subs = geo.sublattice(self.bonds[0][0])
      if not geo.check_sublattice(self.bonds[0][0],subs):
         LG.warning('Sublattice was wrong. Calculated again')
         subs = geo.sublattice(self.bonds[0][0])
      self.subs = np.array(subs)
      for i in range(len(self.subs)):
         self.elements[i].sublattice = self.subs[i]
      self.SUBS = []
      for i in range(len(self.INDS)):
         self.SUBS.append(self.subs[self.INDS[i]])
      self.SUBS = np.array(self.SUBS)
   def get_layer(self):
      """
        This function adds a sublattice attribute to the base elements
      """
      LG.info('Calculating Layer for each element')
      try:
         lays = [e.layer for e in self.elements]
      except AttributeError: 
         lays = geo.layer(self.pos)
         lays = 2*np.array(lays) - 1   #TODO general map to interval [-1,1]
      self.layers = np.array(lays)
      for i in range(len(lays)):
         self.elements[i].layer = lays[i]
      for i in range(len(self.layers)):
         self.elements[i].layer = self.layers[i]
      self.LAYS = []
      for i in range(len(self.INDS)):
         self.LAYS.append(self.layers[self.INDS[i]])
   def get_neig(self,nvec=5,fol='./'):
      LG.info('Looking for neighbors')
      self.bonds = []
      for A in geo.fneig(self.pos,self.latt,fol=fol):
         LG.info('Neighbor: %s'%(A[2]))
         v0,v1 = A[0][0], A[0][1]
         data = [1 for _ in range(len(v0))]
         a = coo_matrix( (data,(v0,v1)), shape=(self.Natoms,self.Natoms) )
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
   def adatom(self,l=1,N=1,at='H',dummy=False,inf=1e6):
      """
        This function chooses N atoms in the center of the cell and adds an
        infinite on-site energy to them killing the hoppings
        N: number of vacancies to add
        d: distance between vacancies
        alpha: angle (respect to X axis) of the vector between vacancies
        ind: Not implemented
      """
      ## Get atoms in layer 1 and not connected in layer 0
      #print(len(self.layers),len(self.subs),len(self.INDS))
      l = max(set(self.layers))
      LG.info('Vacancies introduced in layer: %s'%(l))
      aux = range(max(self.INDS)+1)   #+1 because python starts in 0
      sub_atsA =np.where((self.layers==l)&(self.subs==1),aux,-1)  #A
      sub_atsA = sub_atsA[sub_atsA>0]
      lena = len(self.find_neig(sub_atsA[0]))
      sub_atsB =np.where((self.layers==l)&(self.subs==-1),aux,-1) #B
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
         LG.info('Changing onsite of atom: %s'%(indices[0])) #XXX XXX XXX
         ## TODO generalize for list of adatoms
         pl = self.elements[-1].place + 1
         r = self.pos[indices[0]] + np.array([0,0,1.4])
         LG.info('Adatom in position: %.3f, %.3f, %.3f'%(r[0],r[1],r[2]))
         la = self.elements[-1].layer  # TODO think + 1
         sub_dict = {'A':1,'B':-1}
         tcid_bus = {1:'A',-1:'B'}
         su = -1 * self.elements[-1].sublattice
         if dummy:
            fake_atoms = deepcopy(self.atoms)
            for k,v in fake_atoms[at].items():
               fake_atoms[at][k] = inf
            self.elements.append( Base_Element(pl,at,fake_atoms,r) )
         else: self.elements.append( Base_Element(pl,at,self.atoms,r) )
         self.elements[-1].layer = la  #TODO fix
         self.elements[-1].sublattice = su  #TODO fix
         ## Update basis attributes
         self.evaluate()
         self.get_neig(fol='')   #TODO Fix to read from file?
         self.get_layer()
         self.get_sublattice()
      return indices
   @log_help.log2screen(LG)
   def vacancy(self,l=1,N=1,d=None,alpha=0.,ind=None,inf=1e9):
      """
        This function chooses N atoms in the center of the cell and adds an
        infinite on-site energy to them killing the hoppings
        N: number of vacancies to add
        d: distance between vacancies
        alpha: angle (respect to X axis) of the vector between vacancies
        ind: Not implemented
      """
      ## Get atoms in layer 1 and not connected in layer 0
      #print(len(self.layers),len(self.subs),len(self.INDS))
      l = max(set(self.layers))
      LG.info('Vacancies introduced in layer: %s'%(l))
      aux = range(max(self.INDS)+1)   #+1 because python starts in 0
      sub_atsA =np.where((self.layers==l)&(self.subs==1),aux,-1)  #A
      sub_atsA = sub_atsA[sub_atsA>0]
      lena = len(self.find_neig(sub_atsA[0]))
      sub_atsB =np.where((self.layers==l)&(self.subs==-1),aux,-1) #B
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
      #self.update_basis()
      f = open(f_ind,'w')
      f.write('#ind   n_atom   atom   orb   spin   sublatt   lay\n')
      for i in range(len(self.INDS)):
         f.write(str(self.AUX_INDS[i])+'   '+str(self.INDS[i]))
         f.write('   '+str(self.ATS[i])+'   '+str(self.ORBS[i]))
         f.write('   '+str(self.SPIN[i])+'   '+str(self.SUBS[i]))
         f.write('   '+str(self.LAYS[i])+'\n')
      f.close()
      self.save_xyz(f_xyz)
   def save_xyz(self,fname='base.xyz'):
      """ Save the base positions and lattice vector to a xyz file """
      pos = [E.position for E in self.elements]
      ats = [E.element for E in self.elements]
      IO.write.xyz(pos,self.latt,at=ats,sub=self.subs,fname=fname)
   def dospin(self):
      LG.info('Modifying basis for spin')
      self.LAYS = alg.m2spin(self.LAYS)
      self.ORBS = alg.m2spin(self.ORBS,Delt='')
      self.SUBS = alg.m2spin(self.SUBS)
      self.ATS = alg.m2spin(self.ATS)
      self.AUX_INDS = alg.m2spin(self.AUX_INDS)
      self.INDS = alg.m2spin(self.INDS)
      self.SPIN = np.array([(-1)**i for i in range(len(self.ORBS))])
      if len(self.LAYS) != len(self.ORBS) != len(self.SUBS) !=\
         len(self.ATS) != len(self.AUX_INDS) != len(self.INDS):
         LG.critical('Problem Spin-doubling')
         print('Layers:', len(self.LAYS))
         print('Orbitals:',len(self.ORBS))
         print('Sublattice:',len(self.SUBS))
         print('Atoms:',len(self.ATS))
         print('Indices',len(self.INDS))
         print('Aux Indices:',len(self.AUX_INDS))
         exit()
      self.ndim = len(self.ATS)  # Hamiltonian dimension
      LG.info('New basis dimension: %s'%(self.ndim))

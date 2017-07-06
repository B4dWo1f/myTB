#!/usr/bin/python3
# -*- coding: UTF-8 -*-

from configparser import ConfigParser, ExtendedInterpolation
from os.path import expanduser
import os
import numpy as np
## LOG
import logging
import log_help
LG = logging.getLogger(__name__)



class fol_param(object):
   def __init__(self,out='./out/',slf='./slf/',ham='./ham/'):
      self.out = out
      self.slf = slf
      self.ham = ham
      self.create()
   def __str__(self):
      msg = 'Folder structure\n'
      msg += '        Outputs: %s\n'%(self.out)
      msg += '  Self-Energies: %s\n'%(self.slf)
      msg += '   Hamiltonians: %s\n'%(self.ham)
      return msg
   def create(self):
      for f in [self.out,self.slf,self.ham]:
         LG.warning('Creating folder: %s'%(f))
         os.system('mkdir -p %s'%(f))

class ham_param(object):
   def __init__(self,lzee=np.array((0,0,0)),lSO=0.,lmass=0.,lelec=0.):
      self.lzee = lzee
      self.lSO = lSO
      self.lmass = lmass
      self.lelec = lelec
   def __str__(self):
      msg = 'Hamiltonian parameters\n'
      msg += '             Zeeman: '
      msg += '  (%s, %s, %s)\n'%(self.lzee[0],self.lzee[1],self.lzee[2])
      msg += '         Spin Orbit: %s\n'%(self.lSO)
      msg += '  Sublatt imbalance: %s\n'%(self.lmass)
      msg += '     Electric Field: %s\n'%(self.lelec)
      return msg

class calc_param(object):
   def __init__(self,bands=False,spectrum=False,nk=100):
      self.bands = bands
      self.spectrum = spectrum
      self.nk = nk
   def __str__(self):
      msg = 'Calculations\n'
      msg += '     Bands: %s\n'%(self.bands)
      msg += '  Spectrum: %s\n'%(self.spectrum)
      return msg

class vacancy_param(object):
   def __init__(self,N,d,alpha):
      self.N = N
      self.d = d
      self.alpha = alpha
   def __str__(self):
      if self.N == 0: msg = 'No vacancies to be introduced\n'
      elif self.N == 1:
         msg = '%s Vacancy to be introduced\n'%(self.N)
      elif self.N > 1:
         msg = '%s Vacancies to beintroduced'%(self.N)
         msg += ' separated %s Angs with an angle %s\n'%(self.d,self.alpha)
      return msg

class sys_param(vacancy_param):
   def __init__(self,xyz_file,pasivate,dists,vacancy_param): #,defect,dists):
      self.xyz_file = xyz_file
      self.pasivate = pasivate
      #self.defect = defect
      self.dists = dists
      self.vac = vacancy_param
   def __str__(self):
      msg = 'Physical system parameters\n'
      msg += '      Pasivate Unit Cell: %s\n'%(self.pasivate)
      #msg += '          Type of defect: %s\n'%(self.defect)
      msg += '                XYZ file: %s\n'%(self.xyz_file)
      msg += '  Distance between atoms:\n'
      for key, value in self.dists.items():
         msg += '     %s: %s\n'%(key,value)
      return msg


def setup(fname='SK1.ini'):
   """
    Parse an ini file returning the parameter classes
    TODO: Log this function
   """
   config = ConfigParser(inline_comment_prefixes='#')
   config._interpolation = ExtendedInterpolation()
   config.read(fname)

   out = expanduser(config['I/O']['output'])
   slf = expanduser(config['I/O']['selfes'])
   ham = expanduser(config['I/O']['hamils'])
   FP = fol_param(out,slf,ham)

   lzee = np.array(eval(config['hamiltonian']['lzee']))
   lSO = eval(config['hamiltonian']['lSO'])
   lmass = eval(config['hamiltonian']['lmass'])
   lelec = eval(config['hamiltonian']['lelec'])
   HP = ham_param(lzee,lSO,lmass,lelec)

   bands = eval(config['calculations']['bands'].capitalize())
   spectrum = eval(config['calculations']['spectrum'].capitalize())
   nk = int(config['calculations']['nk'])
   CP = calc_param(bands,spectrum,nk)

   xyz_file = config['system']['xyz_file']
   pasivate = eval(config['system']['pasivate'].capitalize())
   #defect = config['system']['defect']
   dist = eval(config['system']['dist'])
   ##SP = sys_param(xyz_file,pasivate,dist)
   N = int(config['vacancy']['N'])
   d = float(config['vacancy']['d'])
   alpha = float(config['vacancy']['alpha'])
   vp = vacancy_param(N,d,alpha)
   SP = sys_param(xyz_file,pasivate,dist,vp)

   keys,values = [],[]
   for key in config['atoms']:
      keys.append(key.capitalize())
      values.append(eval(config['atoms'][key]))
   atoms = dict(zip(keys, values))

   keys,values = [],[]
   for key in config['hopping']:
      a = eval(config['hopping'][key])
      keys.append('-'.join(sorted(key.title().split('-'))) )
      values.append(a)
   hoppings = dict(zip(keys, values))
   #### FIX of the hoppings
   #    C-C = {Vsss:...} ---> C-C = {1:{Vsss:...,Vsps:...}}
   for k,v in hoppings.items():
      #k should be a string the hopping type (C-C, C-H, Mo-S...)
      #v should be a dictionary with the SK hopping parameters of the form
      #   v = {1:{params...}, 2:{params}}
      try:
         for ik,iv in v.items():
            if type(ik) != int:
               aux = {1:hoppings[k]}
               hoppings[k] = aux
               break
            else: pass
      except AttributeError: pass
   return FP,HP,CP,SP,atoms,hoppings

import os
def compile_fortran(fname):
   root_fname = '.'.join(fname.split('.')[0:-1])
   LG.info('Compilando fortran con f2py')
   os.system('f2py -c -m %s %s'%(root_fname,fname))
   LG.info('   ...Compilado fortran con f2py')
   os.system('cp %s .%s'%(fname,fname))
   LG.warning('Hidden copy to avoid re-compiling')


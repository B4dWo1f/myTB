#!/usr/bin/python3
# -*- coding: UTF-8 -*-

from configparser import ConfigParser, ExtendedInterpolation
from os.path import expanduser
import os
import numpy as np
import IO
## LOG
import logging
import log_help
LG = logging.getLogger(__name__)



class fol_param(object):
   def __init__(self,out='./out/',slf='./slf/',ham='./ham/'):
      self.out = out
      self.slf = slf
      self.ham = ham
      comm = []
      for o,s,h in zip(out.split('/'),slf.split('/'),ham.split('/')):
         if o==s==h: comm.append(o)
         else:break
      self.root = '/'.join(comm)
      self.create()
   def __str__(self):
      msg = 'Folder structure\n'
      msg += 'root foder: %s\n'%(self.root)
      msg += '        Outputs: %s\n'%(self.out.replace(self.root,''))
      msg += '  Self-Energies: %s\n'%(self.slf.replace(self.root,''))
      msg += '   Hamiltonians: %s\n'%(self.ham.replace(self.root,''))
      return msg
   def create(self):
      for f in [self.out,self.slf,self.ham]:
         LG.debug('Creating folder: %s'%(f.replace(self.root,'')))
         os.system('mkdir -p %s'%(f))

class ham_param(object):
   def __init__(self,lzee=np.array((0,0,0)),lSO=0.,lmass=0.,lelec=0.,lrashba=0.):
      self.lzee = lzee
      self.lSO = lSO
      self.lmass = lmass
      self.lelec = lelec
      self.lrashba = lelec * lrashba #XXX testing
                                     # From Phys. Rev. B 74, 165310 (2006)
                                     #    50V/300nm --> 0.011e-3 eV
   def __str__(self):
      msg = 'Hamiltonian parameters\n'
      msg += '             Zeeman: '
      msg += '  (%s, %s, %s)\n'%(self.lzee[0],self.lzee[1],self.lzee[2])
      msg += '         Spin Orbit: %s\n'%(self.lSO)
      msg += '  Sublatt imbalance: %s\n'%(self.lmass)
      msg += '     Electric Field: %s\n'%(self.lelec)
      return msg

class calc_param(object):
   def __init__(self,bands=False,spectrum=False,DOS=False,local=False,
                                      nk=100,ns=9,ndos=150,nddos=20,nkdos=100):
      self.bands = bands
      self.spectrum = spectrum
      self.dos = DOS
      self.local = local
      self.ns = ns
      self.nk = nk
      self.ndos = ndos    # energies in which calculate the DOS
      self.nddos = nddos  # Number of eigvalue to take into account
      self.nkdos = nkdos  # Number of k points (per latt vec)
   def __str__(self):
      msg = 'Calculations\n'
      msg += '  Spectrum: %s\n'%(self.spectrum)
      msg += '     Bands: %s\n'%(self.bands)
      msg += '       DOS: %s\n'%(self.dos)
      msg += '      lDOS: %s\n'%(self.local)
      return msg

class adatom_param(object):
   def __init__(self,N,sp3=0.0,hollow=True):
      self.N = N
      self.sp3 = sp3
      self.hollow = hollow
   def __str__(self):
      if self.N == 0: msg = 'No adatoms to be introduced\n'
      elif self.N == 1:
         msg = '%s adatom to be introduced\n'%(self.N)
      return msg

class vacancy_param(object):
   def __init__(self,N,d,alpha,hollow=True):
      self.N = N
      self.d = d
      self.alpha = alpha
      self.hollow = hollow
   def __str__(self):
      if self.N == 0: msg = 'No vacancies to be introduced\n'
      elif self.N == 1:
         msg = '%s Vacancy to be introduced\n'%(self.N)
      elif self.N > 1:
         msg = '%s Vacancies to beintroduced'%(self.N)
         msg += ' separated %s Angs with an angle %s\n'%(self.d,self.alpha)
      return msg

class sys_param(vacancy_param):
   def __init__(self,xyz_file,pasivate,dists,vacancy_param,adatom_param,DOspin,\
                     force0D,periodic):
      self.xyz_file = xyz_file
      self.pasivate = pasivate
      #self.defect = defect
      self.dists = dists
      self.ada = adatom_param
      self.vac = vacancy_param
      self.force0D = force0D
      self.pbc = periodic
      self.DOspin = DOspin
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
          try/except for format errors
   """
   config = ConfigParser(inline_comment_prefixes='#')
   config._interpolation = ExtendedInterpolation()
   config.read(fname)

   lzee  = config['hamiltonian']['lzee'].replace(')','').replace('(','')
   lzee  = np.array( list(map(float,lzee.split(','))) )
   lSO   = float(config['hamiltonian']['lSO'])
   lmass = float(config['hamiltonian']['lmass'])
   lelec = float(config['hamiltonian']['lelec'])
   HP = ham_param(lzee,lSO,lmass,lelec)

   bands = eval(config['calculations']['bands'].capitalize())
   spectrum = eval(config['calculations']['spectrum'].capitalize())
   try: DOS = eval(config['calculations']['dos'].capitalize())
   except NameError: DOS = config['calculations']['dos'].lower()
   DOS_local = eval(config['calculations']['local'].capitalize())
   Ndos = eval(config['calculations']['Ndos'])
   Nddos = eval(config['calculations']['Nddos'])
   Nkdos = eval(config['calculations']['Nkdos'])
   ns = int(config['calculations']['Ns'])
   nk = int(config['calculations']['nk'])
   CP = calc_param(bands=bands,spectrum=spectrum,DOS=DOS,local=DOS_local,
                                 ns=ns,nk=nk,ndos=Ndos,nddos=Nddos,nkdos=Nkdos)

   sys = config['system']['sys']
   n = int(config['system']['n'])
   l = int(config['system']['l'])
   xyz_file = config['system']['xyz_file']
   #xyz_dir = config['system']['xyz_dir']    #XXX future implementation
   pasivate = eval(config['system']['pasivate'].capitalize())
   dist = eval(config['system']['dist'])
   force0D = eval(config['system']['force0D'].capitalize())
   periodic = eval(config['system']['periodic'].capitalize())
   DOspin = eval(config['system']['DOspin'].capitalize())

   Nv = int(config['vacancy']['N'])
   d = eval(config['vacancy']['d'])
   alpha = float(config['vacancy']['alpha'])
   hollow = eval(config['vacancy']['hollow'])
   vp = vacancy_param(Nv, d, alpha, hollow=hollow)

   Na = int(config['adatom']['na'])
   sp3 = float(config['adatom']['sp3'])
   hollow = bool(config['adatom']['hollow'])
   ap = adatom_param(Na,sp3,hollow)

   ## System parameters
   SP = sys_param(xyz_file,pasivate,dist,vp,ap,DOspin,force0D,periodic)

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
   HP.hoppings = hoppings
   #### FIX the spin doubling if needed but not requested
   if np.linalg.norm(HP.lzee)!= 0. or HP.lSO != 0.:
      LG.warning('Need spin doubling but it wasn\'t requested')
      SP.DOspin = True

   tail = ''
   try: tail += '%sorb/'%len(atoms['C'])
   except KeyError: pass
   tail += config['system']['sys'] + '/'
   tail += 'n%s_l%s/'%(config['system']['n'],config['system']['l'])
   tail += 'nv%s_na%s/'%(Nv,Na)
   tail += 'd%s/'%(d)
   tail += 'alpha%s/'%(alpha)
   tail += 'e%s/'%(lelec)

   basedir = expanduser(config['I/O']['basedir'])
   if basedir[-1] != '/': basedir += '/'
   out = basedir + 'OUTS/'+ tail
   slf = basedir + 'SELFEs/' + '/'.join(tail.split('/')[:-4]) + '/'
   ham = basedir + 'HAMILs/' + '/'.join(tail.split('/')[:-4]) + '/'
   FP = fol_param(out,slf,ham)

   return FP,HP,CP,SP,atoms #,hoppings


#def get_sys(sys,n,l,pasiv=False,xyz_dir='../cells'):
#   if pasiv: xyz_file = xyz_dir + '%s_n%s_l%s_H.xyz'%(sys,n,l)
#   else: xyz_file = xyz_dir + '%s_n%s_l%s.xyz'%(sys,n,l)
#   try: ## Try to load XYZ file
#      ats,pos,latt,sub = IO.read.xyz(xyz_file)
#   except FileNotFoundError:
#      import islands
#      #sysname = {'ac':islands.armchair, 'simple':islands.simple}
#      #ats,pos,latt,sub = sysname[sys](n)
#      #if pasiv: A.pasivate()
#      #if l>1: A.multilayer(l)
#      cell = islands.get_cell(sys,n,l,pasiv)
#      print(cell)
#      exit()
#   return ats,pos,latt,sub

import os
def compile_fortran(fname):
   def doit(fname):
      LG.debug('Backup file (.%s) not found'%(fname))
      LG.info('Compilando fortran con f2py')
      os.system('f2py3.5 -c -m %s %s'%(root_fname,fname))
      LG.info('   ...Compilado fortran con f2py')
      os.system('cp %s .%s'%(fname,fname))
      LG.warning('Hidden copy to avoid re-compiling')
   root_fname = '.'.join(fname.split('.')[0:-1])
   if not os.path.exists('.%s'%(fname)): doit(fname)
   else:
      LG.debug('Backup file (.%s) is present'%(fname))
      diff_for = os.popen('diff %s .%s'%(fname,fname)).read()
      diff_for = diff_for.lstrip().rstrip()
      diff_for.splitlines()
      so = os.popen('ls %s.*so 2> /dev/null'%(root_fname)).read()
      if len(diff_for) > 1 or len(so) == 0: doit(fname)
      else: LG.info('%s is already compiled'%(fname))

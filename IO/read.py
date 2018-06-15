#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import numpy as np
from scipy.sparse import coo_matrix
import sys
import log_help
import logging
LG = logging.getLogger(__name__)

@log_help.log2screen(LG)
def xyz(archivo):
   """
     Reads the lattice information from an extended xyz file. The file is
    assumed to have the following structure:
        N atoms
        [latt vec 1][latt vec 2]...   # eg: [1,0,0][0,1,0]
        C   0   0   0   1/-1   # atom   X   Y   Z   sublattice
     If the lattice vectors are not specified it will return an empty list, so
    the program will calculate the system considering it an island.
   """
   ## Number of atoms
   lines = open(archivo,"r").readlines()
   nat = int(lines[0]) # number of atoms
   LG.debug('Expecting %s atoms'%(nat))
   ## Lattice vectors
   try:
      lines[1] = lines[1].split('#')[0] # clean comments
      vecs = lines[1].replace(' ','').lstrip().rstrip().split('][')
      vecs = [v.replace('[','').replace(']','') for v in vecs]
      vecs = [np.array(list(map(float,v.split(',')))) for v in vecs]
   except:
      LG.warning('Unable to get lattice vectors. Empty list.')
      vecs = []
   ## Atoms, atomic postions and sublattice
   # atoms
   atoms = np.loadtxt(archivo,skiprows=2,usecols=(0,),dtype=bytes)
   atoms = np.array([str(a,'utf-8') for a in atoms])
   # positions & sublattice
   try: 
      pos = np.loadtxt(archivo,skiprows=2,usecols=(1,2,3,4))
      sub = np.asarray(pos[:,-1],int)
      pos = pos[:,0:3]
   except:
      pos = np.loadtxt(archivo,skiprows=2,usecols=(1,2,3))
      sub = np.array([])
   if atoms.shape[0] != pos.shape[0]:
      LG.critical('Different number of atoms and positions')
   if pos.shape[0] != sub.shape[0]:
      LG.critical('Sublayer not read from file')
   LG.info('Read %s atoms, with %s lattice vectors'%(len(atoms),len(vecs)))
   return atoms,pos,vecs,sub

def bands(fname):
   X,Y,Z = np.loadtxt(fname,unpack=True)
   return X,Y,Z

def read_matv(fname):
   """
   read from a plain text with format
   #Ndim  v
   i   j   real   imag
   """
   n = open(fname,'r').readlines()[0].lstrip().rstrip().replace('#','')
   n = int(n.split()[0])
   v = np.array(eval(''.join(n.split()[1:])))  # Useless?
   I,J,R,Im = np.loadtxt(fname,skiprows=1,unpack=True)
   if all(Im==0):
      data = R
      ty = float
   else:
      data = R + 1j*Im
      ty=complex
   M = coo_matrix( (data, (I, J)), shape=(n,n),dtype=ty)
   return M,v


def json_read(fname):
   """
     Returns the dictionary contained in the JSON file
   """
   lines = open(fname,'r').readlines()
   if len(lines) == 1:
      return eval(lines[0].lstrip().rstrip())
   else:
      sys.exit('ERROR while reading JSON file (%s)'%(fname))


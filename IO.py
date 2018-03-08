#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import numpy as np
import sys
import log_help
import logging
LG = logging.getLogger(__name__)

@log_help.log2screen(LG)
def xyz(archivo):
   """
     Reads the lattice information from an ENHANCED xyz file. The file is
    assumed to have the following structure:
        N atoms
        [latt vec 1][latt vec 2]...   # eg: [1,0,0][0,1,0]
        C   0   0   0   # atom   X   Y   Z
     If the lattice vectors are not specified it will return a 0 vector, so
    the program will calculate the system considering it an island.
   """
   #LG = logging.getLogger('IO.xyz')
   lines = open(archivo,"r").readlines()
   nat = int(lines[0]) # number of atoms
   LG.debug('Expecting %s atoms'%(nat))
   try:    ## try to get lattice vectors
      lines[1] = lines[1].split('#')[0] # clean comments
      vecs = lines[1].replace(' ','').lstrip().rstrip().split('][')
      vecs = [v.replace('[','').replace(']','') for v in vecs]
      vecs = [np.array(list(map(float,v.split(',')))) for v in vecs]
   except:
      LG.warning('Unable to get lattice vectors. Empty list.')
      vecs = []
   atoms = np.loadtxt(archivo,skiprows=2,usecols=(0,),dtype=str)
   atoms = np.array([a.replace('b\'','') for a in atoms])
   atoms = np.array([a.replace('\'','') for a in atoms])
   pos = np.loadtxt(archivo,skiprows=2,usecols=(1,2,3,4))
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


def pos2xyz(pos,latt,at='C',sub=[],fname='lattice.xyz'):
   """
     at has to be a string or a list/array
   """
   #LG = logging.getLogger('IO.pos2xyz')
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
         f.write(a+'   %s   %s   %s'%(r[0],r[1],r[2]))
         try: f.write('   '+sub[i]+'\n')
         except: f.write('\n')
      LG.debug('Written all the atomic positions')
   f.close()

def mat(fname,M,v=(0,0,0)):
   M = coo_matrix(M)    #XXX check type first
   dim = M.shape[0]
   f = open(fname,'w')
   f.write('# '+str(dim)+' (%s,%s,%s)\n'%(v[0],v[1],v[2]))
   for r,c,d in zip(M.row,M.col,M.data):
      f.write(str(r)+'   '+str(c)+'   '+\
      str(d.real)+'   '+str(d.imag)+'\n')
   f.close()

from scipy.sparse import coo_matrix
def read_matv(fname,ty=complex):
   """
   read from a plain text with format
   #Ndim  v
   i   j   real   imag
   """
   n = open(fname,'r').readlines()[0].lstrip().rstrip().replace('#','')
   v = np.array(eval(''.join(n.split()[1:])))
   n = int(n.split()[0])
   I,J,R,Im = np.loadtxt(fname,skiprows=1,unpack=True)
   if ty == complex: data = R + 1j*Im
   else: data = R
   M = coo_matrix( (data, (I, J)), shape=(n,n),dtype=ty)
   return M,v

def save_matrix(fname,M,sbin=False,fmt='%.1d'):
   """
      fname: [str] name of the file
      M: [np.matrix or array] Matrix to be saved
      fmt: [str] format to write the matrix. %.1d = int
                                             %n.mf = n tot space per entry with
                                                     m decimal digits
      sbin: [bool] overrides everything and saves using numpy.save
   """
   #LG = logging.getLogger('IO.save_matrix')
   if sbin: np.save(fname,M)
   else: np.savetxt(fname,M,fmt=fmt)



def decide(v):
   """
     Returns the proper form of a given value depending on its type
     For lists returns the string '[A,B,...]', for numbers just numbers...
   """
   if isinstance(v,list): return '['+','.join([decide(iv) for iv in v])+']'
   elif isinstance(v,str): return '\''+v+'\''
   elif isinstance(v,float) or isinstance(v,int): return str(v)
   elif isinstance(v,np.ndarray):
      return 'np.array(['+','.join([decide(iv) for iv in v])+'])'
   elif v == None: return None
   else:
      msg = json_write(v.__dict__,None,False)
      return msg[0:-1]
      print('---')
      print(v)
      print(type(v))
      print('---')
      sys.exit('ERROR trying to convert value to string ---> dict')
      return v


def json_write(dic,fname='file.json',check=True):
   """
     Saves a dictionary into a json file. The final JSON file should contain a
     properly formatted python dictionary
   """
   json = '{'
   for k,v in zip(dic.keys(),dic.values()):
      entry = '\''+k+'\':'
      if decide(v) == None: continue
      entry += decide(v)
      json += entry + ', '
   json = json[:-2]+'}\n'
   if fname == None: return json
   with open(fname,'w') as f: f.write(json)
   if check:
      try: json_read(fname)
      except SyntaxError:
         os.system('cat %s'%(fname))
         os.system('rm %s'%(fname))
         sys.exit('ERROR: when saving JSON file (%s)'%(fname))


def json_read(fname):
   """
     Returns the dictionary contained in the JSON file
   """
   lines = open(fname,'r').readlines()
   if len(lines) == 1:
      return eval(lines[0].lstrip().rstrip())
   else:
      sys.exit('ERROR while reading JSON file (%s)'%(fname))

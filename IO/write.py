#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import numpy as np
from scipy.sparse import coo_matrix
import sys
import log_help
import logging
LG = logging.getLogger(__name__)


def xyz(pos,latt,at='C',sub=[],fname='lattice.xyz'):
   """
     Write the unit cell information to a file. If "at" is a string, the same
     string will be used for every atomic position.
   """
   #LG = logging.getLogger('IO.pos2xyz')
   if isinstance(at,str):
      LG.debug('Only one atom provided. Using %s for all the atoms'%(at))
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
         try: f.write('   %s\n'%(sub[i]))
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

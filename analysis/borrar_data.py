#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import os

for n in reversed([20,30,40,45,50,55,60]):
   print('-'*80)
   print('Doing N=%s'%(n))
   com = 'python3 spectrum_min.py ../../../Documents/data/artificial_lattices/confinement/OUTS/1orb/ac/n%s_l2/nv1_na0/dNone/alpha0.0/ && mv datos.dat datos_n%s.dat'%(n,n)
   print(com)
   os.system(com)

#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import os
import numpy as np


f_in = 'inputs.in'
f_done = 'done.out'
exe = 'main.py'
f_tmp =  'SK1.template.ini'
f_ini = 'SK1.ini'
STOP_file = 'STOP'


def modify(template,inp,modified):
   # alpha, dist, elec
   x,y,z = inp[0],inp[1],inp[2]
   com = 'sed \'s/XXelecXX/%s/\' %s '%(z,template)
   com += ' | sed \'s/XXdXX/%s/\''%(y)
   com += ' | sed \'s/XXalphaXX/%s/\''%(x)
   com += ' > %s'%(modified)
   print('-'*80)
   print(com)
   os.system(com)


while True:
   inputs = open(f_in,'r').read()
   inputs = inputs.splitlines()
   if len(inputs[0]) == 0: break
   inp = list(map(float,inputs[0].split()))
   modify(f_tmp,inp,f_ini)
   com = 'python3 %s %s'%(exe,f_ini)
   print(com)
   os.system(com)
   print('\n\n')
   with open(f_done,'a') as f:
      f.write(inputs[0]+'\n')
   f.close()
   with open(f_in,'w') as f:
      f.write('\n'.join(inputs[1:]))
      f.write('\n')
   f.close()
   if os.path.isfile(STOP_file): exit()

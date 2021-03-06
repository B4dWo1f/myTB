#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import os
import numpy as np
from time import time


f_in = 'inputs.in'
f_done = 'done.out'
exe = 'main.py'
f_tmp =  'SK.template.ini'
f_ini = 'SK.ini'
STOP_file = 'STOP'


def modify(template,inp,modified):
   # alpha, dist, elec
   #x,y,z = inp[0],inp[1],inp[2]
   x = inp[0]
   #y = inp[1]
   com = 'sed \'s/XXelecXX/%s/g\' %s '%(x,template)
   #com += ' | sed \'s/XXVspsXX/%s/\''%(y)
   #com += ' | sed \'s/XXalphaXX/%s/\''%(x)
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
   # com = 'python3 %s %s && echo "hey, done"'%(exe,f_ini)
   com = f'./{exe} {f_ini} && echo "hey, done"'
   print(com)
   told = time()
   os.system(com)
   print(f'\nDone in {time()-told:.3f}s')
   print('--------- main done ---------\n')
   with open(f_done,'a') as f:
      f.write(inputs[0]+'\n')
   f.close()
   with open(f_in,'w') as f:
      f.write('\n'.join(inputs[1:]))
      f.write('\n')
   f.close()
   if os.path.isfile(STOP_file): exit()


#try:
#   import mailator
#   mailator.send_mail('Calculation done',subj='FINISHED')
#except: pass

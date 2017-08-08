#!/usr/bin/python3
# -*- coding: UTF-8 -*-

"""
 Script to run all the known examples
"""

import os

inps = ['G2d1o.ini','G2d1o1v.ini','G2dbi1o.ini','G2dbi1o1ve0.ini',
        'G2dbi1o1ve0.2.ini',
        'G2d4o.ini','G2d4o1a.ini','G2dbi4o.ini','G2dbi4o1ae0.ini',
        'G2dbi4o1ae0.2.ini']
lab = ['Graphene 2D 1orbital',
       'Graphene 2D 1orbital 1 vacancy',
       'Graphene 2D bilayer 1orbital',
       'Graphene 2D bilayer 1orbital 1 vacancy elec = 0',
       'Graphene 2D bilayer 1orbital 1 vacancy elec = 0.2',
       'Graphene 2D 4orbital',
       'Graphene 2D 4orbital 1 adatom',
       'Graphene 2D bilayer 4orbital',
       'Graphene 2D bilayer 4orbital 1 adatom elec = 0',
       'Graphene 2D bilayer 4orbital 1 adatom elec = 0.2']


from time import sleep

w = 80  # width of screen
for i in range(len(inps)):
   l = '---- '+lab[i]+' ----'
   while len(l) < w:
      l = ' '+l+' '
   if len(l) > w: l = l[0:-1]
   print('='*w)
   print(l)
   print('='*w)
   #sleep(10)
   os.system('python3 main.py checks/%s'%(inps[i]))
   os.system('rm -r /tmp/OUTs/ /tmp/HAMILs/ /tmp/SELFEs/')
   print('\n\n')

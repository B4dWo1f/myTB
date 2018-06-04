#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import os
import logging
logging.basicConfig(level=logging.DEBUG,
                  format='%(asctime)s-%(name)-s-%(levelname)-s-%(message)s',
                  datefmt='%Y/%m/%d-%H:%M',
                  filename='install.log', filemode='w')
sh = logging.StreamHandler()
sh.setLevel(logging.INFO)
fmt = logging.Formatter('%(name)s: %(levelname)s %(message)s')
sh.setFormatter(fmt)
logging.getLogger('').addHandler(sh)
LG = logging.getLogger('main')


def compile_fortran(fname):
   """
     This function checks if a certain fortran file has already been compiled
     and compiles it otherwise
   """
   def doit(fname):
      """
        Actual compilation with f2py
      """
      LG.info('Compilando fortran con f2py')
      os.system('f2py3.5 -c -m %s %s'%(root_fname,fname))
      LG.info('   ...Compilado fortran con f2py')
      os.system('cp %s .%s'%(fname,fname))
      LG.warning('Hidden copy to avoid re-compiling')
   root_fname = '.'.join(fname.split('.')[0:-1])
   if not os.path.exists('.%s'%(fname)):
      LG.debug('Backup file (.%s) not found'%(fname))
      doit(fname)
   else:
      LG.debug('Backup file (.%s) is present'%(fname))
      diff_for = os.popen('diff %s .%s'%(fname,fname)).read()
      diff_for = diff_for.lstrip().rstrip()
      diff_for.splitlines()
      so = os.popen('ls %s.*so 2> /dev/null'%(root_fname)).read()
      if len(diff_for) > 1 or len(so) == 0: doit(fname)
      else: LG.info('%s is already compiled'%(fname))
   ## Final check
   com = 'import %s'%(root_fname)
   try: exec(com)
   except ImportError:
      LG.critical('Parallel-Python was not installed')
      exit()
   


msg = 'This script will install all neccessary dependencies to install myTB'

## System dependencies
com = 'sudo apt-get install gfortran gfortran-multilib unzip'
LG.debug(com)


packages = ['python3-numpy','python3-matplotlib','python3-scipy','python3-pint']
names = ['Numpy','MatPlotLib','Scipy','Pint']
for nm,pkg in zip(names,packages):
   com = 'sudo apt-get install %s'%(pkg)
   LG.info('Installing %s'%(nm))
   LG.debug(com)
   os.system(com)

## Extra
url = 'https://www.parallelpython.com/downloads/pp/pp-1.6.6.zip'
com = 'mkdir -p /tmp/PP'
LG.debug(com)
com = 'wget -P /tmp/PP "%s"'%(url)
LG.debug(com)
com = 'unzip /tmp/PP/pp-1.6.6.zip -d /tmp/PP'
LG.debug(com)
com = 'cd /tmp/PP/pp-1.6.6 && sudo python3 setup.py install'
LG.debug(com)
LG.info('Installed Parallel python')


## Checks
# numpy
try: import numpy as np
except ImportError:
   LG.critical('Numpy was not installed')
   exit()

# matplotlib
try: import matplotlib.pyplot as plt
except ImportError:
   LG.critical('MatPlotLib was not installed')
   exit()

# scipy
try: import scipy
except ImportError:
   LG.critical('Scipy was not installed')
   exit()

# pint
try: import pint
except ImportError:
   LG.critical('Pint was not installed')
   exit()

# parallel-python
try: import pp
except ImportError:
   LG.critical('Parallel-Python was not installed')
   exit()


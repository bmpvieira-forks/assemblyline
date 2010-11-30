'''
Created on Nov 29, 2010

@author: mkiyer
'''
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("lib.cNode", ["lib/cNode.pyx"])]

setup(name = 'assemblyline',
      cmdclass = {'build_ext': build_ext},
      ext_modules = ext_modules)

# simple setup.py, adapted from the one that sage
# automatically produces with 'load diricihet_conrey.pyx"
# and the cython documentation

# Build using 'python setup.py'
import distutils.sysconfig, os, sys
#from distutils.core import setup, Extension
from setuptools import setup, Extension
from Cython.Distutils import build_ext

if not os.environ.has_key('SAGE_ROOT'):
    print "    ERROR: The environment variable SAGE_ROOT must be defined."
    sys.exit(1)
else:
    SAGE_ROOT  = os.environ['SAGE_ROOT']
    SAGE_LOCAL = os.environ['SAGE_LOCAL']
    SAGE_SRC = os.environ['SAGE_SRC']

extra_link_args =  ['-L' + SAGE_LOCAL + '/lib']
extra_compile_args = ['-w', '-O2']

ext_modules = [Extension('dirichlet_conrey', sources=['dirichlet_conrey.pyx', ],
                     library_dirs=[SAGE_LOCAL + '/lib/'],
                     include_dirs=[SAGE_SRC + '/sage/ext'],
                     extra_compile_args = extra_compile_args,
                     extra_link_args = extra_link_args)]
                     

include_dirs = [SAGE_LOCAL + "/include/csage/",
                SAGE_LOCAL + "/include/",
                SAGE_LOCAL + "/include/python2.7",
                SAGE_LOCAL + "/lib/python2.7/site-packages/numpy/core/include/",
                SAGE_SRC + "/sage/ext/",
                SAGE_SRC,
                ]

setup(name='DirichletConrey',
      version='0.111',
      description='desc',
      author='J. W. Bober',
      author_email='jwbober@gmail.com',
      url='http://github.com/jwbober/conrey-dirichlet-characters',
      ext_modules = ext_modules,
      include_dirs = include_dirs,
      cmdclass = {'build_ext' : build_ext})

    

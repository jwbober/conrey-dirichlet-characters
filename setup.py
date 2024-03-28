# simple setup.py, adapted from the one that sage
# automatically produces with 'load diricihet_conrey.pyx"
# and the cython documentation

# Build using 'python setup.py'
import os, sys
#from distutils.core import setup, Extension
from setuptools import setup, Extension
from Cython.Distutils import build_ext
from numpy import get_include
from distutils.sysconfig import get_python_inc

if not 'SAGE_ROOT' in os.environ:
    print("    ERROR: The environment variable SAGE_ROOT must be defined.")
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

SAGE_PYTHON_VERSION = os.getenv('SAGE_PYTHON_VERSION',3)

for e in ext_modules:
    e.cython_directives = {'language_level': "{}".format(SAGE_PYTHON_VERSION)}

include_dirs = [
        SAGE_LOCAL + "/include/csage/",
        SAGE_LOCAL + "/include/",
        get_python_inc(), # python include dir, per distutils
        get_include(),    # numpy include dir, per numpy
        SAGE_SRC + "/sage/ext/",
        SAGE_SRC,
]


setup(name='DirichletConrey',
      version='0.112',
      description='desc',
      author='J. W. Bober',
      author_email='jwbober@gmail.com',
      url='http://github.com/jwbober/conrey-dirichlet-characters',
      ext_modules = ext_modules,
      include_dirs = include_dirs,
      cmdclass = {'build_ext' : build_ext})



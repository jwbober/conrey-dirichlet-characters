# simple setup.py, adapted from the one that sage
# automatically produces with 'load diricihet_conrey.pyx"
# and the cython documentation

# Build using 'python setup.py'
import distutils.sysconfig, os, sys
from distutils.core import setup, Extension
from Cython.Distutils import build_ext

if not os.environ.has_key('SAGE_ROOT'):
    print "    ERROR: The environment variable SAGE_ROOT must be defined."
    sys.exit(1)
else:
    SAGE_ROOT  = os.environ['SAGE_ROOT']
    SAGE_LOCAL = SAGE_ROOT + '/local/'

extra_link_args =  ['-L' + SAGE_LOCAL + '/lib']
extra_compile_args = ['-w', '-O2']

ext_modules = [Extension('dirichlet_conrey', sources=['dirichlet_conrey.pyx', ],
                     library_dirs=[SAGE_LOCAL + '/lib/'],
                     include_dirs=[SAGE_ROOT + '/devel/sage/sage/ext'],
                     libraries=['csage'],
                     extra_compile_args = extra_compile_args,
                     extra_link_args = extra_link_args)]
                     
setup(ext_modules = ext_modules,
      include_dirs = ['/home/bober/sage-5.0.beta5/local/include/csage/', '/home/bober/sage-5.0.beta5/local/include/', '/home/bober/sage-5.0.beta5/local/include/python2.7/', '/home/bober/sage-5.0.beta5/local/lib/python2.7/site-packages/numpy/core/include', '/home/bober/sage-5.0.beta5/devel/sage/sage/ext/', '/home/bober/sage-5.0.beta5/devel/sage/', '/home/bober/sage-5.0.beta5/devel/sage/sage/gsl/', '.'],
      cmdclass = {'build_ext' : build_ext})
    

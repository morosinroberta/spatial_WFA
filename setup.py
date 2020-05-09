from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import os
import numpy
from distutils import sysconfig
#import numpy.distutils.intelccompiler
import numpy.distutils.ccompiler

os.environ["CC"] = "g++"
os.environ["CXX"] = "g++"

from distutils import sysconfig
sysconfig.get_config_vars()['CFLAGS'] = ''
sysconfig.get_config_vars()['OPT'] = ''
sysconfig.get_config_vars()['PY_CFLAGS'] = ''
sysconfig.get_config_vars()['PY_CORE_CFLAGS'] = ''
sysconfig.get_config_vars()['CC'] = 'gcc'
sysconfig.get_config_vars()['CXX'] = 'g++'
sysconfig.get_config_vars()['BASECFLAGS'] = ''
sysconfig.get_config_vars()['CCSHARED'] = ''
sysconfig.get_config_vars()['LDSHARED'] = 'g++ -shared'
sysconfig.get_config_vars()['CPP'] = 'g++'
sysconfig.get_config_vars()['CPPFLAGS'] = ''
sysconfig.get_config_vars()['BLDSHARED'] = ''
sysconfig.get_config_vars()['CONFIGURE_LDFLAGS'] = ''
sysconfig.get_config_vars()['LDFLAGS'] = ''
sysconfig.get_config_vars()['PY_LDFLAGS'] = ''



comp_flags=['-Ofast','-std=c++17','-march=native','-fPIC','-fopt-info-vec']#,
root_dir = '/usr/'
extension = Extension("spatWFA",
                      sources=["WFA.pyx"], 
                      include_dirs=["./", root_dir+"/include/",numpy.get_include()],
                      language="c++",
                      extra_compile_args=comp_flags,
                      extra_link_args=['-shared'],
                      library_dirs=[root_dir+'/lib/','./'],
                      libraries=[])

extension.cython_directives = {'language_level': "3"}

setup(
    name = 'WFA Spat',
    version = '1.0',
    author = 'Jaime de la Cruz Rodriguez & Roberta Morosin (ISP-SU 2020)',
    ext_modules=[extension],
    cmdclass = {'build_ext': build_ext}
)


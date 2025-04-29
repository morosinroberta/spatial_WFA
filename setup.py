from setuptools import setup
from setuptools.extension import Extension
from Cython.Distutils import build_ext
import os
import numpy
import platform as plt
import sys
import pathlib

p = pathlib.Path(sys.executable)
root_dir = str(pathlib.Path(*p.parts[0:-2]))


if(plt.system() == 'Darwin'):
    root_dir = "/opt/local/" # change to the root folder where eigen3 is installed.
    CC = 'clang'
    CXX= 'clang++'
    link_opts = ["-stdlib=libc++","-bundle","-undefined","dynamic_lookup"]
    comp_flags = ["-mcpu=native","-mtune=native"]
else:
    root_dir = '/usr/'
    CC = 'gcc'
    CXX= 'g++'
    link_opts = ["-shared"]
    comp_flags = ["-march=native"]



comp_flags=['-O3', '-flto','-g0','-fstrict-aliasing',\
            '-std=c++14','-fPIC','-fopenmp', '-I./src', "-DNPY_NO_DEPRECATED_API", '-DNDEBUG', \
            '-pedantic', '-Wall']

    
os.environ["CC"] = CC
os.environ["CXX"] = CXX



extension = Extension("spatWFA",
                      sources=["WFA.pyx"], 
                      include_dirs=["./", root_dir+"/include/", numpy.get_include()],
                      language="c++",
                      extra_compile_args=comp_flags,
                      extra_link_args=comp_flags+link_opts,
                      library_dirs=[root_dir+'/lib/','./'],
                      libraries=[])

extension.cython_directives = {'language_level': "3"}

setup(
    name = 'WFA Spat',
    version = '1.0',
    author = 'R. Morosin, J. de la Cruz Rodriguez, G. Vissers & R. Yadav (ISP-SU 2020)',
    ext_modules=[extension],
    cmdclass = {'build_ext': build_ext}
)


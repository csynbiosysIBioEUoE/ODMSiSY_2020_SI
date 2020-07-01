# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 15:05:15 2020

@author: s1778490
"""

from distutils.core import setup, Extension
from Cython.Build import cythonize
#import numpy

setup(
    ext_modules = cythonize("Model1_Cython.pyx")
)

setup(
    ext_modules = cythonize("Model2_Cython.pyx")
)

setup(
    ext_modules = cythonize("Model3_Cython.pyx")
)

setup(
    ext_modules = cythonize("SolveALLCy.pyx")
)

setup(
    ext_modules = cythonize("SolveALLCyGen.pyx")
)

#setup(
#    ext_modules = cythonize("BhattacharyyaDistance2.pyx")
#)

setup(
    ext_modules = cythonize("Utility2S.pyx")
)


setup(
    ext_modules = cythonize("fitGP.pyx")
)

setup(
    ext_modules = cythonize("gaussianprocess.pyx")
)


setup(
    ext_modules = cythonize("stats.pyx")
)




#setup(
#    ext_modules = cythonize("Model1_CythonTest3.pyx"),
#    include_dirs=[numpy.get_include()]
#)



################################# RUN: python setup.py build_ext --inplace
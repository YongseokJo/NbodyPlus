"""
Build script for the potential_energy C++ extension.

Usage:
    python setup_potential.py build_ext --inplace

Requirements:
    pip install pybind11
"""

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys
import os

class get_pybind_include:
    """Helper to defer pybind11 include path until build time."""
    def __str__(self):
        import pybind11
        return pybind11.get_include()


ext_modules = [
    Extension(
        'potential_energy',
        ['potential_energy.cpp'],
        include_dirs=[get_pybind_include()],
        language='c++',
        extra_compile_args=['-O3', '-std=c++17', '-fopenmp', '-fPIC'],
        extra_link_args=['-fopenmp'],
    ),
]


class BuildExt(build_ext):
    """Custom build extension to add compiler flags."""
    c_opts = {
        'unix': ['-O3', '-std=c++17', '-fopenmp'],
    }
    l_opts = {
        'unix': ['-fopenmp'],
    }

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        link_opts = self.l_opts.get(ct, [])

        for ext in self.extensions:
            ext.extra_compile_args = opts
            ext.extra_link_args = link_opts
        build_ext.build_extensions(self)


setup(
    name='potential_energy',
    version='1.0.0',
    author='ABYSS',
    description='High-performance potential energy computation',
    ext_modules=ext_modules,
    cmdclass={'build_ext': BuildExt},
    zip_safe=False,
)

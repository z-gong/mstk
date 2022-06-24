import numpy
from setuptools import setup, Extension

include_dirs = [numpy.get_include()]

libdcd = Extension("mstools.trajectory.libdcd",
                   ['mstools/trajectory/libdcd/libdcd.pyx'],
                   include_dirs=include_dirs + ['mstools/trajectory/libdcd/include'])

libxdr = Extension("mstools.trajectory.libxdr",
                   ['mstools/trajectory/libxdr/libxdr.pyx',
                    'mstools/trajectory/libxdr/src/xdrfile.c',
                    'mstools/trajectory/libxdr/src/xdrfile_xtc.c',
                    'mstools/trajectory/libxdr/src/xdrfile_trr.c',
                    'mstools/trajectory/libxdr/src/trr_seek.c',
                    'mstools/trajectory/libxdr/src/xtc_seek.c', ],
                   include_dirs=include_dirs + ['mstools/trajectory/libxdr/include']
                   )

extensions = [libdcd, libxdr]

setup(
    name="mstools",
    py_modules=["mstools"],
    version="0.1.0",
    install_requires=[],
    include_package_data=True,
    ext_modules=extensions,
    zip_safe=False
)

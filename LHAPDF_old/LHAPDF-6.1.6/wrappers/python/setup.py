#! /usr/bin/env python

import os
from distutils.core import setup
from glob import glob
from distutils.extension import Extension

incdir_src = os.path.abspath("../../include")
incdir_build = os.path.abspath("../../include")
libdir = os.path.abspath("../../src/.libs")

ext = Extension("lhapdf",
                ["lhapdf.cpp"],
                include_dirs = [incdir_src, incdir_build],
                extra_compile_args= " -I/home/home4/institut_1b/dmeuser/master_code/LHAPDF_old/LHAPDF-6.1.6/../local/include".split(),
                library_dirs = [libdir],
                language = "C++",
                libraries = ["stdc++", "LHAPDF"])

setup(name = "LHAPDF",
      version = "6.1.6",
      ext_modules = [ext])

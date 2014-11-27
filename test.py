#!/usr/bin/python

import numpy as np
import cffi
import libcloudphxx as libcl

# CFFI stuff
ffi = cffi.FFI()
lib = ffi.dlopen('libdales4.so')

# C functions
ffi.cdef("void save_ptr(char*,void*);")

# Fortran functions
ffi.cdef("void main();")

# ...
lib.main()

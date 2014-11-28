#!/usr/bin/python

import numpy, os, shutil, subprocess, cffi, libcloudphxx

# CFFI stuff
ffi = cffi.FFI()
clib = ffi.dlopen('ptrutil.so')
flib = ffi.dlopen('libdales4.so')

# C functions
ffi.cdef("void save_ptr(char*,void*);")

# Fortran functions
ffi.cdef("void main(int, char**);")


# executing DALES with BOMEX set-up
bomexdir = "./dales/cases/bomex/"
testdir = "./test/"
argfile = "namoptions.001"

os.mkdir(testdir)
for f in ("lscale.inp.001", "prof.inp.001"):
  os.symlink('../' + bomexdir + f, testdir + f)
shutil.copy(bomexdir + argfile, testdir)

subprocess.call(['sed', '-i', '', '-e', 's/runtime    =  28800/runtime    =  100/', testdir + argfile])
subprocess.call(['sed', '-i', '', '-e', 's/ladaptive  = .true./ladaptive  = .false./', testdir + argfile])

@ffi.callback("void()")
def micro_step():
  print "BQQ"

os.chdir(testdir)
clib.save_ptr("/tmp/micro_step.ptr", micro_step)
flib.main(2, [ ffi.new("char[]", ""), ffi.new("char[]", argfile) ])

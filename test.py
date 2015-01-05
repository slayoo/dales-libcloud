#!/usr/bin/python

import numpy, os, shutil, subprocess, cffi, libcloudphxx, traceback
from params import params

ptrfname = "/tmp/micro_step-" + str(os.getuid()) + "-" + str(os.getpid()) + ".ptr"

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

subprocess.call(['sed', '-i', '-e', 's/runtime    =  28800/runtime    =  100/', testdir + argfile])
subprocess.call(['sed', '-i', '-e', 's/ladaptive  = .true./ladaptive  = .false./', testdir + argfile])

def ptr2np(ptr, size_1, size_2 = 1, size_3 = 1):
  return numpy.frombuffer(
    ffi.buffer(ptr, size_1*size_2*size_3*numpy.dtype(numpy.float64).itemsize),
    dtype=numpy.float64
  ).reshape(size_1, size_2, size_3, order='F').squeeze()
 
# global variables - to be initialised in the first timestep
prtcls = False
first_timestep = True
arrays = {}
dx, dy, dz, dt = 0,0,0,0


@ffi.callback("bool( double*, int     , double*, int,   int,   int,   double*, int,   int,   int,   double*, int,   int,   int,   double*, int,    int,    int,    double*, int,     int,     int    )")
def micro_step(      rhobf  , s1_rhobf, u0,      s1_u0, s2_u0, s3_u0, v0,      s1_v0, s2_v0, s3_v0, w0,      s1_w0, s2_w0, s3_w0, qt0,     s1_qt0, s2_qt0, s3_qt0, thl0,    s1_thl0, s2_thl0, s3_thl0):
  try:
    global prtcls, first_timestep, arrays, dx, dy, dz, dt

    # exposing DALES data through numpy (no copying)
    rhobf = ptr2np(rhobf, s1_rhobf)
    u0    = ptr2np(u0,    s1_u0,   s2_u0,   s3_u0  )
    v0    = ptr2np(v0,    s1_v0,   s2_v0,   s3_v0  )
    w0    = ptr2np(w0,    s1_w0,   s2_w0,   s3_w0  )
    qt0   = ptr2np(qt0,   s1_qt0,  s2_qt0,  s3_qt0 )
    thl0  = ptr2np(thl0,  s1_thl0, s2_thl0, s3_thl0)

    if first_timestep:

      # first, removing the no-longer-needed pointer file
      os.unlink(ptrfname)

      # sanity checks 
      assert u0.shape[0] == v0.shape[0] == w0.shape[0] == qt0.shape[0] == thl0.shape[0]
      assert u0.shape[1] == v0.shape[1] == w0.shape[1] == qt0.shape[1] == thl0.shape[1]
      assert u0.shape[2] == v0.shape[2] == w0.shape[2] == qt0.shape[2] == thl0.shape[2]
      assert rhobf.shape[0] == qt0.shape[2]

      nx, ny, nz = qt0.shape[0]-2, qt0.shape[1]-2, qt0.shape[2]-1

      # TODO! pass dx, dy, dz, dt
      dx, dy, dz = 1, 1, 1
      dt = 1

      # allocating arrays
      arrays["rhod"]    = numpy.empty((nz,))
      arrays["th_d"]    = numpy.empty((nx,   ny,   nz  ))
      arrays["rv"]      = numpy.empty((nx,   ny,   nz  ))
      arrays["rhod_Cx"] = numpy.empty((nx+1, ny,   nz  ))
      arrays["rhod_Cy"] = numpy.empty((nx,   ny+1, nz  ))
      arrays["rhod_Cz"] = numpy.empty((nx,   ny,   nz+1))

      # initialising libcloudph++
      params["opts_init"].dt = dt
      params["opts_init"].nx, params["opts_init"].ny, params["opts_init"].nz = nx, ny, nz
      params["opts_init"].x1, params["opts_init"].y1, params["opts_init"].z1 = dx*nx, dy*ny, dz*nz #TODO: double check
      prtcls = libcloudphxx.lgrngn.factory(params["backend"], params["opts_init"])

    # converting data from DALES for use with the library
    # - DALES has an unused top level
    # - DALES has unit-length halo in x and y
    arrays["rhod"][    :] = rhobf[            0:-1] # TODO: rho  -> rho_d conversion
    arrays["th_d"][:,:,:] = thl0[ 1:-1, 1:-1, 0:-1] # TODO: th_l -> th_d  conversion
    arrays["rv"  ][:,:,:] = qt0[  1:-1, 1:-1, 0:-1] # TODO: qt   -> rv    conversion

    # this assumes DALES u0, v0 and w0 have periodic condition in the halo slabs
    arrays["rhod_Cx"][:, :, :] = u0[1:,   1:-1, 0:-1] * dt / dx * arrays["rhod"].reshape(1,1,arrays["rhod"].shape[0])
    arrays["rhod_Cy"][:, :, :] = v0[1:-1, 1:,   0:-1] * dt / dy * arrays["rhod"].reshape(1,1,arrays["rhod"].shape[0])
    arrays["rhod_Cz"][:, :, :] = w0[1:-1, 1:-1, 0:  ] * dt / dz #* arrays["rhod"].reshape(1,1,arrays["rhod"].shape[0]) #TODO <- rho at half levels needed

    # initialising the particles 
    if first_timestep:
      print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
      print "+ Python: initialising particle properties (call to C++) +"
      print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
      prtcls.init(arrays["th_d"], arrays["rv"], arrays["rhod"])
      #diagnostics(prtcls, 1, size_x, size_z) # writing down state at t=0

    # the timestepping
    print "++++++++++++++++++++++++++++++++++++++"
    print "+ Python: timestepping (call to C++) +"
    print "++++++++++++++++++++++++++++++++++++++"
    prtcls.step_sync(params["opts"], arrays["th_d"], arrays["rv"], arrays["rhod"]) #TODO: Courants
    prtcls.step_async(params["opts"])

    first_timestep = False
  except:
    traceback.print_exc()
    return False
  else:
    return True

# storing the pointer
clib.save_ptr(ptrfname, micro_step)

# running DALES
os.chdir(testdir)
flib.main(2, [ ffi.new("char[]", ""), ffi.new("char[]", argfile) ])

#!/usr/bin/python

import numpy, os, shutil, subprocess, cffi, libcloudphxx, traceback, math
from params import params
from diag import diag

# skip NaN checks
numpy.seterr(all='ignore')

ptrfname = "/tmp/micro_step-" + str(os.getuid()) + "-" + str(os.getpid()) + ".ptr"

# CFFI stuff
ffi = cffi.FFI()
clib = ffi.dlopen('ptrutil.so')
flib = ffi.dlopen('libdales4.so')

# C functions
ffi.cdef("void save_ptr(char*,void*);")

# executing DALES with BOMEX set-up
bomexdir = "./dales/cases/bomex/"
testdir = "./test/"
argfile = "namoptions"

os.mkdir(testdir)
for f in ("lscale.inp.001", "prof.inp.001"):
  os.symlink('../' + bomexdir + f, testdir + f)
shutil.copy(bomexdir + argfile + '.001', testdir + argfile) # stripping .001 as this will be the default handled when compiled with ifort which generates no-argument main()

subprocess.call(['sed', '-i', '-e', 's/runtime    =  28800/runtime    =  ' + str(params["runtime"]) + '/', testdir + argfile]) 
subprocess.call(['sed', '-i', '-e', 's/ladaptive  = .true./ladaptive  = .false.\\\ndtmax = 2/', testdir + argfile])
subprocess.call(['sed', '-i', '-e', 's/lfielddump  = .false./lfielddump  = .true./', testdir + argfile])

def ptr2np(ptr, size_1, size_2 = 1, size_3 = 1):
  return numpy.frombuffer(
    ffi.buffer(ptr, size_1*size_2*size_3*numpy.dtype(numpy.float64).itemsize),
    dtype=numpy.float64
  ).reshape(size_1, size_2, size_3, order='F').squeeze()
 
# global variables - to be initialised in the first timestep
prtcls = False
first_timestep = True
arrays = {}

# constants from DALES #TODO: allow overriding them in libcloud? or override them in DALES?
L = 2.5e6
c_pd = 1004.
p_1000 = 100000.
R_d = 287.0
R_v = 461.5
eps = R_d / R_v
rho_w = 1e3

def th_std2dry(th, rv):
  from libcloudphxx.common import R_v, R_d, c_pd
  return th * (1 + rv * R_v / R_d)**(R_d/c_pd)

def q2r(q):
  return q / (1 - q)

def r2q(r):
  return r / (1 + r)

@ffi.callback("bool(double, double, double, double, double*, int,     double*, int,   int,   int,   double*, int,   int,   int,   double*, int,   int,   int,   double*, int,    int,    int,    double*, int,     int,     int    )")
def micro_step(     dt,     dx,     dy,     dz,     exnf,    s1_exnf, u0,      s1_u0, s2_u0, s3_u0, v0,      s1_v0, s2_v0, s3_v0, w0,      s1_w0, s2_w0, s3_w0, qt0,     s1_qt0, s2_qt0, s3_qt0, thl0,    s1_thl0, s2_thl0, s3_thl0):
  try:
    global prtcls, first_timestep, arrays

    # exposing DALES data through numpy (no copying)
    # - DALES has an unused top level
    # - DALES has unit-length halo in x and y

    exnf  = ptr2np(exnf,  s1_exnf)
    u0    = ptr2np(u0,    s1_u0,   s2_u0,   s3_u0  )
    v0    = ptr2np(v0,    s1_v0,   s2_v0,   s3_v0  )
    w0    = ptr2np(w0,    s1_w0,   s2_w0,   s3_w0  )
    qt0   = ptr2np(qt0,   s1_qt0,  s2_qt0,  s3_qt0 ) 
    thl0  = ptr2np(thl0,  s1_thl0, s2_thl0, s3_thl0)

    if first_timestep:

      # first, removing the no-longer-needed pointer file
      os.unlink(ptrfname)

      # sanity checks for array shapes
      assert u0.shape[0] == v0.shape[0] == w0.shape[0] == qt0.shape[0] == thl0.shape[0]
      assert u0.shape[1] == v0.shape[1] == w0.shape[1] == qt0.shape[1] == thl0.shape[1]
      assert u0.shape[2] == v0.shape[2] == w0.shape[2] == qt0.shape[2] == thl0.shape[2]
      assert exnf.shape[0] == qt0.shape[2]

      # sub-second timestep to cope with condensation
      assert dt / params["opts_init"].sstp_cond < 1

      nx, ny, nz = qt0.shape[0]-2, qt0.shape[1]-2, qt0.shape[2]-1

      # allocating arrays
      arrays["rhod"]    = numpy.empty((nz,))
      arrays["th_d"]    = numpy.empty((nx,   ny,   nz  ))
      arrays["rv"]      = numpy.empty((nx,   ny,   nz  ))
      arrays["Cx"]      = numpy.empty((nx+1, ny,   nz  ))
      arrays["Cy"]      = numpy.empty((nx,   ny+1, nz  ))
      arrays["Cz"]      = numpy.empty((nx,   ny,   nz+1))

      # initialising libcloudph++
      params["opts_init"].dt = dt
      params["opts_init"].dx, params["opts_init"].dy, params["opts_init"].dz = dx, dy, dz
      params["opts_init"].nx, params["opts_init"].ny, params["opts_init"].nz = nx, ny, nz
      params["opts_init"].x1, params["opts_init"].y1, params["opts_init"].z1 = dx*nx, dy*ny, dz*nz

      try:
        print("Trying with CUDA backend..."),
        prtcls = libcloudphxx.lgrngn.factory(libcloudphxx.lgrngn.backend_t.CUDA, params["opts_init"])
        print (" OK!")
      except:
        print (" KO!")
        try:
          print("Trying with OpenMP backend..."),
          prtcls = libcloudphxx.lgrngn.factory(libcloudphxx.lgrngn.backend_t.OpenMP, params["opts_init"])
          print (" OK!")
        except:
          print (" KO!")
          print("Trying with serial backend..."),
          prtcls = libcloudphxx.lgrngn.factory(libcloudphxx.lgrngn.backend_t.serial, params["opts_init"])
          print (" OK!")

      # calculating rho_d profile (constant-in-time)
      # assuming there is no liquid water at t==0
      rv = q2r(qt0[ 1:-1, 1:-1, 0:-1].mean(axis=0).mean(axis=0))
      pd = exnf[0:-1]**(c_pd / R_d) * p_1000 * eps / (rv + eps)
      T  = exnf[0:-1] * thl0[1:-1, 1:-1, 0:-1].mean(axis=0).mean(axis=0)
      arrays["rhod"][:] = pd / R_d / T

    else:
      # we assume adaptivity is turned off
      assert dt == params["opts_init"].dt

      #TODO: assert for invariant rho_d

      #TODO: assert for no streaching (constant dz)

    prtcls.diag_all() # TODO: DALES q_t does not include q_r (nor q_a)
    prtcls.diag_wet_mom(3)
    rc = numpy.frombuffer(prtcls.outbuf()).reshape(arrays["rv"].shape) * 4./3. * math.pi * rho_w
    
    arrays["rv"  ][:,:,:] = q2r(qt0[1:-1, 1:-1, 0:-1]) - rc
    arrays["th_d"] = th_std2dry(thl0[1:-1, 1:-1, 0:-1] + r2q(rc) / exnf[0:-1] * L / c_pd, arrays["rv"])

    # assert if temperature is properly recovered from rho_d and th_d?
    # (also assuming thl == theta hence only at t == 0)
    if __debug__ and first_timestep:
      mxdiff = 0
      for i in range(arrays["th_d"].shape[0]):
	for j in range(arrays["th_d"].shape[1]):
	  for k in range(arrays["th_d"].shape[2]):
	    mxdiff = max(
              mxdiff, 
              exnf[k] * thl0[i+1, j+1, k] - libcloudphxx.common.T(arrays["th_d"][i,j,k], arrays["rhod"][k])
            )
      assert mxdiff < 0.05

    # this assumes DALES u0, v0 and w0 have periodic condition in the halo slabs
    arrays["Cx"][:, :, :] = u0[1:,   1:-1, 0:-1] * dt / dx
    arrays["Cy"][:, :, :] = v0[1:-1, 1:,   0:-1] * dt / dy
    arrays["Cz"][:, :, :] = w0[1:-1, 1:-1, 0:  ] * dt / dz 

    if first_timestep:
      # initialising particle properties
      prtcls.init(arrays["th_d"], arrays["rv"], arrays["rhod"])
      # writing down state at t=0
      diag(prtcls)

    # the timestepping
    prtcls.step_sync(
      params["opts"], 
      arrays["th_d"], 
      arrays["rv"], 
      arrays["Cx"],
      arrays["Cy"],
      arrays["Cz"]
    ) 
    prtcls.step_async(params["opts"]) #TODO: handle the async logic

    diag(prtcls) 

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
ok = False
for main in [
  'MAIN__', # ifort
  'main'    # gfortran
]:
  ffi.cdef("void "+main+"(int, char**);") 
  try:
    getattr(flib, main)(2, [ ffi.new("char[]", ""), ffi.new("char[]", argfile) ])
  except:
    pass
  else:
    ok = True
if not ok: raise Exception("main call failed")

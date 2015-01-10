from scipy.io import netcdf
import numpy, pdb

t = 0
f = None
dtav = 60 # TODO: read it from file or pass it?
shape = None
outfreq = None

def diagnostics(prtcls):
  global t,f,shape,outfreq

  # file definition
  if t == 0:
    shape = (prtcls.opts_init.nx, prtcls.opts_init.ny, prtcls.opts_init.nz)
    outfreq = int(dtav / prtcls.opts_init.dt)

    f = netcdf.netcdf_file('libcloud.nc', 'w')

    f.createDimension('time', 0) 
    f.createDimension('xt', shape[0]) 
    f.createDimension('yt', shape[1]) 
    f.createDimension('zt', shape[2]) 

    f.createVariable('time', numpy.float32, ('time',))
    f.variables['time'].longname = 'Time'
    f.variables['time'].units = 's'

    f.createVariable('sd_conc', numpy.float32, ('time','zt','yt','xt',))
    f.variables['sd_conc'].longname = 'super-droplet concentration'
    f.variables['sd_conc'].units = '1/dx/dy/dz'

  # filling the file with data
  if t % outfreq == 0:
    # Time
    f.variables['time'][t / outfreq] = t * prtcls.opts_init.dt

    # sd_conc
    prtcls.diag_sd_conc()
    f.variables['sd_conc'][t / outfreq,:,:,:] = numpy.swapaxes(numpy.frombuffer(prtcls.outbuf()).reshape(shape),0,2)

  print t
  t = t+1

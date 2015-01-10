from scipy.io import netcdf
import numpy, math

t = 0
f = None
dtav = 60 # TODO: read it from file or pass it?
shape = None
outfreq = None

# aerosol / cloud / rain thresholds
th_ac = 0.5e-6
th_cr = 25e-6

def getbuf(prtcls):
  return numpy.swapaxes(
    numpy.frombuffer(
      prtcls.outbuf()
    ).reshape((
      prtcls.opts_init.nx,
      prtcls.opts_init.ny,
      prtcls.opts_init.nz
    )),
    0,2
  )

def newvar(f, name, longname, unit):
  f.createVariable(name, numpy.float32, ('time','zt','yt','xt'))
  f.variables[name].longname = longname
  f.variables[name].units = unit

def diagnostics(prtcls):
  global t,f,shape,outfreq

  # file definition
  if t == 0:
    outfreq = int(dtav / prtcls.opts_init.dt)

    f = netcdf.netcdf_file('libcloud.nc', 'w')

    f.createDimension('time', 0) 
    f.createDimension('xt', prtcls.opts_init.nx) 
    f.createDimension('yt', prtcls.opts_init.ny) 
    f.createDimension('zt', prtcls.opts_init.nz) 

    f.createVariable('time', numpy.float32, ('time',))
    f.variables['time'].longname = 'Time'
    f.variables['time'].units = 's'

    newvar(f, 'sd_conc', 'super-droplet concentration',                       '1/dx/dy/dz') # TODO: int?
    newvar(f, 'na',      'aerosol concentration (per mass of dry air)',       '1/kg'      )
    newvar(f, 'nc',      'cloud droplet concentration (per mass of dry air)', '1/kg'      )
    newvar(f, 'nr',      'rain drop concentration (per mass of dry air)',     '1/kg'      )
    newvar(f, 'qc',      'cloud water mixing ratio',                          'kg/kg'     )
    newvar(f, 'r_eff',   'cloud droplet effective radius',                    'm'         )
    newvar(f, 'qr',      'rain water mixing ratio',                           'kg/kg'     )

  # filling the file with data
  if t % outfreq == 0:
    rec = t / outfreq

    # Time
    f.variables['time'][rec] = t * prtcls.opts_init.dt

    # sd_conc
    prtcls.diag_sd_conc()
    f.variables['sd_conc'][rec,:,:,:] = getbuf(prtcls)

    # na
    prtcls.diag_wet_rng(0, th_ac)
    prtcls.diag_wet_mom(0)
    f.variables['na'][rec,:,:,:] = getbuf(prtcls)

    # ql & nc & r_eff
    prtcls.diag_wet_rng(th_ac, th_cr)
    prtcls.diag_wet_mom(0)
    f.variables['nc'][rec,:,:,:] = getbuf(prtcls)
    prtcls.diag_wet_mom(3)
    f.variables['qc'][rec,:,:,:] = getbuf(prtcls)
    prtcls.diag_wet_mom(2)
    f.variables['r_eff'][rec,:,:,:] = f.variables['qc'][rec,:,:,:] / getbuf(prtcls)
    f.variables['qc'][rec,:,:,:] *= 4./3 * math.pi # * rho_w

    # qr & nr
    prtcls.diag_wet_rng(th_cr, 1)
    prtcls.diag_wet_mom(0)
    f.variables['nr'][rec,:,:,:] = getbuf(prtcls)
    prtcls.diag_wet_mom(3)
    f.variables['qr'][rec,:,:,:] = getbuf(prtcls) * 4./3 * math.pi # * rho_w

  # incrementing timestep counter
  t = t+1

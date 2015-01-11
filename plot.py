import numpy, Gnuplot
from scipy.io import netcdf

def path(n, t, var, dz):
  return numpy.sum(numpy.squeeze(n.variables[var][t,:,:,:]), axis=0) * dz

nd = netcdf.netcdf_file('test/fielddump.000.001.nc', 'r')
nl = netcdf.netcdf_file('test/libcloud.nc', 'r')

x = nd.variables['xt'][:]
y = nd.variables['yt'][:]
z = nd.variables['zt'][:]

dx = x[1] - x[0]
dy = y[1] - y[0]
dz = z[1] - z[0]

g = Gnuplot.Gnuplot(persist=1)
g('set xrange [0:' + str(x[-1]+dx/2) + ']')
g('set yrange [0:' + str(y[-1]+dy/2) + ']')
g('set xlabel "Y"')
g('set ylabel "X"')
g('set cbrange [100:50000]')
g('set logscale cb')
g('set view map')
g('set multiplot layout 2,1')

for t in range(nd.variables['time'].shape[0]):
  #TODO: it should work with binary=1!
  dd = Gnuplot.GridData(path(nd, t, 'ql', dz), y, x, with_='image', binary=0) 
  dl = Gnuplot.GridData(path(nl, t, 'sd_conc', dz), y, x, with_='image', binary=0)
  g.splot(dd)
  g.splot(dl)

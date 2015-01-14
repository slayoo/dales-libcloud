# This Python file uses the following encoding: utf-8

import numpy, Gnuplot
from scipy.io import netcdf

rng_lwp = (2e-5,   1      )
rng_ref = (1.5e-6, 21e-6  )
rng_lwc = (.05e-3, .75e-3 )
rng_nc  = (10e6,   160e6  )


def path(n, t, var, mlt, dz): #TODO: multiply by rho_d
  path = numpy.sum(numpy.squeeze(n.variables[var][t,:,:,:]), axis=0) * dz * mlt
  return path 

def prof(n, t, var, mlt):
  vals = numpy.nan_to_num(numpy.squeeze(n.variables[var][t,:,:,:])) * mlt
  if var == "qc" or var == "ql": 
    vals[vals < rng_lwc[0]] = numpy.nan
  if var == "nc": 
    vals[numpy.squeeze(numpy.nan_to_num(n.variables["r_eff"][t,:,:,:]))<rng_ref[0]]= numpy.nan
  if var == "r_eff": vals[vals < mlt*rng_ref[0]] = numpy.nan
  return numpy.nanmean(numpy.nanmean(vals, axis=1), axis=1)

def hist2d(n, z, var, rng, bns, mlt):
  vals = numpy.squeeze(n.variables[var][t,:,:,:]).copy()
  vals[numpy.isnan(vals)]=-1
  if var != "r_eff" and "r_eff" in n.variables: 
    vals[numpy.squeeze(numpy.nan_to_num(n.variables["r_eff"][t,:,:,:]))<rng_ref[0]]=-2
  bins, xedges, yedges = numpy.histogram2d(
    vals.ravel(),
    numpy.repeat(z, vals.shape[1]*vals.shape[2]),
    range=rng, bins=bns
  )
  xedges = xedges[0:-1] + (xedges[1] - xedges[0]) / 2
  xedges = xedges / mlt
  yedges = yedges[0:-1] + (yedges[1] - yedges[0]) / 2
  return bins, xedges, yedges


nd = netcdf.netcdf_file('test.skua.3600/fielddump.000.001.nc', 'r')
nl = netcdf.netcdf_file('test.skua.3600/libcloud.nc', 'r')

x = nd.variables['xt'][:] / 1e3
y = nd.variables['yt'][:] / 1e3
z = nd.variables['zt'][:] / 1e3

dx = x[1] - x[0]
dy = y[1] - y[0]
dz = z[1] - z[0]

g = Gnuplot.Gnuplot()# persist=1)

g('set term pdf enhanced rounded size 24cm,16cm')

for t in range(nd.variables['time'].shape[0]):
  #if t != 44: continue

  g('reset')
  g('set xrange [0:' + str(x[-1]+dx/2) + ']')
  g('set yrange [0:' + str(y[-1]+dy/2) + ']')
  g('set xlabel "Y [km]"')
  g('set ylabel "X [km]"')
  g('set size square')
  g('set cbrange [' + str(rng_lwp[0]) + ':' + str(rng_lwp[1]) + ']')
  g('set cbtics offset -1')
  g('set palette defined (0 "white", .25 "blue", 1 "yellow")')
  g('set logscale cb')
  g('set view map')
  g('set grid')

  g('set output "plot/' + str("%03d" % t) + '.pdf"')
  g('set multiplot layout 2,3')

  g('set object 1 rectangle from screen .001,.003 to screen .338,.997 fillcolor rgb"#ffffff" behind')
  g('set object 2 rectangle from screen .341,.003 to screen .999,.997 fillcolor rgb"#ffffff" behind')

  g('set label 1 "DALES (BOMEX, bulk μ-physics, t=' + str("%d" % (nd.variables['time'][t]/60)) + 'm)" at screen .004,.98 left font ",16"')
  g('set label 2 "DALES-piggybacked libcloduph++ Lagrangian/Monte-Carlo μ-physics on a GPU" at screen .346,.98 left font ",16"')

  g('set title "LWP [kg/m^2]"')
  dlwp = path(nd, t, 'ql', 1e-5, dz*1e3)
  g.splot(Gnuplot.GridData(dlwp, y, x, with_='image', binary=0))

  g('unset object 1')
  g('unset object 2')
  g('unset label 1')
  g('unset label 2')

  llwp = path(nl, t, 'qc', 1e3, dz*1e3) #TODO: get rid of 1e3 after correcting diag.py
  g.splot(Gnuplot.GridData(llwp, y, x, with_='image', binary=0))

  #print numpy.corrcoef(llwp.flat, dlwp.flat)

  g('set logscale x')
  g('set logscale y2')
  g('set xrange [' + str(rng_lwp[0]) + ':' + str(rng_lwp[1]) + ']')
  g('set y2range [' + str(rng_lwp[0]) + ':' + str(rng_lwp[1]) + ']')
  g('set xlabel "DALES"')
  g('set y2label "libcloudph++" offset -2')
  g('unset ytics')
  g('set ylabel " "')
  g('set y2tics')
  g('set grid y2tics')
  g('set title "LWP [kg/m^2] (scatter plot)"')
  g.plot(
    Gnuplot.Data(dlwp.ravel(), llwp.ravel(), with_='points lt 4 pt 5 ps .3', axes='x1y2'),
    Gnuplot.Data(rng_lwp, rng_lwp, with_='lines lt rgb "black"', axes='x1y2')
  )


  g('unset logscale')
  g('unset title')
  g('set ytics')
  g('unset y2tics')
  g('unset y2label')

  g('set ylabel "Z [km]"')
  g('set yrange [0:1.8]')
  g('set yrange [.4:1.2]')

  g('set palette defined (0 "white", .2 "white", .55 "orange", .66 "green", 1 "olive")')
  g('set cblabel "freq of occurance, log sclale"')
  g('unset cbtics')
  g('unset cbrange')
  g('set logscale cb')

  g('set key bottom right')
  g('set xlabel "q_c [g/kg]"')
  g('set xrange [' + str(rng_lwc[0] * 1e3) + ':' + str(rng_lwc[1] * 1e3) + ']')
  g('set title "2D histogram (q_c > ' + str(rng_lwc[0]*1e3) + ')"')
  bins, xedges, yedges = hist2d(nd, z, "ql", 
    [
      [rng_lwc[0]*1e5, rng_lwc[1]*1e5],
      [z[0]-dz/2,  z[-1]+dz/2]
    ],
    [21,80],
    1e5/1e3
  )
  g.splot(
    Gnuplot.GridData(bins, xedges, yedges, with_='image', binary=0),
    Gnuplot.Data(prof(nl, t, "qc", 1e3) * 1e3, z, numpy.zeros(z.shape), with_='lines lw 2', title='libcloudph++ mean'),
    Gnuplot.Data(prof(nd, t, "ql", 1e-5) * 1e3, z, numpy.zeros(z.shape), with_='lines lt 6 lw 4', title='DALES mean')
  )

  g('set title "2D histogram (r_{eff} > ' + str(rng_ref[0]*1e6) + ')"')

  g('set xlabel "r_{eff} [μm]"')
  g('set xtics 5')
  g('set xrange [' + str(rng_ref[0]*1e6) + ':' + str(rng_ref[1]*1e6) + ']')


  bins, xedges, yedges = hist2d(nl, z, "r_eff", 
    [
      [rng_ref[0], rng_ref[1]],
      [z[0]-dz/2,  z[-1]+dz/2]
    ],
    [21,80],
    1e-6
  )
  g.splot(
    Gnuplot.GridData(bins, xedges, yedges, with_='image', binary=0),
    Gnuplot.Data(prof(nl, t, "r_eff", 1e6), z, numpy.zeros(z.shape), with_='lines lw 2 ', title='libcloudph++ mean'),
  )

  g('set xlabel "n_c [cm^{-3}]"') #TODO: in fact it is still mg-1
  g('set xrange [' + str(rng_nc[0]/1e6) + ':' + str(rng_nc[1]/1e6) + ']')
  g('set xtics auto')
  bins, xedges, yedges = hist2d(nl, z, "nc", 
    [[rng_nc[0],rng_nc[1]],[z[0]-dz/2, z[-1]+dz/2]],
    [21,80],
    1e6
  )
  g.splot(
    Gnuplot.GridData(bins, xedges, yedges, with_='image', binary=0),
    Gnuplot.Data(prof(nl, t, "nc", 1e-6), z, numpy.zeros(z.shape), with_='lines lw 2', title='libcloudph++ mean'),
  )

  g('unset multiplot')

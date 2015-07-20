# This Python file uses the following encoding: utf-8

import numpy, Gnuplot
from scipy.io import netcdf

rng_lwp = numpy.array([.05,   .5      ])
rng_ref = numpy.array([1.5e-6, 21e-6  ])
rng_lwc = numpy.array([.05e-3, .75e-3 ])
rng_nc  = numpy.array([10e6,   1600e6  ])


def path(n, t, var, mlt, dz, rho): 
  path = numpy.sum(numpy.squeeze(n.variables[var][t,:,:,:]*rho.reshape((rho.shape[0], 1, 1))), axis=0) * dz * mlt
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


root = 'test.skua_n1e3'
nd = netcdf.netcdf_file(root + '/fielddump.000.001.nc', 'r')
nl = netcdf.netcdf_file(root + '/libcloud.nc', 'r')
np = netcdf.netcdf_file(root + '/profiles.001.nc', 'r')

x = nd.variables['xt'][:] / 1e3
y = nd.variables['yt'][:] / 1e3
z = nd.variables['zt'][:] / 1e3

dx = x[1] - x[0]
dy = y[1] - y[0]
dz = z[1] - z[0]

g = Gnuplot.Gnuplot()

g('set term pdf enhanced rounded')

for t in range(nd.variables['time'].shape[0]):
  if t != 45: continue
  rhof = np.variables['dn0'][t/10,:] # TODO: get rid of /10 factor!

  g('reset')
  g('set xrange [0:' + str(x[-1]+dx/2) + ']')
  g('set yrange [0:' + str(y[-1]+dy/2) + ']')
  g('set xlabel "Y [km]"')
  g('set ylabel "X [km]"')
  g('set size square')
  g('set cbrange [' + str(rng_lwp[0]) + ':' + str(rng_lwp[1]) + ']')
  g('set cbtics offset -1')
  g('set palette defined (0 "white", .25 "blue", 1 "yellow")')
  g('set view map')
  g('set grid')


  g('set title "DALES LWP [kg/m^2]"')
  dlwp = path(nd, t, 'ql', 1e-5, dz*1e3, rhof)
  g('set output "plot/' + str("%03d" % t) + '_lwpmap_dales.pdf"')
  g.splot(Gnuplot.GridData(dlwp, y, x, with_='image', binary=0))


  g('set title "libcloudph++ LWP [kg/m^2]"')
  llwp = path(nl, t, 'qc', 1, dz*1e3, rhof)
  g('set output "plot/' + str("%03d" % t) + '_lwpmap_libcloud.pdf"')
  g.splot(Gnuplot.GridData(llwp, y, x, with_='image', binary=0))

  g('set xrange [' + str(rng_lwp[0]) + ':' + str(rng_lwp[1]) + ']')
  g('set y2range [' + str(rng_lwp[0]) + ':' + str(rng_lwp[1]) + ']')
  g('set xlabel "DALES"')
  g('set y2label "libcloudph++" offset -1')
  g('unset ytics')
  g('set ylabel " "')
  g('set y2tics')
  g('set grid y2tics')
  g('set title "LWP [kg/m^2] (scatter plot)"')
  g('set output "plot/' + str("%03d" % t) + '_lwpcorr.pdf"')
  g.plot(
    Gnuplot.Data(dlwp.ravel(), llwp.ravel(), with_='points lt 4 pt 5 ps .3', axes='x1y2'),
    Gnuplot.Data(rng_lwp, rng_lwp, with_='lines lt rgb "black"', axes='x1y2')
  )

  g('set title "q_c [g/kg] (scatter plot)"')
  g('set output "plot/' + str("%03d" % t) + '_qccorr.pdf"')
  g('set xrange [' + str(rng_lwc[0] * 1e3) + ':' + str(rng_lwc[1] * 1e3) + ']')
  g('set y2range [' + str(rng_lwc[0] * 1e3) + ':' + str(rng_lwc[1] * 1e3) + ']')
  dqc = 1e3*1e-5*nd.variables['ql'][t,:,:,:]
  lqc = 1e3*nl.variables['qc'][t,:,:,:]
  g.plot(
    Gnuplot.Data(dqc.ravel(), lqc.ravel(), with_='points lt 4 pt 5 ps .3', axes='x1y2'),
    Gnuplot.Data(1e3*rng_lwc, 1e3*rng_lwc, with_='lines lt rgb "black"', axes='x1y2')
  )

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
  g('set output "plot/' + str("%03d" % t) + '_prof_qc.pdf"')
  g.splot(
    Gnuplot.GridData(bins, xedges, yedges, with_='image', binary=0),
    Gnuplot.Data(prof(nl, t, "qc", 1) * 1e3, z, numpy.zeros(z.shape), with_='lines lw 2', title='libcloudph++ mean'),
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
  g('set output "plot/' + str("%03d" % t) + '_prof_reff.pdf"')
  g.splot(
    Gnuplot.GridData(bins, xedges, yedges, with_='image', binary=0),
    Gnuplot.Data(prof(nl, t, "r_eff", 1e6), z, numpy.zeros(z.shape), with_='lines lw 2 ', title='libcloudph++ mean'),
  )

  g('set xlabel "n_c [mg^{-1}]"') 
  g('set xrange [' + str(rng_nc[0]/1e6) + ':' + str(rng_nc[1]/1e6) + ']')
  g('set xtics 250')
  bins, xedges, yedges = hist2d(nl, z, "nc", 
    [[rng_nc[0],rng_nc[1]],[z[0]-dz/2, z[-1]+dz/2]],
    [21,80],
    1e6
  )
  g('set output "plot/' + str("%03d" % t) + '_prof_nc.pdf"')
  g.splot(
    Gnuplot.GridData(bins, xedges, yedges, with_='image', binary=0),
    Gnuplot.Data(prof(nl, t, "nc", 1e-6), z, numpy.zeros(z.shape), with_='lines lw 2', title='libcloudph++ mean'),
  )

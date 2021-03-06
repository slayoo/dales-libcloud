from libcloudphxx import lgrngn

def lognormal(lnr):
  from math import exp, log, sqrt, pi
  n_tot = 100e6
  meanr = .04e-6
  gstdv = 1.4

  return n_tot * exp(
    -pow((lnr - log(meanr)), 2) / 2 / pow(log(gstdv),2)
  ) / log(gstdv) / sqrt(2*pi);


opts_init = lgrngn.opts_init_t()
opts_init.dry_distros = { .61 : lognormal }
opts_init.sd_conc_mean = .1 # keep this setting low for Travis!
opts_init.sstp_cond = 40 # keep this setting low for Travis!
opts_init.sstp_coal = 0
opts_init.kernel = lgrngn.kernel_t.geometric

opts = lgrngn.opts_t()
opts.cond = True
opts.coal = False
opts.adve = True 
opts.sedi = False
opts.chem = False

# the default DALES imicro=0 setting disables rain

params = {
  "opts"      : opts,
  "opts_init" : opts_init,
  "runtime"   : 60 # DALES override
}


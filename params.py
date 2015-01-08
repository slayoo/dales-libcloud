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
opts_init.sd_conc_mean = 10
opts_init.dry_distros = { .61 : lognormal }

opts = lgrngn.opts_t()
opts.sstp_cond = 100
opts.sstp_coal = 1
opts.cond = True
opts.coal = True
opts.adve = True
opts.sedi = True
opts.chem = False
opts.kernel = lgrngn.kernel_t.geometric

params = {
  "backend"   : lgrngn.backend_t.serial,
  "opts"      : opts,
  "opts_init" : opts_init
}

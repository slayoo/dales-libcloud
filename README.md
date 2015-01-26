dales-libcloudph++ coupling example
===================================

Credits (alph. order):
  - Sylwester Arabas (Py-C++ coupling, DALES-specific coding), 
  - Dorota Jarecka (Py-Fortran coupling), 
  - Harm Jonker (original idea, DALES internals consultancy)

Summary:
  - dales.diff       (CMake: shared lib & hardcoded mpi compiler; fielddump @t=0; modlibcloud call)
  - modlibcloud.f90  (load fprt; call Python passing dt,dx,rhob,exnf,uvw,qt,th; skip rk3 & spinup; error handing)
  - ptrutil.c        (file store/load f-ction pointer - impossible in Fortran?; cubmersone in Python)
  - main.py          (Py->C->F->Py logic; state var conversion (with sanity checks); grid geometry; lib & diag calls; DALES namelist overrides)
  - params.py        (process toggling; aerosol; numerical params - sd_conc, substeps)
  - diag.py          (fielddump-compatible separate nc file for cloud properties; only place where C/F memory layout conv needed additional calc.)
  - .travis.yml      (automated in-the-cloud test on a fresh Ubuntu install; incl. DALES clone, patch & compilation)
  - plot.py          (wip, uses Python-gnuplot)

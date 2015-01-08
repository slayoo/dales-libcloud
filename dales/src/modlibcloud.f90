! autor: Sylweter Arabas
! license: GPL

module modlibcloud

  use modglobal, only : rk3step, rdt, dx, dy, dz, ntimee
  use modfields, only : rhobf, rhobh, rhof, u0, v0, w0, qt0, thl0
  use iso_c_binding, only: c_funptr, c_f_procpointer, c_null_char, c_double

  interface
    ! C function used to load a pointer to Python function
    subroutine load_ptr(fname, ptr) bind(c)
      use iso_c_binding, only: c_funptr, c_char
      character(kind=c_char), dimension(*), intent(in) :: fname
      type(c_funptr), intent(out) :: ptr
    end

    ! Python function used to call C++ library
    function micro_step_py(                &
      rdt, dx, dy, dz,                     &
      rhobf, s1_rhobf,                     &
      rhobh, s1_rhobh,                     &
      rhof,  s1_rhof,                      &
      u0,    s1_u0,   s2_u0,   s3_u0,      &
      v0,    s1_v0,   s2_v0,   s3_v0,      &
      w0,    s1_w0,   s2_w0,   s3_w0,      &
      qt0,   s1_qt0,  s2_qt0,  s3_qt0,     &
      thl0,  s1_thl0, s2_thl0, s3_thl0     &
    ) bind(c) 

      use iso_c_binding, only: c_double, c_int, c_bool

      logical(c_bool) :: micro_step_py

      integer(c_int), intent(in), value :: &
        s1_rhobf,                          &
        s1_rhobh,                          &
        s1_rhof,                           &
        s1_u0,   s2_u0,   s3_u0,           &
        s1_v0,   s2_v0,   s3_v0,           &
        s1_w0,   s2_w0,   s3_w0,           &
        s1_qt0,  s2_qt0,  s3_qt0,          &
        s1_thl0, s2_thl0, s3_thl0

      real(c_double), intent(in), value :: &
        rdt, dx, dy, dz                       

      real(c_double), intent(in) ::        &
        rhobf(s1_rhobf),                   &
        rhobh(s1_rhobh),                   &
        rhof(s1_rhof),                     &
        u0(s1_u0, s2_u0, s3_u0),           &
        v0(s1_v0, s2_v0, s3_v0),           &
        w0(s1_w0, s2_w0, s3_w0),           &
        qt0(s1_qt0, s2_qt0, s3_qt0),       &
        thl0(s1_thl0, s2_thl0, s3_thl0)
    end
  end interface

  type(c_funptr) :: cptr
  procedure(micro_step_py), pointer :: fptr => NULL()
  real a_real

  contains

  subroutine piggyback_libcloudphxx_lgrngn
    implicit none

    character(10) :: uid, pid

    if (associated(fptr) .eqv. .false.) then 
      ! assert for numerical precision
      if (sizeof(a_real) .ne. c_double) stop("DALES does not use double precision!")

      ! load pointer to Python micro_step() routine (getuid() and getpid() are GNU extensions)
      write (uid, "(I10.0)") getuid()
      write (pid, "(I10.0)") getpid()
      call load_ptr("/tmp/micro_step-" // trim(adjustl(uid)) // "-" // trim(adjustl(pid)) // ".ptr" // c_null_char,cptr)
 
      ! associate the C pointer with the F pointer
      call c_f_procpointer(cptr, fptr)
    end if

    if (rk3step /= 3) return 
    if (ntimee == 0) return ! spinup

    if (.not. fptr(                                        &
      rdt, dx, dy, dz,                                     &
      rhobf, size(rhobf, 1),                               &
      rhobh, size(rhobh, 1),                               &
      rhof,  size(rhof,  1),                               &
      u0,    size(u0,    1), size(u0,   2), size(u0,   3), &
      v0,    size(v0,    1), size(v0,   2), size(v0,   3), &
      w0,    size(w0,    1), size(w0,   2), size(w0,   3), &
      qt0,   size(qt0,   1), size(qt0,  2), size(qt0,  3), &
      thl0,  size(thl0,  1), size(thl0, 2), size(thl0, 3)  &
    )) stop("Error in Python!!!")
  end
end module

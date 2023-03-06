module common_parameter
  use kind_parameters
  implicit none 

  integer(ikind) ,parameter :: npar=1040000, nplink=75  !! NB, 200 for Quintic kernel, h=1.3dx,sup=3h
  !! NOTE: npmax should be npar*nplink
  integer(ikind) ,parameter :: npmax=npar*nplink
  integer(ikind) ,parameter :: dims = 2

  !! integer,parameter::nmax=500000

  real(rkind), parameter :: pi=3.14159265358979323846d0
!  real(rkind), parameter :: pi=4.0_rkind * DATAN(1.0_rkind)

end module common_parameter

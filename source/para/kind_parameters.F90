module kind_parameters

  use iso_c_binding

  implicit none

  private
  public :: rkind,ikind


  integer ,parameter :: rkind =c_double !selected_real_kind(digits)
  integer ,parameter :: ikind =c_int    !selected_int_kind(decades)

end module

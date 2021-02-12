module sphtools
  !! This module contains Smoothing kernels, gradients, and similar functions
  use kind_parameters
  use common_2d
  implicit none

  public :: fac, Wab

contains
!! ------------------------------------------------------------------------------------------------
  function fac(qq) result(factemp)
    ! gradient of kernel function
    !! NOTE: using intrinsic floor function allows if statements based on integer equalities, 
    !! which is faster than double precision inequalities.
  
    real(rkind), intent(in) :: qq
    real(rkind) ::  factemp
    integer(ikind) :: iselect
    real(rkind) ::  q3, q2, q1, q34, q24, q14

    iselect = floor(qq)
    if(iselect.eq.0)then
       q3=3_rkind-qq
       q2=2_rkind-qq
       q1=1_rkind-qq
       q34=q3**4
       q24=q2**4
       q14=q1**4
       factemp=ad_7h*(-5.0*q34+30.0*q24-75.0*q14) 
    elseif(iselect.eq.1)then
       q3=3_rkind-qq
       q2=2_rkind-qq
       q34=q3**4
       q24=q2**4
       factemp=ad_7h*(-5.0*q34+30.0*q24)
    elseif(iselect.eq.2)then
       q3=3_rkind-qq
       q34=q3**4
       factemp=ad_7h*(-5.0*q34)    
    else
       factemp=0_rkind
    endif

  end function fac
!! ------------------------------------------------------------------------------------------------
  function Wab(qq) result(fval)
  !! kernel function!
    real(rkind) ::  qq
    real(rkind) :: fval
    real(rkind) ::  q3, q2, q1, q35, q25, q15
    integer(ikind) :: iselect
    
    iselect = floor(qq)

    if(iselect.eq.0)then
       q3=3.0_rkind-qq
       q2=2.0_rkind-qq
       q1=1.0_rkind-qq
       q35=q3**5
       q25=q2**5
       q15=q1**5
       fval=ad_7*(q35-6.0_rkind*q25+15_rkind*q15)
    elseif(iselect.eq.1)then
       q3=3.0_rkind-qq
       q2=2.0_rkind-qq
       q35=q3**5
       q25=q2**5
       fval=ad_7*(q35-6.0_rkind*q25)
    elseif(iselect.eq.2)then
       q3=3.0_rkind-qq
       q35=q3**5
       fval=ad_7*q35
    else
       fval=0.0_rkind
    endif
  end function Wab
!! ------------------------------------------------------------------------------------------------
  function shifting_cohesion(qq) result(cval)
    !! A force based on Akinci's cohesion force, used in shifting sometimes
    real(rkind) :: qq
    real(rkind) :: cval
    real(rkind) :: q3,q13
    q3 = qq**3
    q13 = (1.0_rkind-qq)**3
    if(qq.gt.0.0_rkind .and. qq.le.0.5_rkind)then
       cval = 1.0/32.0
    elseif(qq.gt.0.5_rkind.and.qq.le.1.0_rkind)then
       cval = q13*q3
    else
       cval = 0.0_rkind
    end if
  end function
!! ------------------------------------------------------------------------------------------------
end module sphtools

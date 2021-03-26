#ifdef use_labfm
module labfm
  use kind_parameters
  use common_parameter
  use common_2d
  use sphtools
  use omp_lib
  implicit none
  
  private
  public :: labfm_weights

  real(rkind), parameter :: oosqrt2=1.0d0/dsqrt(2.0d0)
  real(rkind), parameter :: eta3 = 1.0d-14
  
contains
  subroutine labfm_weights
     integer(ikind) :: i,j,k
     real(rkind) :: rad,qq,x,y,xx,yy,det,tmp,kval,rad2,rad3
     real(rkind),dimension(dims) :: gradw,rij

     !! Linear system to find ABF coefficients
     real(rkind),dimension(:,:),allocatable :: amatGx,amatGy,amatL
     real(rkind),dimension(:),allocatable :: bvecGx,bvecGy,bvecL,gvec,xvec
     integer(ikind) :: i1,i2,nsize,nsizeG
     real(rkind) :: ff1,hh,xs,ys
     
     logical :: nearsurf

!! Choice of ABF type:: 1=original, 2=Hermite polynomials, 3=Legendre polynomials
!! ABFs 2, 3 and 4 are multiplied by an RBF (Wab(qq) set in sphtools).
#define ABF 2

     !! Set desired order::
     k=4
     nsizeG=(k*k+3*k)/2   !!  5,9,14,20,27,35,44... for k=2,3,4,5,6,7,8...

     !! Left hand sides and arrays for interparticle weights
     allocate(amatGx(nsizeG,nsizeG),amatGy(nsizeG,nsizeG),amatL(nsizeG,nsizeG))
     amatGx=0.0d0;amatGy=0.0d0;amatL=0.0d0
 
     !! Right hand sides, vectors of monomials and ABFs
     allocate(bvecGx(nsizeG),bvecGy(nsizeG),bvecL(nsizeG),gvec(nsizeG),xvec(nsizeG))
     bvecGx=0.0d0;bvecGy=0.0d0;bvecL=0.0d0;gvec=0.0d0;xvec=0.0d0

     !! No parallelism for individual linear system solving...
!     call openblas_set_num_threads(1)

     !$OMP PARALLEL DO PRIVATE(nsize,amatGx,k,j,rij,rad,qq,x,y,xx,yy,nearsurf, &
     !$OMP ff1,gvec,xvec,i1,i2,amatL,amatGy,bvecGx,bvecGy,bvecL,xs,ys)
     do i=1,npfb
        if(n_surf(i).eq.1) cycle
        ij_w_G(i,:,:)=0.0d0;ij_w_L(i,:)=0.0d0
        kgcm(i,1,1)=1.0d0;kgcm(i,2,2)=1.0d0;kgcm(i,1,2)=0.0d0;kgcm(i,2,1)=0.0d0        
        nsize = nsizeG
        amatGx=0.0d0
        nearsurf=.false.
        do k=1,ij_count(i)
           j = ij_link(i,k)
           if(n_surf(j).eq.1) nearsurf=.true. 
           rij(:) = r(i,:) - r(j,:)         
           x = -rij(1);y = -rij(2)

           !! Different types of ABF need different arguments (xx,yy)
           !! to account for domain of orthogonality
           rad = sqrt(x*x + y*y);qq=rad/h
#if ABF==1   
           ff1 = fac(qq)/h
           xx=x;yy=y
#elif ABF==2     
           ff1 = Wab(qq)
           xx=x/h;yy=y/h   !! Hermite   
#elif ABF==3
           ff1 = Wab(qq)
           xx=x/h/2.0d0;yy=y/h/2.0d0  !! Legendre
#elif ABF==4
           ff1 = Wab(qq)
           xx=x/h;yy=y/h    !! Laguerre   
#endif      
           !! Populate the ABF and monomial arrays
           gvec(1:5) = abfs12(rad,xx,yy,ff1);xvec(1:5) = monomials12(x,y)
           gvec(6:9) = abfs3(rad,xx,yy,ff1);xvec(6:9) = monomials3(x,y)
           gvec(10:14) = abfs4(rad,xx,yy,ff1);xvec(10:14) = monomials4(x,y)
          
                           
           !! Build the LHS - it is the same for all three operators
           do i1=1,nsize
              amatGy(i1,:) = xvec(i1)*gvec(:)   !! Contribution to LHS for this interaction
           end do    
           amatGx(:,:) = amatGx(:,:) + amatGy(:,:)                     
        end do   
        
        if(nearsurf)then
           amatGx(:,6:14)=0.0d0  ! right block
           amatGx(6:14,:)=0.0d0  ! bottom block
           do i1=6,14
              amatGx(i1,i1)=1.0d0
           end do
        endif

        amatGy = amatGx;amatL=amatGx !! copying LHS

        !! Build RHS for ddx and ddy
        bvecGx=0.0d0;bvecGx(1)=1.0d0
        bvecGy=0.0d0;bvecGy(2)=1.0d0        

        !! Solve system for grad coefficients
        i1=0;i2=0         
        call dgesv(nsize,1,amatGx,nsize,i1,bvecGx,nsize,i2)   
        i1=0;i2=0;nsize=nsizeG    
        call dgesv(nsize,1,amatGy,nsize,i1,bvecGy,nsize,i2)    

        !! Solve system for Lap coefficients
        bvecL(:)=0.0d0;bvecL(3)=1.0d0;bvecL(5)=1.0d0;i1=0;i2=0;nsize=nsizeG
        call dgesv(nsize,1,amatL,nsize,i1,bvecL,nsize,i2)     

        !! Another loop of neighbours to calculate interparticle weights
        do k=1,ij_count(i)
           j = ij_link(i,k) 
           rij(:) = r(i,:) - r(j,:)          
           x=-rij(1);y=-rij(2)

           !! Calculate arguments (diff ABFs need args over diff ranges)
           rad = sqrt(x*x + y*y);qq=rad/h
#if ABF==1   
           ff1 = fac(qq)/h
           xx=x;yy=y
#elif ABF==2     
           ff1 = Wab(qq)
           xx=x/h;yy=y/h    !! Hermite   
#elif ABF==3
           ff1 = Wab(qq)
           xx=x/h/2.0d0;yy=y/h/2.0d0  !! Legendre
#elif ABF==4
           ff1 = Wab(qq)
           xx=x/h;yy=y/h    !! Laguerre   
#endif          
           !! Populate the ABF array        
           gvec(1:5) = abfs12(rad,xx,yy,ff1)
           gvec(6:9) = abfs3(rad,xx,yy,ff1)
           gvec(10:14) = abfs4(rad,xx,yy,ff1)


           !! Weights for operators        
           ij_w_G(i,k,1) = dot_product(bvecGx,gvec)
           ij_w_G(i,k,2) = dot_product(bvecGy,gvec)
           ij_w_L(i,k) = -dot_product(bvecL,gvec)  
        end do

        wait(1) !! Temporary fix for a weird openblas/lapack/openmp bug 
        !! The bug may be due to a race condition or something, but not really sure
        !! What's the best way to fix this??
     end do
     !$OMP END PARALLEL DO
     deallocate(amatGx,amatGy,amatL,bvecGx,bvecGy,bvecL,gvec,xvec)     
     return
  end subroutine labfm_weights 
!! ------------------------------------------------------------------------------------------------
!! FUNCTIONS TO CALCULATE THE TAYLOR MONOMIALS
!! ------------------------------------------------------------------------------------------------
  function monomials4(x,y) result(cxvec)
     real(rkind),intent(in) :: x,y
     real(rkind),dimension(5) :: cxvec
     real(rkind) :: x2,y2,x3,y3,x4,y4
     x2=x*x;y2=y*y;x3=x2*x;y3=y2*y;x4=x3*x;y4=y3*y
     cxvec(1) = (1.0/24.0)*x4;cxvec(2)=(1.0/6.0)*x3*y;
     cxvec(3) = 0.25d0*x2*y2;cxvec(4)=(1.0/6.0)*x*y3;cxvec(5)=(1.0/24.0)*y4
  end function monomials4
  function monomials3(x,y) result(cxvec)
     real(rkind),intent(in) :: x,y
     real(rkind),dimension(4) :: cxvec
     real(rkind) :: x2,y2,x3,y3
     x2=x*x;y2=y*y;x3=x2*x;y3=y2*y 
     cxvec(1) = (1.0/6.0)*x3;cxvec(2) = 0.5d0*x2*y;cxvec(3) = 0.5d0*x*y2;cxvec(4) = (1.0/6.0)*y3;
  end function monomials3
  function monomials12(x,y) result(cxvec)
     real(rkind),intent(in) :: x,y
     real(rkind),dimension(5) :: cxvec
     real(rkind) :: x2,y2
     x2=x*x;y2=y*y
     cxvec(1) = x;cxvec(2)=y;cxvec(3)=0.5d0*x2;cxvec(4)=x*y;cxvec(5)=0.5d0*y2
  end function monomials12
!! ------------------------------------------------------------------------------------------------
#if ABF==1
!! ------------------------------------------------------------------------------------------------
!! ABFs generated from partial derivatives of an RBF
!! The "original" LABFM
!! ------------------------------------------------------------------------------------------------
  function abfs4(rad,x,y,ff1) result(ggvec)     !! FOUR
     real(rkind),intent(in) :: rad,x,y,ff1
     real(rkind),dimension(5) :: ggvec
     real(rkind) :: rad3,rad2,r7
     real(rkind) :: x2,y2,x3,y3,x4,y4
     rad2 = max(rad*rad,eta2);rad3=max(rad*rad*rad,eta3);r7=rad2*rad2*rad3
     x2=x*x;y2=y*y;x3=x2*x;y3=y2*y;x4=x3*x;y4=y3*y
     ggvec(1) = (12.0*x2*y2-3.0*y4)*ff1/r7                
     ggvec(2) = (9.0*x*y3-6.0*x3*y)*ff1/r7
     ggvec(3) = (2.0*x4-11.0*x2*y2+2.0*y4)*ff1/r7
     ggvec(4) = (9.0*x3*y-6.0*x*y3)*ff1/r7
     ggvec(5) = (12.0*x2*y2-3.0*x4)*ff1/r7
  end function abfs4
  function abfs3(rad,x,y,ff1) result(ggvec)     !! THREE
     real(rkind),intent(in) :: rad,x,y,ff1
     real(rkind),dimension(4) :: ggvec
     real(rkind) :: rad3,rad2,r5
     real(rkind) :: x2,y2,x3,y3
     rad2 = max(rad*rad,eta2);rad3=max(rad*rad*rad,eta3);r5=rad3*rad2
     x2=x*x;y2=y*y;x3=x2*x;y3=y2*y
     ggvec(1) = -3.0*y2*x*ff1/r5                 
     ggvec(2) = (2.0*x2*y - y3)*ff1/r5
     ggvec(3) = (2.0*x*y2 - x3)*ff1/r5
     ggvec(4) = -3.0*x2*y*ff1/r5
  end function abfs3
  function abfs12(rad,x,y,ff1) result(ggvec)     !! TWO AND ONE
     real(rkind),intent(in) :: rad,x,y,ff1
     real(rkind),dimension(5) :: ggvec
     real(rkind) :: rad3,rad2
     real(rkind) :: x2,y2
     rad2 = max(rad*rad,eta2);rad3=max(rad*rad*rad,eta3)
     x2=x*x;y2=y*y 
     ggvec(1) = ff1*x/max(rad,eta2) 
     ggvec(2) = ff1*y/max(rad,eta2)
     ggvec(3) = y2*ff1/rad3    
     ggvec(4) = -x*y*ff1/rad3
     ggvec(5) = x2*ff1/rad3
  end function abfs12
!! ------------------------------------------------------------------------------------------------
#elif ABF==2
!! HERMITE ABFs: Bivariate Hermite polynomials (probabalistic kind)
!! Formula based on Area, Dimitrov & Godoy (2015), J. Math. Anal. Appl. 421(1):830-841
!! with a=c=1, b=0.
!! Generated with: H_{p,q}(x,y) = H_{p}(x/sqrt(2))*H_{q}(y/sqrt(2))*2^((p+q)/2)
!!
!! NB sqrt2 and oosqrt are set in kind_parameters
!! NB ff1 multiplies the Hermite polynomial by an RBF to improve accuracy
!! ------------------------------------------------------------------------------------------------
  function abfs4(dummy,x,y,ff1) result(ggvec)      !! FOUR
     real(rkind),intent(in) :: x,y,ff1,dummy
     real(rkind),dimension(5) :: ggvec
     real(rkind) :: xx,yy
     xx = oosqrt2*x;yy=oosqrt2*y
     ggvec(1) = ff1*Hermite4(xx)*0.25d0
     ggvec(2) = ff1*Hermite3(xx)*Hermite1(yy)*0.25d0
     ggvec(3) = ff1*Hermite2(xx)*Hermite2(yy)*0.25d0
     ggvec(4) = ff1*Hermite1(xx)*Hermite3(yy)*0.25d0
     ggvec(5) = ff1*Hermite4(yy)*0.25d0
  end function abfs4
  function abfs3(dummy,x,y,ff1) result(ggvec)     !!! THREE
     real(rkind),intent(in) :: x,y,ff1,dummy
     real(rkind),dimension(4) :: ggvec
     real(rkind) :: xx,yy
     xx = oosqrt2*x;yy=oosqrt2*y
     ggvec(1) = ff1*Hermite3(xx)*oosqrt2*0.5d0
     ggvec(2) = ff1*Hermite2(xx)*Hermite1(yy)*oosqrt2*0.5d0
     ggvec(3) = ff1*Hermite1(xx)*Hermite2(yy)*oosqrt2*0.5d0
     ggvec(4) = ff1*Hermite3(yy)*oosqrt2*0.5d0
  end function abfs3
  function abfs12(dummy,x,y,ff1) result(ggvec)     !! TWO AND ONE
     real(rkind),intent(in) :: x,y,ff1,dummy
     real(rkind),dimension(5) :: ggvec
     real(rkind) :: xx,yy
     xx = x*oosqrt2;yy=y*oosqrt2
     ggvec(1) = ff1*Hermite1(xx)*oosqrt2
     ggvec(2) = ff1*Hermite1(yy)*oosqrt2
     ggvec(3) = ff1*Hermite2(xx)*0.5d0
     ggvec(4) = ff1*Hermite1(xx)*Hermite1(yy)*0.5d0
     ggvec(5) = ff1*Hermite2(yy)*0.5d0
  end function abfs12
!! ------------------------------------------------------------------------------------------------
#elif ABF==3
!! LEGENDRE ABFs: Legendre polynomials
!! NB ff1 multiplies the Legendre polynomial by an RBF to improve accuracy
!! ------------------------------------------------------------------------------------------------
  function abfs4(dummy,x,y,ff1) result(ggvec)      !! FOUR
     real(rkind),intent(in) :: x,y,dummy,ff1
     real(rkind),dimension(5) :: ggvec
     ggvec(1) = ff1*Legendre4(x)
     ggvec(2) = ff1*Legendre3(x)*Legendre1(y)
     ggvec(3) = ff1*Legendre2(x)*Legendre2(y)
     ggvec(4) = ff1*Legendre1(x)*Legendre3(y)
     ggvec(5) = ff1*Legendre4(y)
  end function abfs4
  function abfs3(dummy,x,y,ff1) result(ggvec)     !!! THREE
     real(rkind),intent(in) :: x,y,dummy,ff1
     real(rkind),dimension(4) :: ggvec
     ggvec(1) = ff1*Legendre3(x)
     ggvec(2) = ff1*Legendre2(x)*Legendre1(y)
     ggvec(3) = ff1*Legendre1(x)*Legendre2(y)
     ggvec(4) = ff1*Legendre3(y)
  end function abfs3
  function abfs12(dummy,x,y,ff1) result(ggvec)     !! TWO AND ONE
     real(rkind),intent(in) :: x,y,dummy,ff1
     real(rkind),dimension(5) :: ggvec
     ggvec(1) = ff1*Legendre1(x)
     ggvec(2) = ff1*Legendre1(y)
     ggvec(3) = ff1*Legendre2(x)
     ggvec(4) = ff1*Legendre1(x)*Legendre1(y)
     ggvec(5) = ff1*Legendre2(y)
  end function abfs12
!! ------------------------------------------------------------------------------------------------
#elif ABF==4  
!! LAGUERRE ABFs: Laguerre polynomials
!! NB ff1 multiplies the Laguerre polynomial by an RBF to improve accuracy
!! ------------------------------------------------------------------------------------------------
  function abfs4(dummy,x,y,ff1) result(ggvec)      !! FOUR
     real(rkind),intent(in) :: x,y,ff1,dummy
     real(rkind),dimension(5) :: ggvec
     real(rkind) :: xx,yy
     xx = oosqrt2*x;yy=oosqrt2*y
     ggvec(1) = ff1*Laguerre4(xx)*0.25d0
     ggvec(2) = ff1*Laguerre3(xx)*Laguerre1(yy)*0.25d0
     ggvec(3) = ff1*Laguerre2(xx)*Laguerre2(yy)*0.25d0
     ggvec(4) = ff1*Laguerre1(xx)*Laguerre3(yy)*0.25d0
     ggvec(5) = ff1*Laguerre4(yy)*0.25d0
  end function abfs4
  function abfs3(dummy,x,y,ff1) result(ggvec)     !!! THREE
     real(rkind),intent(in) :: x,y,ff1,dummy
     real(rkind),dimension(4) :: ggvec
     real(rkind) :: xx,yy
     xx = oosqrt2*x;yy=oosqrt2*y
     ggvec(1) = ff1*Laguerre3(xx)*oosqrt2*0.5d0
     ggvec(2) = ff1*Laguerre2(xx)*Laguerre1(yy)*oosqrt2*0.5d0
     ggvec(3) = ff1*Laguerre1(xx)*Laguerre2(yy)*oosqrt2*0.5d0
     ggvec(4) = ff1*Laguerre3(yy)*oosqrt2*0.5d0
  end function abfs3
  function abfs12(dummy,x,y,ff1) result(ggvec)     !! TWO AND ONE
     real(rkind),intent(in) :: x,y,ff1,dummy
     real(rkind),dimension(5) :: ggvec
     real(rkind) :: xx,yy
     xx = x*oosqrt2;yy=y*oosqrt2
     ggvec(1) = ff1*Laguerre1(xx)*oosqrt2
     ggvec(2) = ff1*Laguerre1(yy)*oosqrt2
     ggvec(3) = ff1*Laguerre2(xx)*0.5d0
     ggvec(4) = ff1*Laguerre1(xx)*Laguerre1(yy)*0.5d0
     ggvec(5) = ff1*Laguerre2(yy)*0.5d0
  end function abfs12
#endif
!! ------------------------------------------------------------------------------------------------
!! Below this line are functions which return univariate polynomials from which some ABFs above
!! are constructed
!! ------------------------------------------------------------------------------------------------
  !! Univariate Hermite Polynomials (Physicists kind)
  function Hermite1(z) result(Hres)
     real(rkind),intent(in) :: z
     real(rkind) :: Hres
     Hres = 2.0d0*z
  end function Hermite1
  function Hermite2(z) result(Hres)
     real(rkind),intent(in) :: z
     real(rkind) :: Hres
     Hres = 4.0d0*z*z - 2.0d0
  end function Hermite2
  function Hermite3(z) result(Hres)
     real(rkind),intent(in) :: z
     real(rkind) :: Hres
     Hres = 8.0d0*z*z*z - 12.0d0*z
  end function Hermite3
  function Hermite4(z) result(Hres)
     real(rkind),intent(in) :: z
     real(rkind) :: Hres
     Hres = 16.0d0*z*z*z*z - 48.0d0*z*z + 12.0d0
  end function Hermite4
!! ------------------------------------------------------------------------------------------------
!! Univariate Legendre polynomials
!! ------------------------------------------------------------------------------------------------
  function Legendre1(z) result(Lres)
     real(rkind),intent(in) :: z
     real(rkind) :: Lres
     Lres = z
  end function Legendre1
  function Legendre2(z) result(Lres)
     real(rkind),intent(in) :: z
     real(rkind) :: Lres
     Lres = 0.5d0*(3.0d0*z*z-1.0d0)
  end function Legendre2
  function Legendre3(z) result(Lres)
     real(rkind),intent(in) :: z
     real(rkind) :: Lres
     Lres = 0.5d0*(5.0d0*z*z*z-3.0d0*z)
  end function Legendre3
  function Legendre4(z) result(Lres)
     real(rkind),intent(in) :: z
     real(rkind) :: Lres
     real(rkind) :: z2
     z2=z*z
     Lres = 0.125d0*(35.0d0*z2*z2 - 30.0d0*z2 + 3.0d0)
  end function Legendre4
!! ------------------------------------------------------------------------------------------------
!! Univariate Laguerre polynomials
!! ------------------------------------------------------------------------------------------------
  function Laguerre1(z) result(Lres)
     real(rkind),intent(in) :: z
     real(rkind) :: Lres
     Lres = 1.0d0 - z
  end function Laguerre1
  function Laguerre2(z) result(Lres)
     real(rkind),intent(in) :: z
     real(rkind) :: Lres
     Lres = 0.5d0*(z*z-4.0d0*z+2.0d0)
  end function Laguerre2
  function Laguerre3(z) result(Lres)
     real(rkind),intent(in) :: z
     real(rkind) :: Lres
     Lres = (-z*z*z+9.0d0*z*z-18.0d0*z+6.0d0)/6.0d0
  end function Laguerre3
  function Laguerre4(z) result(Lres)
     real(rkind),intent(in) :: z
     real(rkind) :: Lres
     real(rkind) :: z2
     z2=z*z
     Lres = (z2*z2 - 16.0d0*z2*z + 72.0d0*z2 - 96.0d0*z + 24.0d0)/24.0d0
  end function Laguerre4
!! ------------------------------------------------------------------------------------------------
end module labfm
#endif

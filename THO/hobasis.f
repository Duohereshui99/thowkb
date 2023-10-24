!!!!!!两种3dbasis的写法
      module hobasis 
                use parameter
                contains

      function normho(nu,n,l)
        implicit none
        real*8:: nu,pi,normho
        integer:: n,l,m,p
        pi=acos(-1d0)

        normho=sqrt(sqrt(2*nu**3/pi)*2**(n+2*l+3)*fact(n)
     &  *nu**l/doublefact(2*n+2*l+1))
        if(normho<1e-6) then
           write(*,*)'nu,n,l,normho',nu,n,l,normho
           write(*,*)'fact(n)',fact(n)
           write(*,*)'fact(n+l)=',fact(n+l)
           write(*,*)'fact(2n+2l+1)=',fact(2*n+2*l+1)
           stop
        endif
      end function normho
      !!cccccccccccccc

          function fact(n)           
            implicit none
            integer n
            real*8 fact, dgamma,x
            x=dfloat(n+1)
            fact=dgamma(x)
          end function fact
ccccccc
      real*8 function doublefact(n)!双阶乘
            implicit none
            integer::n,i
            real::s
            s=1.0
            if(mod(n,2)==0) then
                  do i=n,2,-2
                        s=s*i
                  end do
            else
                  do i=n,1,-2
                        s=s*i
                  end do
            end if
            doublefact=s
      end function
ccccccc
      real*8 function ho3d(n,l,nu,r)            !!标准形式的三维BASIS
      implicit none                           !!修改后的nu取mu*omega/(2hbar)
      integer l,n,p
      real*8::r,tp,norma,nu
      norma=normho(nu,n,l)
!        if (norma<1e-6) then
!           write(*,*)'ho3d: Norm=0!!!for  nu,n,l',nu,n,l
!        endif
      ho3d=norma*r**l*dexp(-nu*r**2)*
     & general_laguerre(2*nu*r**2,n,l+0.5d0) 
      end function ho3d
ccccccc
ccccccc
cccccccccccccccccccccccccccccccccccccccccccccc    
c     Generalized Laguerre function L(n,l+1/2,x)
c *** -----------------------------------------------
      function laguerre(n,l,x)
        implicit none
        integer n,l,p
        real*8:: x,eps,tp,laguerre
        parameter(eps=1e-6)

        if(x<eps) x=eps
        tp=(-x)**n/fact(n)

        do p=1,n
!           write(*,*)'p=',p
           tp=tp-(n+l+1.5-p)*(n+1-p)/p/x*tp
        enddo
        laguerre=tp
        end function laguerre
ccccccc

      recursive function general_laguerre(x, k, alpha) result(L)
      real*8, intent(in) :: x
      integer, intent(in) :: k
      real*8, intent(in) :: alpha
      real*8 :: L

      if (k == 0) then
         L = 1.0
      else if (k == 1) then
         L = 1.0 + alpha - x
      else
         L = ((2*k - 1 + alpha - x) * general_laguerre(x, k-1, alpha) 
     &    - (k - 1 + alpha) * general_laguerre(x, k-2, alpha)) / k
      end if

      end function general_laguerre
        end module

cccccccccccccc
!         module HOBASIS_3D
!             use parameter
!             use sp_f
!             contains
!             function fact(n)           
!                   implicit none
!                   integer n
!                   real*8 fact, dgamma,x
!                   x=dfloat(n+1)
!                   fact=dgamma(x)
!                 end function fact
! ccccccc
!             real*8 function doublefact(n)!双阶乘
!                   implicit none
!                   integer::n,i
!                   real::s
!                   s=1.0
!                   if(mod(n,2)==0) then
!                         do i=n,2,-2
!                               s=s*i
!                         end do
!                   else
!                         do i=n,1,-2
!                               s=s*i
!                         end do
!                   end if
!                   doublefact=s
!             end function
! ccccccc!gsl库超几何函数形式的三维basis
!        real*8 function HOBASIS_3D(n,l,alpha,z)
!             implicit none
!             integer::n,l
!             real*8::z,result,alpha
!             result=alpha**(3.d0/2.d0)*sqrt(2**((l+2-n)*1.d0)
!      & *doublefact(2*l+2*n+1)*1.d0
!      & /sqrt(PI)/(fact(n)*1.d0)/(doublefact(2*l+1))**2*1.d0)
!      & *(alpha*z)**l*exp(-alpha**2*z*z/2.0)
!      & *hg(-n*1.0d0,(l+1.5)*1.0d0,alpha**2*z*z*1.0d0)
!             HOBASIS_3D=result
!             end function
!       end module
! ccccccc



!       module sp_f       !特殊函数,如laguerre多项式,以及HYPERgeometry
!             use iso_c_binding
!             ! Declare external C function
!             interface

! !       real*8 function Laguerre(n,a,x) !n是整数,a和x是实数
! !      &      bind(C, name="gsl_sf_laguerre_n") 
! !             import :: c_double,c_int
! !             integer(c_int),value::n
! !             real(c_double),value::a,x
! !       end function
! ccccccc
!       real*8 function hg(a,b,x)   !超几何函数
!      &      bind(C, name="gsl_sf_hyperg_1F1") 
!             import :: c_double
!             real(c_double),value::a,b,x
!       end function
! ccccccc
!             end interface
!       end module
!求导插值积分等算法
ccccccc
      module algorithm      !一些用到的函数,把它们module在一个func里面
        contains
ccccccc
! interpolation function for uniform grids 
!Y:求第Y个格点处的函数值,一般是小数,
!F:函数值数列, N:有N个格点
!该函数是直接copy的

      function FFR4(Y,F,N)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 F(N),P,P1,P2,Q,X,FFR4
      REAL*8 Y
      PARAMETER(X=.16666666666667)
      P=Y
      I=P
ccccccc
      IF(I.LE.0) GO TO 2
      IF(I.GE.N-2) GO TO 4
    1 P=P-I
      P1=P-1.
      P2=P-2.
      Q=P+1.
      FFR4=(-P2*F(I)+Q*F(I+3))*P*P1*X+(P1*F(I+1)-P*F(I+2))*Q*P2*.5
      RETURN
    2 IF(I.LT.0) GO TO 3
      I=1
      GO TO 1
    3 FFR4=F(1)
      RETURN
    4 IF(I.GT.N-2) GO TO 5
      I=N-3
      GO TO 1
    5 FFR4=F(N)
      RETURN
      end function
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!高斯-勒让德积分。 N:积分格点数； x1,x2:积分下上限    
!输出x,w:格点位置和权重。
!gauleg相当于只对自变量区间作用,做一个重新离散化
!gauleg的核心是将这个新离散化的自变量数组给一组新的权重
!用这个新给的权重计算更稳定,这是重排的原因
!这个权重相当于dx(n)
!该函数是copy的
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE gauleg(N,x1,x2,X,W)
        IMPLICIT NONE
        INTEGER N
        REAL*8 x1,x2,X(N),W(N)
        REAL*8 z1,z,xm,xl,pp,p3,p2,p1,pi,tol
        INTEGER m,i,j

        pi=acos(-1.0)
        tol=1.E-12

        m=(n+1)/2
        xm=0.5*(x2+x1)
        xl=0.5*(x2-x1)

         DO 10 i=1,m
         z=cos(pi*(i-0.25)/(N+0.5))

 20      CONTINUE
         p1=1.0E0
         p2=0.0E0
         DO 30 j=1,N
          p3=p2
          p2=p1
          p1=((2*j-1)*z*p2-(j-1)*p3)/j
 30      CONTINUE
         pp=N*(z*p1-p2)/(z*z-1.0E0)
         z1=z
         z=z1-p1/pp
         IF( abs(z1-z) .GT. tol) GOTO 20 ! Scheifenende

         X(i) = xm - xl*z
         X(n+1-i) = xm + xl*z
         W(i) = 2.E0*xl/((1.0-z*z)*pp*pp)
         W(n+1-i) = W(i)
 10     CONTINUE
        END SUBROUTINE gauleg 
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
! 下面算二阶导数用了一个简单的五点差分
! 形式选取了subroutine的形式
! 其输入值依次为:函数（值）数组y,
! 以及计划输出的二阶导数数组d2y,
! 数组的size,n；以及步长dx（这里是均匀步长）
! 其中只有d2y是一个空的数组用来往里面填写并输出  
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine second_derivative(y,d2y,n,dx)
       implicit none
       integer,intent(in)::n
       real*8,intent(in)::dx
       real*8,dimension(1:n),intent(in)::y
       real*8,dimension(1:n),intent(out)::d2y
       integer::i 
       d2y(1)=(35.d0/12.d0*y(1)-26.d0/3.d0*y(2)+19.d0/2.d0*y(3)
     &   -14.d0/3.d0*y(4)+11.d0/12.d0*y(5))/(dx**2)
       d2y(2)=(11.d0/12.d0*y(1)-5.d0/3.d0*y(2)+1.d0/2.d0*y(3)
     &  +1.d0/3.d0*y(4)-1.d0/12.d0*y(5))/(dx**2)
       d2y(n-1)=(-1.d0/12.d0*y(N-4)+1.d0/3.d0*y(N-3)+1.d0/2.d0*y(N-2)
     &  -5.d0/3.d0*y(N-1)+11.d0/12.d0*y(N))/(dx**2)
       d2y(n)=(11.d0/12.d0*y(N-4)-14.d0/3.d0*y(N-3)+19.d0/2.d0*y(N-2)
     &  -26.d0/3.d0*y(N-1)+35.d0/12.d0*y(N))/(dx**2)
       do i=3,n-2
       d2y(i)=(-y(i-2)+16.d0*y(i-1)-30.d0*y(i)+ 
     & 16.d0*y(i+1)-y(i+2))/(12.d0*dx**2)
          !  write(*,*) d2y(i)
       end do
       end subroutine second_derivative
ccccccc


c *** Calculate d^2u(r)/dr^2 using five points derivative formula
c     f(ndim)=function to make derivative
c     h      =step
c     j      =point for derivative
      function deriv1(f,h,ndim,j)
        implicit none
        integer ndim,j
        real*8 f(ndim),h,deriv1

        if ((j.eq.1).or.(j.eq.2)) then
           deriv1=(-f(j+2)+4d0*f(j+1)-3d0*f(j))/2d0/h
        else if (j.eq.ndim-1) then
           deriv1=(3d0*f(j)-4d0*f(j-1)+f(j-2))/2d0/h
        else if (j.eq.ndim) then
           deriv1=0 !!!CHECK
        else ! five points formula
           deriv1=(f(j-2)-8*f(j-1)+8*f(j+1)-f(j+2))/h/12.
         end if
      end function deriv1
ccccccc
      end module
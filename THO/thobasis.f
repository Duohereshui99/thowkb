        module thobasis
            use parameter
            use hobasis
            contains
ccccccc       
      real*8 function LSTFUN(gamma,m,r) !LST变换,相当于公式里面的s(r)
      implicit none
      real*8::gamma,m,r
      LSTFUN=(1/((1.d0/r)**m+
     & (1.d0/gamma/sqrt(r))**m))**(1.d0/m)
      end function
ccccccc
      real*8 function D1LSTFUN(gamma,m,r)     !s(r)的一阶导数的解析表达式
      implicit none
      real*8::gamma,m,r
      real*8::x,y       !计算的中间量
      x=(1.d0/r)**(m+1.d0)+1.d0/(2.0*r)*(1.d0/(gamma*sqrt(r)))**m
      y=(1.d0/r)**m+(1.d0/(gamma*sqrt(r)))**m
      D1LSTFUN=LSTFUN(gamma,m,r)*x/y
      end function
ccccccc
      real*8 function THOFUNC(n,l,alpha,gamma,m,r)    !alpha是HOBASIS的无量纲参数
      implicit none
      integer::n,l
      real*8::alpha,gamma,m,r
      THOFUNC=sqrt(D1LSTFUN(gamma,m,r))
     & *ho3d(n,l,alpha,LSTFUN(gamma,m,r))
     & *LSTFUN(gamma,m,r)/r !再乘s(r)
      end function
ccccccc
      real*8 function f(beta,r)          !!算THO相移用到的好性质的函数f(r)
      implicit none
      real*8::r,beta
      f=1-exp(-beta*r**2)
      end function
ccccccc
      real*8 function d2f(beta,r)        !!算THO相移用到的好性质的函数f(r)的二阶导数
      implicit none
      real*8::r,beta
      d2f=(2*beta-4*beta**2*r**2)*exp(-beta*r**2)
      end function
ccccccc
        end module

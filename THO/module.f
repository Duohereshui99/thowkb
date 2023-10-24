ccccccc
      module mesh                              !mesh:积分的一些常量
            implicit none
            integer::n_int!=200       !积分格点数
            integer::n_diff!=12000    !微分格点数  
            real*8::hcm!=0.005        !微分格点的步长   !!n_diff*hcm=总的范围
      end module
ccccccc
      module system   !系统的一些量,这里是一些np的两体信息
            implicit none          
            real*8::z_1!=1d0           !带电量
            real*8::z_2!=0d0
            real*8::z12                !带电量乘积
            real*8::mass_1!=1.0078     !质量数
            real*8::mass_2!=1.0087
            real*8::mu                          !约化质量
            integer::L!=0              !角动量量子数
            real*8::Q               !Q value
            real*8::F0
            real*8::gamma0          !width
            real*8::t_half0         !t half
            real*8::P                !   P factor
      end module
ccccccc
      module lstpar     !THO方法LST变换以及变分计算中用到的一些参数
            implicit none
            real*8::gamma!=1.5       !!LST的gamma参数
            real*8::m!=4             !!LST的m参数
            integer::n_basis!=14       !基的个数N=n_basis+1
            !real*8,parameter::omega=0.16414517233728546         !!谐振子势参数ω
            real*8::b!=1.6           !!b为理论上LST中用于变分的参数
            real*8,parameter::beta=0.01       !!beta为计算THO相移f(r)函数的参数
            real*8::alpha                     !alpha参数,3DHOBASIS的无量纲参数
      end module
ccccccc
      module potential !!势相关的参数
            implicit none
            real*8::v0!=-72.150       !!gausspot阱深
            real*8::r0!=0             !!gausspot中心
            real*8::a!=1.484          !!gausspot宽度
      end module
ccccccc
      module scatmatch                          !!散射match的部分
            use mesh
            use system
            implicit none
            complex*16::hlp,hln,dhlp,dhln             !hankel及其一阶导数
            complex*16::Sl                            !S矩阵
            complex*16::dy                            !match处trial wf的导数
            complex*16::c                             !match常数
            real*8::k2                                !选取能量的波矢的模方 
            real*8::eta                               !Sommerfeld常数
            real*8,parameter::XLMIN=0d0               !COUL90输入参数
            integer::IFAIL                            !COUL90参数  
            integer,parameter::KFN=0                  !COUL90参数
            integer::n_match                          !match处的index 
      end module
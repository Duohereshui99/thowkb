ccccccc
      program main
            use parameter     !使用相关的物理常数
            use system        !将两体系统的信息包括在内
            use algorithm     !求导插值积分等算法
            use mesh          !格点信息
            use coulfunc      !!coul90
            use lstpar        !!LST变换的参数
            use variable      !!变量数组模块
            use scatmatch     !!散射match变量模块
            use potential     !!势的相关参数变量
            use potentialfunction !!势的具体函数
            use thobasis
            implicit none
            real*8::r,s,t           !求和变量
            real*8::t1,t2           !cputime,查看运行时长
            integer::i,j,k,ll,u     !计数量,在积分的时候用
ccccccc
            namelist /systems/ z_1,z_2,mass_1,mass_2,L
            namelist /lst/ gamma,m,b,n_basis
            namelist /gaussianpot/ v0,r0,a
            namelist /meshes/ n_int,n_diff,hcm
ccccccc
            open(777,file='test.in')
            read(777,nml=systems)
            read(777,nml=lst)
            read(777,nml=gaussianpot)
            read(777,nml=meshes)
ccccccc
ccccccc
            call cpu_time(t1)
ccccccc
            allocate(psi(0:n_basis,1:n_int))         !basis及其二阶导数
            allocate(d2psi(0:n_basis,1:n_int))     
            allocate(psi_1(1:n_diff))       !插值后的basis及其二阶导数
            allocate(d2psi_1(1:n_diff))
            allocate(vpot(1:n_int),vpot_1(1:n_diff))                 !势
            allocate(H(0:n_basis,0:n_basis))        !Hamiltonian矩阵
            allocate(rr(1:n_int),rrw(1:n_int))               !积分点和权重 
            allocate(wr(0:n_basis),wi(0:n_basis))              !特征值的实部和虚部
            allocate(vr(0:n_basis,0:n_basis))       !每个特征值对应的特征向量
            allocate(PHI(1:n_int))                  !最终的波函数(基态bound)
            allocate(z(0:n_basis))              !特征值复数形式
            allocate(fr(1:n_diff))                  !均匀格点f(r)=1-exp(-beta*r**2)
            allocate(d2fr(1:n_diff))               !均匀格点f(r)的二阶导数
            allocate(fr_1(1:n_int))                !gauss积分点f(r)=1-exp(-beta*r**2)
            allocate(d2fr_1(1:n_int))               !gauss积分点f(r)的二阶导数
            allocate(FF(1:n_diff))                 !均匀格点Coulomb函数F
            allocate(GG(1:n_diff))                  !均匀格点Coulomb函数G
            allocate(FF_1(1:n_int))                 !gauss积分点Coulomb函数F
            allocate(GG_1(1:n_int)) 
            allocate(scatwf(1:n_int))              !散射波函数
            allocate(FC(0:L),GC(0:L))              !coul90的FC,GC
            allocate(FCP(0:L),GCP(0:L))
            allocate(fff(n_int))
            allocate(rrr(20))
ccccccc
            z12=z_1*z_2
            mu=amu*(mass_2*mass_1)/(mass_1+mass_2)              !约化质量
            alpha=1d0/2d0/b**2!sqrt(mu*omega/hbarc)             !HOBASIS的无量纲参数α=mu*omega/hbar
            !1/b**2=mu*omega/hbar,而basis的标准形式下的系数是mu*omega/2/hbar
            write(*,*) 'omega=',hbarc/mu/b**2
            write(*,*) 'nbasis=:',n_basis
!用前面的subroutine来获得basis及其二阶导数,还有势
!以及积分格点和权重
            call gauleg(n_int,0.d0,hcm*n_diff,rr,rrw)
ccccccc
            do i=1,n_diff  !势的用来做FFR4的数组,一开始的size取n_diff
            r=i*hcm 
            vpot_1(i) = vpot1(v0,a,r0,z12,mu,l,r)  
            end do 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!下面是两层循环,相当于i行j列的矩阵,第i行表示第i个基
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do i=0,n_basis
            do j=1,n_diff     !U(r)=rR(r)
            r=j*hcm
            psi_1(j)=THOFUNC(i,L,alpha,gamma,m,r)*r
            end do
            call second_derivative(psi_1,d2psi_1,n_diff,hcm)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
!下面可以直接用它做插值
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            do j=1,n_int
            psi(i,j)=FFR4(rr(j)/hcm,psi_1,n_diff)
            d2psi(i,j)=FFR4(rr(j)/hcm,d2psi_1,n_diff)
            end do 
      end do
ccccccc
            do i=1,n_int
                  vpot(i)=FFR4(rr(i)/hcm,vpot_1,n_diff)
            end do
ccccccc


!积分算矩阵元
            do i=0,n_basis
                  do j=0,n_basis                 
                  do k=1,n_int
                     H(i,j)=H(i,j)+psi(i,k)*vpot(k)*psi(j,k)*rrw(k)    !vpot项
     &               +(-hbarc*hbarc/2/mu)*rrw(k)*psi(i,k)*d2psi(j,k)   !动能项
                  end do
                  end do
            end do


            do i=0,n_basis
                  write(22,*) H(i,:)      !H矩阵在fort.22里面
            end do
            call eigenvalue(n_basis+1,H,wr,wi,vr)
            
            ll=0;  !循环遍历查找最小特征值的下标,对应特征向量列数的index
            do i=0,n_basis
                  if(wr(i)<wr(ll)) then
                        ll=i
                  end if
            end do

            write(*,*)"egs:",minval(wr)
            write(*,*)"index of egs:",ll 


ccccccc
            u=0     !循环遍历查找最小正特征值的下标,即找连续态的index
            do i=0,n_basis
                  if(wr(i)<wr(u).and.wr(i)>0.and.wr(i)<10) then
                        u=i
                  end if
            end do
            !!u表示最小正特征值或者是我们手动选取的下标
            u=10
            write(*,*)'minimum positive eigenvalue selected:',wr(u)
            write(*,*)'index of the positive eigenvalue:',u
ccccccc
            Q=wr(u)
            write(*,*)'Q=',Q
             
            z=cmplx(wr,wi)    !将两个实数数组转化为complex型

            do i=0,n_basis
                  write(23,*) z(i)        !复特征值在fort.23里面
                  write(24,*) vr(:,i)     !特征向量,按列write,fort.24
            end do
            
            do i=1,n_int
                  s=0
                  do j=0,n_basis
                        s=s+vr(j,u)*psi(j,i)    !第二个指标u对应了wr的index,即第几列特征向量
                  end do
                  PHI(i)=s                      !这里的index i对应了rr的index
                  write(25,*) rr(i),PHI(i)      !最终径向波函数U=rR,fort.25
                  write(26,*) rr(i),vpot(i)     !势,fort.26
            end do                              !PHI(i)对应rr(i)
            t=0
            do i=1,n_int
               t=t+PHI(i)**2*rrw(i)   
            end do
            write(*,*)'norm of ps:',t
ccccccc
       do i=1,n_int        !!k^2=2mu/h^2(Q-V(r))=k^2,|k|=sqrt(abs(k^2))
          fff(i)=kr(Q,mu,v0,a,r0,z12,l,rr(i))
        end do
        k=1
        do i=1,n_int-1      !!where Q=V(r),r(i) -> rr(r(i)) -> r_i
          if(fff(i)*fff(i+1)<0) then 
            rrr(k)=i
            k=k+1
          end if
        end do
        write(*,*) 'the number of radial points where Q=V(r):',k-1
        write(*,*) rr(rrr(1)),rr(rrr(2)),rr(rrr(3))

        F0=0d0 
          do i=rrr(1),rrr(2)
            s=0d0 
            do j=rrr(1),i           !!int_r1^r dr'
              s=s+rrw(j)*sqrt(abs(fff(j)))
            end do                                    
            F0=F0+rrw(i)*(cos(s-pi/4d0))**2/sqrt(abs(fff(i)))     !!int_r1^r2 dr
          end do
          !!F goes reciprocally to get normalization factor F
          F0=1d0/F0
          write(*,*) 'normalization factor F:',F0

          t=0d0
          do i=rrr(2),rrr(3)
            t=t+sqrt(abs(fff(i)))*rrw(i)
          end do

          P=0d0
          do i=rrr(1),rrr(2)
            P=P+rrw(i)*PHI(i)**2*rrw(i)
          end do
          write(*,*) 'P factor:' ,P
          write(*,*) 't:' ,t
          gamma=P*F0*hbarc**2/4d0/mu*exp(-2d0*t)
          t_half0=hbarc*log(2d0)/gamma                   !!fm
          t_half0=t_half0*ratio                             !!s

          write(*,*) 'gamma:',gamma,'MeV'
          write(*,*) 'half life T:',t_half0,'s'
          write(*,*) 'half life T:',t_half0/3600d0/24d0/365d0,'yr'




!!计算半径的方均根
            ! s=0
            ! do i=1,n_int
            !       s=s+PHI(i)**2*rrw(i)*rr(i)**2
            ! end do
            ! write(*,*)'rms=',sqrt(s)
ccccccc

ccccccc
            ! do i=1,n_int                        !!从势消失处定义match点,判据为势比较小的地方
            !       if(abs(vpot(i))<3e-3) then 
            !             n_match=i
            !             exit
            !       end if
            ! end do
ccccccc           
!比如我们选取最小正特征能量来算散射波函数,首先算它的k^2
ccccccc
            k2=2d0*mu*wr(u)/hbarc**2      !!k^2,这里取的指标是u,需和前面的u对应
                                          !!k用sqrt(k2)表示
            eta=mu*Z_1*Z_2*e2/hbarc**2/sqrt(k2) !eta
!              CALL COUL90(sqrt(k2)*rr(n_match),eta,XLMIN,l,FC,GC,FCP,GCP
!      &      ,KFN,IFAIL)
! ccccccc
!              if(IFAIL /= 0 ) then
!                    write(*,*)"error when call coulomb function!"
!              end if
! ccccccc
!              hlp=complex(GC(l),FC(l))
!              hln=complex(GC(l),-FC(l))
!              dhlp=complex(GCP(l),FCP(l))
!              dhln=complex(GCP(l),-FCP(l))
! ccccccc
! !!match处的trial wf的导数先用前或后向差分方便算一下,如果误差比较大再去改别的算法          
!              dy=(PHI(n_match+1)-PHI(n_match))/
!      &       (rr(n_match+1)-rr(n_match))
! ccccccc  S矩阵元S_{l}
!              Sl=(dy*hln-PHI(n_match)*dhln*sqrt(k2))/
!      &      (dy*hlp-PHI(n_match)*dhlp*sqrt(k2))
!              write(*,*) 'L=',L,'Sl=',Sl
!       !       write(999,*) wr(u),Sl
!       !      end do

!              c=(hln-Sl*hlp)*(0d0,1d0)/2d0/PHI(n_match)

!              scatwf=c*PHI

!              do i=1,n_int
!                    write(27,*) rr(i),real(scatwf(i))
!                    write(28,*) rr(i),aimag(scatwf(i))
!                    write(29,*) rr(i),abs(scatwf(i))
!              end do
ccccccc
!             !!下面用给定的公式计算相移
!             !!先写好f(r)及其二阶导数在gauss格点处的值
!             do i=1,n_diff
!                   fr(i)=f(beta,i*hcm)
!                   d2fr(i)=d2f(beta,i*hcm)                  
!             end do      
! ccccccc
!             do i=1,n_int
!                   fr_1(i)=FFR4(rr(i)/hcm,fr,n_diff)
!                   d2fr_1(i)=FFR4(rr(i)/hcm,d2fr,n_diff)
!             end do

!             do i=1,n_diff                       !给COULOMB函数的均匀格点值
!             CALL COUL90(sqrt(k2)*i*hcm,eta,XLMIN,l,FC,GC,FCP,GCP
!      &      ,KFN,IFAIL)
!             FF(i)=FC(l)
!             GG(i)=GC(l)
!             end do

!             do i=1,n_int                       !给COULOMB函数的gauss格点值
!             FF_1(i)=FFR4(rr(i)/hcm,FF,n_diff)
!             GG_1(i)=FFR4(rr(i)/hcm,GG,n_diff)
!             end do

!             s=0                                 !!计算tanδ的分子
!             do i=1,n_int                        !!wr选的指标是wr(u)
!                  s=s+PHI(i)*wr(u)*fr_1(i)*FF_1(i)*rrw(i)
!      &           -PHI(i)*(-hbarc*hbarc/2/mu)*d2fr_1(i)*FF_1(i)*rrw(i)
!      &           -PHI(i)*vpot(i)*fr_1(i)*FF_1(i)*rrw(i)
!             end do
!             t=0
!             do i=1,n_int                        !!计算tanδ的分母
!                  t=t+PHI(i)*wr(u)*fr_1(i)*GG_1(i)*rrw(i)
!      &           -PHI(i)*(-hbarc*hbarc/2/mu)*d2fr_1(i)*GG_1(i)*rrw(i)
!      &           -PHI(i)*vpot(i)*fr_1(i)*GG_1(i)*rrw(i)
!             end do
!             write(*,*)'tanδ=',-s/t
!             write(*,*)'δ=',atan(-s/t)
!             write(*,*)'Sl=',exp(-2*(0d0,1d0)*atan(-s/t))
            deallocate(fff,rrr)
            deallocate(psi,d2psi,vpot,H,rr,rrw,wr,wi,vr)
            deallocate(vpot_1)
            deallocate(PHI,z,fr,d2fr)
            deallocate(psi_1,d2psi_1)
        !    deallocate(fr,d2fr,fr_1,d2fr_1)
            deallocate(FF,GG,FF_1,GG_1)
            deallocate(scatwf,FC,GC,FCP,GCP)
            call cpu_time(t2)
            write(*,*) 'cputime:',t2-t1
     
      end program
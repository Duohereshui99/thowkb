ccccccc
      module variable                           !!数组分配内存(静态)的变量
            use mesh
            use system
            implicit none
            real*8,allocatable::psi(:,:),d2psi(:,:)  !basis及其二阶导数
            real*8,allocatable::psi_1(:),d2psi_1(:)  !插值计算的中间量
            real*8,allocatable::vpot(:),vpot_1(:)                 !势
            real*8,allocatable::H(:,:)        !Hamiltonian矩阵
            real*8,allocatable::rr(:),rrw(:)               !积分点和权重 
            real*8,allocatable::wr(:),wi(:)              !特征值的实部和虚部
            real*8,allocatable::vr(:,:)       !每个特征值对应的特征向量
            real*8,allocatable::PHI(:)                  !最终的波函数(基态bound)
            complex*16,allocatable::z(:)              !特征值复数形式
            real*8,allocatable::fr(:)                  !均匀格点f(r)=1-exp(-beta*r**2)
            real*8,allocatable::d2fr(:)                !均匀格点f(r)的二阶导数
            real*8,allocatable::fr_1(:)                 !gauss积分点f(r)=1-exp(-beta*r**2)
            real*8,allocatable::d2fr_1(:)               !gauss积分点f(r)的二阶导数
            real*8,allocatable::FF(:)                  !均匀格点Coulomb函数F
            real*8,allocatable::GG(:)                  !均匀格点Coulomb函数G
            real*8,allocatable::FF_1(:)                 !gauss积分点Coulomb函数F
            real*8,allocatable::GG_1(:)                 !gauss积分点Coulomb函数G
            complex*16,allocatable::scatwf(:)           !散射波函数  
            real*8,allocatable::FC(:),GC(:)
            real*8,allocatable::FCP(:),GCP(:)      !COUL函数F,G及其一阶导数   
            real*8,allocatable::fff(:)   !!f=k(r)=2mu/h^2(Q-V(r))=k^2
            integer,allocatable::rrr(:)    !index of r_1,r_2,r_3, ...
      end module
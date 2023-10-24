cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!下调用LAPACK中的dgeev函数来计算实矩阵的特征值和特征向量
!输入量n为方阵的dimension,wr和wi分别是特征值的实部和虚部(一般为实数)
!vr为特征向量方阵
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine eigenvalue(n,a,wr,wi,vr)
            implicit none
            integer::n
            real*8::wr(n),wi(n),a(n,n),vr(n,n)
            real*8,allocatable::work(:)
            integer::lwork,info      
            allocate(work(3*n))     !work为工作区,lwork为工作区大小
            call dgeev('N','V',n,a,n,wr,wi,1,1,vr,n,work,-1,info)
            !在上面第一次call的时候令lwork=-1,
            !此时需要的工作区大小lwork存储在work(1)处
            !从而重新分配工作区大小
            lwork=work(1)
            deallocate(work)
            allocate(work(lwork))
            call dgeev('N','V',n,a,n,wr,wi,1,1,vr,n,work,lwork,info)
            deallocate(work)
      end subroutine
ccccccc
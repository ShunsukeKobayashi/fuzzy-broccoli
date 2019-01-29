module m_Test
  use m_Thomas
contains
  subroutine testIPF()

    implicit none
    real(8),allocatable :: dataA(:,:),dataB(:,:)
    integer :: i,j
    allocate(dataA(2,2),dataB(2,2))
    do i=1,2
      do j=1,2
      dataB(i,j)=i+2*(j-1)
      dataA(i,j)=2*(i-1)+j
      enddo
    enddo
    write(*,*) dataA,"dataA"
    write(*,*) dataB,"dataB"
    call InverseMatrix(dataA,dataB)
    write(*,*) dataB,"estimated B'"
    write(*,*) MATMUL(dataA,dataB)
    return
  end subroutine testIPF

  subroutine testTDM()
      implicit none
      real(8),allocatable :: A(:,:,:),B(:,:,:),C(:,:,:),R(:,:,:)
      integer :: i,j,k,n,m
      n=2
      m=3
      allocate(A(n,n,m),B(n,n,m),C(n,n,m),R(n,1,m))
      do k=1,m
        do i=1,n
          do j=1,n
            A(i,j,k) = (2*(i-1)+j)*k/2.0
            B(i,j,k) = (i+2*(j-1))*k/3.0
            C(i,j,k) = (2*k*(i-1)+j)/4.0
            R(i,1,k) = i*j*k
          enddo
        enddo
      enddo
      B(:,:,m)=0
      C(:,:,1)=0
      write(*,*) "A"
      do i=1,m
        write(*,*) A(:,:,i)
      enddo
      write(*,*) "B"
      do i=1,m
        write(*,*) B(:,:,i)
      enddo
      write(*,*) "C"
      do i=1,m
        write(*,*) C(:,:,i)
      enddo
      write(*,*) "R"
      do i=1,m
        write(*,*) R(:,:,i)
      enddo
      call n_Thomas(C,A,B,R)
      write(*,*) "A"
      do i=1,m
        write(*,*) A(:,:,i)
      enddo
      write(*,*) "B"
      do i=1,m
        write(*,*) B(:,:,i)
      enddo
      write(*,*) "C"
      do i=1,m
        write(*,*) C(:,:,i)
      enddo
      write(*,*) "X"
      do i=1,m
        write(*,*) R(:,:,i)
      enddo
      return
  end subroutine testTDM

  subroutine meshtest(mesh)
    use m_meshObj
    type(meshObj) :: mesh
    write(*,*) 'meshtest'
    write(*,*) 'mesh%pw',mesh%pw(0:1)
    write(*,*) 'mesh%pg',mesh%pg(0:1)
    write(*,*) 'mesh%C',mesh%C(0:1)
    write(*,*) 'mesh%rhow',mesh%rhow(0:1)
    write(*,*) 'mesh%rhog',mesh%rhog(0:1)
    write(*,*) 'mesh%phiw',mesh%phiw(0:1)
    write(*,*) 'mesh%phig',mesh%phig(0:1)
    write(*,*) 'mesh%Sw',mesh%Sw(0:1)
    write(*,*) 'mesh%Sg',mesh%Sg(0:1)
    write(*,*) 'mesh%chiw',mesh%chiw(0:1)
    write(*,*) 'mesh%chig',mesh%chig(0:1)
    write(*,*) 'mesh%krw',mesh%krw(0:1)
    write(*,*) 'mesh%krg',mesh%krg(0:1)
    return
  end subroutine

end module m_Test

module m_Thomas
  implicit none
  contains

    subroutine InnerProductForm(a)
      implicit none
      real(8),intent(inout) :: a(:,:)
      real(8) :: dajk
      integer :: i,j,k,n

      n=size(a(1,:))
      do k=1,n
        do j=1,k-1
          dajk = a(j,k)
          do i=j+1,n
            a(i,k) = a(i,k)-a(i,j)*dajk
          enddo
        enddo
        a(k,k) = 1.0d0/a(k,k)
        do i=k+1,n
          a(i,k) = a(i,k)*a(k,k)
        enddo
      enddo
      return
    end subroutine InnerProductForm

    subroutine Substitution(a,b)
      implicit none
      real(8),intent(inout) :: a(:,:),b(:,:)
      real(8),allocatable :: l(:,:),u(:,:),y(:,:)
      integer ::  i,j,k,n,m
      n = size(a(1,:))
      m = size(b(1,:))
      allocate(l(n,n),u(n,n),y(n,n))
      !L,Uを生成
      do j = 1,n
        do i = 1,j-1
          u(i,j) = a(i,j)
        enddo
        u(j,j) = 1.0d0/a(j,j)
        l(j,j) = 1
        do i = j+1,n
          l(i,j) = a(i,j)
        enddo
      enddo
      !ynを導出
      do k = 1,m
        y(1,k) = b(1,k)
        do i = 2,n
            y(i,k) = b(i,k)
            do j=1,i-1
              y(i,k) = y(i,k)-l(i,j)*y(j,k)
            enddo
        enddo
      enddo
      !xnを導出 (=b')
      do k = 1,m
        b(n,k) = y(n,k)/u(n,n)
        do i = n-1,1,-1
          b(i,k) = y(i,k)
          do j = i+1,n
            b(i,k) = b(i,k)-u(i,j)*b(j,k)
          enddo
          b(i,k) = b(i,k)/u(i,i)
        enddo
      enddo
      !write(*,*) b,"b"
      return
    end subroutine Substitution

    subroutine InverseMatrix(a,b)
      !B'=A^-1*Bを返す
      implicit none
      real(8),intent(in) :: a(:,:)
      real(8),intent(inout) :: b(:,:)
      real(8),allocatable :: IPF(:,:)
      integer :: n
      n = size(a(1,:))
      IPF = a
      call InnerProductForm(IPF)
      call Substitution(IPF,b)
      return
    end subroutine InverseMatrix

    subroutine n_Thomas(c,a,b,r)
      implicit none
      real(8),intent(inout) :: c(:,:,:),a(:,:,:),b(:,:,:),r(:,:,:)
      !row,column,index
      integer :: i,n

      n = size(r(1,1,:))
      call InverseMatrix(a(:,:,1),b(:,:,1))
      do i = 2,n
        a(:,:,i) = a(:,:,i)-MATMUL(c(:,:,i),b(:,:,i-1))
        call InverseMatrix(a(:,:,i),b(:,:,i))
      enddo
      call InverseMatrix(a(:,:,1),r(:,:,1))
      do i = 2,n
        r(:,:,i) = r(:,:,i) - MATMUL(c(:,:,i),r(:,:,i-1))
        call InverseMatrix(a(:,:,i),r(:,:,i))
      enddo
      do i = n-1,1,-1
        r(:,:,i) = r(:,:,i)- MATMUL(b(:,:,i),r(:,:,i+1))
      enddo
      return
    end subroutine n_Thomas

    subroutine Thomas(c,a,b,r)
      implicit none
      real(8),intent(inout) :: c(:),a(:),b(:),r(:)
      integer :: i,n

      n = size(r)
      b(1) = b(1)/a(1)
      do i = 2,n
        a(i) = a(i)-c(i)*b(i-1)
        b(i) = b(i)/a(i)
      enddo
      r(1) = r(1)/a(1)
      do i = 2,n
        r(i) = (r(i)-c(i)*r(i-1))/a(i)
      enddo
      do i = n-1,1,-1
        r(i) = r(i)-b(i)*r(i+1)
      enddo
      return
    end subroutine Thomas

end module m_Thomas

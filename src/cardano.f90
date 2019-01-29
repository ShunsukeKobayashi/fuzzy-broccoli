module m_cardano
  implicit none
contains
  subroutine cardano(a,b,c,x)
    !3次方程式をカルダノの公式により解く
    !x**3 + a*x**2 + b*x + c = 0
    implicit none
    integer :: i
    complex(8),intent(in) :: a,b,c
    complex(8),parameter :: omega = (-1d0 + sqrt(cmplx(-3d0)))/2d0
    complex(8) :: u,v,p,q,d
    complex(8),intent(out) :: x(:)

    p = b/3d0 - a**2/9d0
    q = c/2d0 + a**3/27d0 -a*b/6d0
    d = q**2+p**3
    u = (-cmplx(q) + sqrt(cmplx(d)))**(1d0/3d0)

    if(abs(u) /= 0.0d0) then
      !uv = -p
      v = -cmplx(p)/u

      x(1) = u + v -cmplx(a)/3.0d0
      x(2) = u*omega + v*omega**2 -cmplx(a)/3.0d0
      x(3) = u*omega**2 + v*omega -cmplx(a)/3.0d0
    else
      v = (-cmplx(2*q))**(1.0d0/3.0d0)

      x(1) = v -cmplx(a)/3.0d0
      x(2) = v*omega -cmplx(a)/3.0d0
      x(3) = v*omega**2 -cmplx(a)/3.0d0
    endif
  end subroutine cardano
end module m_cardano

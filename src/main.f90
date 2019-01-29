program diffusion
  use m_Test !for debug
  use m_inputCtrl
  use m_inputMesh
  use m_inputMediumProperty
  use m_inputFluidProperty
  use m_Thomas
  use m_residual
  use m_meshObj
  implicit none

  real(8),allocatable :: c(:,:,:),a(:,:,:),b(:,:,:),r(:,:,:)
  type(meshObj),allocatable :: mesh(:)
  type(inputMediumPropertyObj) :: mediumProps
  type(inputFluidPropertyObj) :: fluidProps
  type(inputCtrlObj) :: ctrl
  type(inputMeshObj) :: meshes
  real(8) :: time,tmpMaxval,injection,reductionWeight,flux
  integer :: i,n,nunit,newton
  character :: linebuf*256

  call mediumProps%InitialSet(1)
  call fluidProps%InitialSet(2)
  call meshes%InitialSet(3)
  call ctrl%InitialSet(4)

    n = meshes%numberOfMeshes
    allocate(c(3,3,n),a(3,3,n),b(3,3,n),r(3,1,n),mesh(n))

    call initializeMeshes(mesh,mediumProps,fluidProps,ctrl,meshes%mesh)
    close(nunit)
    write(*,*) 'initialization'
    write(*,*) mesh(1)%rhog(0)
    ctrl%initPw = meshes%mesh(1)%waterPressure
    ctrl%initPg = meshes%mesh(1)%gasPressure
    ctrl%initC = meshes%mesh(1)%saltConcentration
    time = 0
    do while(time < ctrl%period)
      injection = 1d2
      time = time + ctrl%interval
      newton = 0
      do while(newton < ctrl%maxNewton)
        newton = newton+1
        call fillMatrix(c,mesh,mediumProps,fluidProps,ctrl,-1)
        call fillMatrix(a,mesh,mediumProps,fluidProps,ctrl,0)
        call fillMatrix(b,mesh,mediumProps,fluidProps,ctrl,1)
        call fillMatrix(r,mesh,mediumProps,fluidProps,ctrl)
        call fillMatrixBoundary(a,b,c,r,mesh,mediumProps,fluidProps,ctrl,injection)
        call n_Thomas(c,a,b,r)
        !ダンピング比を0~1で黄金分割方　ダンピング比*更新値
        !ダンピングした値で更新してみて、質量保存を計算する
        reductionWeight = 0
        tmpMaxval = max(maxval(r(:,:,:)),abs(minval(r(:,:,:))))
        call reduction(mesh,r,reductionWeight,injection)
        r(:,:,:) = reductionWeight*r(:,:,:)
        !if(tmpMaxval > ctrl%reductionthreshold)then
        !  r(:,:,:) = ctrl%reductionthreshold*r(:,:,:)/tmpMaxval
        !endif
        !write(*,*) reductionWeight,tmpMaxval
        call updateMeshesPressure(mesh,r)
        !write(*,*) 'time',time,':',newton,injection
        if(tmpMaxval < 1d2) then
          if(injection < 1d2)then
            injection = min(1d2,injection*10)
          else
            exit
          endif
        endif
      enddo
      call updateMeshes(mesh,mediumProps,fluidProps,ctrl)
      write(*,*) 'time:',time,newton
    enddo

    open(newunit = nunit, file = "./output/output_diffusion.csv")
    write(nunit,*) 'r,pw,pg,C,sg,R^2/t'
    do i=1,n
      write(linebuf,*) real(mesh(i)%x),',',mesh(i)%pw(0),',',mesh(i)%pg(0),',',mesh(i)%C(0),',',mesh(i)%sg(0),&
      & ',',real(mesh(i)%x)**2/time
      call del_spaces(linebuf)           ! 余分な空白を削除する
      write (nunit, '(a)') trim(linebuf)    ! 出力する
    enddo
    close(nunit)
  !test
  !call testIPF()
  !call testTDM()
  !end test
  stop 'stop!'
contains

subroutine reduction(mesh,r,reductionWeight,injection)
!黄金分割法により最適なダンピング比を求める
implicit none
type(meshObj),intent(in) :: mesh(:)
double precision,intent(in) :: r(:,:,:),injection
double precision,intent(inout) :: reductionWeight
double precision,parameter :: golden_ratio = (5**0.5-1)/2
double precision :: sum_a,sum_b,a,b,lambda_a,lambda_b
double precision :: i,n
a = 0;b = 1
do while(.TRUE.)
  lambda_a = golden_ratio*(b-a) + a
  lambda_b = golden_ratio**2*(b-a) + a
  sum_a = get_sum(mesh,r,lambda_a,injection)
  sum_b = get_sum(mesh,r,lambda_b,injection)
  if(sum_a < sum_b)then
    a = lambda_a
    if(b-lambda_a < 1d-4)then
      reductionWeight = (b+lambda_a)/2
      exit
    endif
  else
    b = lambda_b
    if(lambda_b-a < 1d-4)then
      reductionWeight =(a+lambda_b)/2
      exit
    endif
  endif
enddo
end subroutine

double precision function get_sum(mesh,r,reductionWeight,injection)
  implicit none
  type(meshObj),intent(in) :: mesh(:)
  double precision,intent(in) :: r(:,:,:),injection
  double precision,intent(inout) :: reductionWeight
  type(meshObj),allocatable :: tmpMesh(:)
  double precision,allocatable :: w_residual(:),g_residual(:),c_residual(:)
  integer :: i,n
  n = size(mesh)
  allocate(tmpMesh(n),w_residual(n),g_residual(n),c_residual(n))
  !tmpR(1:n,1:n,1:n) = r(1:n,1:n,1:n)*reductionWeight
  tmpMesh(1:n) = mesh(1:n)
  call updateMeshesPressure(tmpMesh,r*reductionWeight)
  w_residual(1) = boundaryWaterResidual(tmpMesh(1),tmpMesh(2),mediumProps,fluidProps,ctrl)
  g_residual(1) = boundaryGasResidual(injection,tmpMesh(1),tmpMesh(2),mediumProps,fluidProps,ctrl)
  c_residual(1) = 0
  do i=2,n-1
    w_residual(i) = waterResidual(tmpMesh(i-1),tmpMesh(i),tmpMesh(i+1),mediumProps,fluidProps,ctrl)
    g_residual(i) = gasResidual(tmpMesh(i-1),tmpMesh(i),tmpMesh(i+1),mediumProps,fluidProps,ctrl)
    c_residual(i) = saltResidual(tmpMesh(i-1),tmpMesh(i),tmpMesh(i+1),mediumProps,fluidProps,ctrl)
  enddo
  get_sum = sum(abs(w_residual)) + sum(abs(g_residual)) + sum(abs(c_residual))
  deallocate(tmpMesh,w_residual,g_residual,c_residual)
  return
end function get_sum

subroutine del_spaces(s)
  character (*), intent (inout) :: s
  character (len=len(s)) tmp
  integer i, j
  j = 1
  do i = 1, len(s)
    if (s(i:i)==' ') cycle
    tmp(j:j) = s(i:i)
    j = j + 1
  end do
  s = tmp(1:j-1)
end subroutine del_spaces

end program diffusion

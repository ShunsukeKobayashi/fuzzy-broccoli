module m_meshObj
  use m_inputCtrl
  use m_inputMesh
  use m_inputMediumProperty
  use m_inputFluidProperty
  implicit none
  type meshObj
    real(8) :: x,dx,pw(0:1),pg(0:1),C(0:1),rhow(0:1),rhog(0:1),phiw(0:1),phig(0:1)&
    & ,Sw(0:1),Sg(0:1),chiw(0:1),chig(0:1),krw(0:1),krg(0:1)
  end type meshObj
contains
  complex(8) function getPw(mesh,mediumProps,fluidProps,ctrl,kind)
    implicit none
    type(meshObj),intent(in) :: mesh
    type(inputMediumPropertyObj) :: mediumProps
    type(inputFluidPropertyObj) :: fluidProps
    type(inputCtrlObj) :: ctrl
    integer,optional,intent(in) :: kind
    complex(8),parameter :: ii = (0,1d0)
    complex(8) :: pw
    if(present(kind))then
      if(kind == 0)then
        getPw = mesh%pw(1) + ii*ctrl%epsilon
      else
        getPw = mesh%pw(1)
      endif
    else
      getPw = mesh%pw(1)
    endif
    return
  end function getPw

  complex(8) function getPg(mesh,mediumProps,fluidProps,ctrl,kind)
    implicit none
    type(meshObj),intent(in) :: mesh
    type(inputMediumPropertyObj) :: mediumProps
    type(inputFluidPropertyObj) :: fluidProps
    type(inputCtrlObj) :: ctrl
    integer,optional,intent(in) :: kind
    complex(8),parameter :: ii = (0,1d0)
    complex(8) :: pg
    if(present(kind))then
      if(kind == 1)then
        getPg = mesh%pg(1) + ii*ctrl%epsilon
      else
        getPg = mesh%pg(1)
      endif
    else
      getPg = mesh%pg(1)
    endif
    return
  end function getPg

  complex(8) function getC(mesh,mediumProps,fluidProps,ctrl,kind)
    implicit none
    type(meshObj),intent(in) :: mesh
    type(inputMediumPropertyObj) :: mediumProps
    type(inputFluidPropertyObj) :: fluidProps
    type(inputCtrlObj) :: ctrl
    integer,optional,intent(in) :: kind
    complex(8),parameter :: ii = (0,1d0)
    if(present(kind))then
      if(kind == 2)then
        getC = mesh%C(1) + ii*ctrl%epsilon
      else
        getC = mesh%C(1)
      endif
    else
      getC = mesh%C(1)
    endif
    return
  end function getC

  complex(8) function getOsmoticPressure(mesh,mediumProps,fluidProps,ctrl,kind)
    use m_physicalConstant
    implicit none
    type(meshObj),intent(in) :: mesh
    type(inputMediumPropertyObj) :: mediumProps
    type(inputFluidPropertyObj) :: fluidProps
    type(inputCtrlObj) :: ctrl
    integer,optional,intent(in) :: kind
    complex(8) :: rhow,C
    rhow = getRhoW(mesh,mediumProps,fluidProps,ctrl,kind)
    C = getC(mesh,mediumProps,fluidProps,ctrl,kind)
    getOsmoticPressure = rhow*C*R*fluidProps%temperature/fluidProps%molarConcentration
    return
  end function getOsmoticPressure

  complex(8) function getRhoG(mesh,mediumProps,fluidProps,ctrl,kind)
    use m_physicalConstant
    implicit none
    type(meshObj),intent(in) :: mesh
    type(inputMediumPropertyObj) :: mediumProps
    type(inputFluidPropertyObj) :: fluidProps
    type(inputCtrlObj) :: ctrl
    integer,intent(in),optional :: kind
    complex(8) :: pg,v
    ! pg = getPg(mesh,mediumProps,fluidProps,ctrl,kind)
    ! getRhoG = fluidProps%molecularWeight/(R*fluidProps%temperature) * pg
    !redlich-kwong equation, TODO 複素数階微分に対応する
    v = getVolume(mesh,mediumProps,fluidProps,ctrl,kind)
    getRhoG = fluidProps%molecularWeight/v
    return
  end function getRhoG

  complex(8) function getVolume(mesh,mediumProps,fluidProps,ctrl,kind)
    !m^3/molで分子容を計算している。
    use m_physicalConstant
    use m_molecularProps
    use m_cardano
    implicit none
    type(meshObj),intent(in) :: mesh
    type(inputMediumPropertyObj) :: mediumProps
    type(inputFluidPropertyObj) :: fluidProps
    type(inputCtrlObj) :: ctrl
    integer,intent(in),optional :: kind
    complex(8) :: s,t,u,a,pg
    double precision :: vliquid,vgas,w1,w2
    integer :: i
    complex(8),allocatable :: x(:)
    allocate(x(3))
    x(:) = 0
    a = getAg(fluidProps)
    pg = getPg(mesh,mediumProps,fluidProps,ctrl,kind)
    s = - R*fluidProps%temperature/pg
    t = - (R*fluidProps%temperature*bg/pg - a/(pg*sqrt(fluidProps%temperature)) + bg**2)
    u = - a*bg/(pg*sqrt(fluidProps%temperature))
    call cardano(s,t,u,x)
    if(imag(x(2)+x(3)) /= 0)then
      !3つの実数解を持つ
      vliquid = minval(real(x))
      vgas = maxval(real(x))
      w1 = mesh%pg(1)*(vgas-vliquid)
      w2 = R*fluidProps%temperature*log((vgas-bg)/(vliquid-bg)) + &
      & a/(sqrt(fluidProps%temperature))*bg*log(((vgas+bg)*vliquid)/((vliquid+bg)*vgas))
      !write(*,*) w1,w2
      if(w2 > w1)then
        getVolume = vgas
      else
        getVolume = vliquid
      endif
    else
      !1つだけ実数解を持つ
      !write(*,*) 'kyosu'
      getVolume = real(x(1))
    endif
    deallocate(x)
    return
  end function getVolume

  double precision function getAg(fluidProps)
    !Redlich-Kwong式のパラメータa_mixを計算する、280-380Kでフィッティングする
    implicit none
    type(inputFluidPropertyObj),intent(in) :: fluidProps
    getAg = 7.54 - 4.02*1d-3*fluidProps%temperature
    return
  end function getAg

  complex(8) function getRhoW(mesh,mediumProps,fluidProps,ctrl,kind)
    implicit none
    type(meshObj),intent(in) :: mesh
    type(inputMediumPropertyObj) :: mediumProps
    type(inputFluidPropertyObj) :: fluidProps
    type(inputCtrlObj) :: ctrl
    integer,intent(in),optional :: kind
    complex(8) :: pw,C
    pw = getPw(mesh,mediumProps,fluidProps,ctrl,kind)
    C = getC(mesh,mediumProps,fluidProps,ctrl,kind)
    getRhoW = fluidProps%rhow0 + (pw - fluidProps%p0)/fluidProps%kw + fluidProps%gamma*C
    return
  end function getRhoW

  complex(8) function getSwe(mesh,mediumProps,fluidProps,ctrl,kind)
    implicit none
    type(meshObj),intent(in) :: mesh
    type(inputMediumPropertyObj) :: mediumProps
    type(inputFluidPropertyObj) :: fluidProps
    type(inputCtrlObj) :: ctrl
    integer,optional :: kind
    complex(8) pg,pw
    real(8) :: m
    m = 1d0 - 1d0/mediumProps%n
    pw = getPw(mesh,mediumProps,fluidProps,ctrl,kind)
    pg = getPg(mesh,mediumProps,fluidProps,ctrl,kind)
    getSwe = (1d0 + (mediumProps%arufa *(pg-pw))**mediumProps%n)**(-m)
    return
  end function getSwe

  complex(8) function getSge(mesh,mediumProps,fluidProps,ctrl,kind)
    implicit none
    type(meshObj),intent(in) :: mesh
    type(inputMediumPropertyObj) :: mediumProps
    type(inputFluidPropertyObj) :: fluidProps
    type(inputCtrlObj) :: ctrl
    integer,optional :: kind
    complex(8) :: Sw
    Sw = getSw(mesh,mediumProps,fluidProps,ctrl,kind)
    getSge = ((1d0-Sw) - fluidProps%gasSaturateRate)/(1d0-fluidProps%gasSaturateRate)
    return
  end function getSge

  complex(8) function getSw(mesh,mediumProps,fluidProps,ctrl,kind)
    implicit none
    type(meshObj),intent(in) :: mesh
    type(inputMediumPropertyObj) :: mediumProps
    type(inputFluidPropertyObj) :: fluidProps
    type(inputCtrlObj) :: ctrl
    integer,optional :: kind
    complex(8) :: Swe
    Swe = getSwe(mesh,mediumProps,fluidProps,ctrl,kind)
    getSw = Swe*(1-fluidProps%waterSaturateRate) + fluidProps%waterSaturateRate
    return
  end function getSw

  complex(8) function getKrW(mesh,mediumProps,fluidProps,ctrl,kind)
    implicit none
    type(meshObj),intent(in) :: mesh
    type(inputMediumPropertyObj) :: mediumProps
    type(inputFluidPropertyObj) :: fluidProps
    type(inputCtrlObj) :: ctrl
    integer,optional :: kind
    complex(8) :: pw,pg
    real(8) :: m
    m = 1d0 - 1d0/mediumProps%n
    pw = getPw(mesh,mediumProps,fluidProps,ctrl,kind)
    pg = getPg(mesh,mediumProps,fluidProps,ctrl,kind)
    getKrW = (1d0 + (mediumProps%arufa * (pg-pw))**mediumProps%n)**(-m*mediumProps%L) * &
    & (1d0 - (1d0 - (1d0 + (mediumProps%arufa*(pg-pw))**mediumProps%n)**(-1d0))**m)**2
    return
  end function getKrW

  complex(8) function getKrG(mesh,mediumProps,fluidProps,ctrl,kind)
    implicit none
    type(meshObj),intent(in) :: mesh
    type(inputMediumPropertyObj) :: mediumProps
    type(inputFluidPropertyObj) :: fluidProps
    type(inputCtrlObj) :: ctrl
    integer,optional :: kind
    complex(8) :: Sge,Sw,s_hat
    real(8) :: m
    !m = 1d0 - 1d0/mediumProps%n
    !Sge = getSge(mesh,mediumProps,fluidProps,ctrl,kind)
    !getKrG = Sge**mediumProps%L*(1d0 - (1d0 - Sge**(1d0/m))**m)**2
    !getKrG = (1d0 - Swe)**2*(1d0 -Swe**(5d0/3d0)) !どちらでも
    sw = getSw(mesh,mediumProps,fluidProps,ctrl,kind)
    S_hat = (sw - 0.30)/0.65
    getKrG = (1-s_hat**2)*(1-s_hat)**2
    return
  end function getKrG

  complex(8) function getPhiW(mesh,mediumProps,fluidProps,ctrl,kind)
    implicit none
    type(meshObj),intent(in) :: mesh
    type(inputMediumPropertyObj) :: mediumProps
    type(inputFluidPropertyObj) :: fluidProps
    type(inputCtrlObj) :: ctrl
    integer,optional :: kind
    complex(8) :: Swe,Sw,pw,pg
    Swe = getSwe(mesh,mediumProps,fluidProps,ctrl,kind)
    Sw = getSw(mesh,mediumProps,fluidProps,ctrl,kind)
    pw = getPw(mesh,mediumProps,fluidProps,ctrl,kind)
    pg = getPg(mesh,mediumProps,fluidProps,ctrl,kind)
    getPhiW = mesh%phiw(0) + (mesh%phiw(0)+mesh%phig(0))*(Sw - mesh%Sw(0)) + mediumProps%compressibility*Swe* &
    & (Swe*pw - mesh%chiw(0)*mesh%pw(0) + (1d0-Swe)*pg - mesh%chig(0)*mesh%pg(0))
    return
  end function getPhiW

  complex(8) function getPhiG(mesh,mediumProps,fluidProps,ctrl,kind)
    implicit none
    type(meshObj),intent(in) :: mesh
    type(inputMediumPropertyObj) :: mediumProps
    type(inputFluidPropertyObj) :: fluidProps
    type(inputCtrlObj) :: ctrl
    integer,optional :: kind
    complex(8) :: Swe,Sw,pw,pg
    Swe = getSwe(mesh,mediumProps,fluidProps,ctrl,kind)
    Sw = getSw(mesh,mediumProps,fluidProps,ctrl,kind)
    pw = getPw(mesh,mediumProps,fluidProps,ctrl,kind)
    pg = getPg(mesh,mediumProps,fluidProps,ctrl,kind)
    getPhiG = mesh%phig(0) + (mesh%phiw(0)+mesh%phig(0))*((1d0-Sw) - mesh%Sg(0)) &
    & + mediumProps%compressibility*(1d0-Swe)*(Swe*pw - mesh%chiw(0)*mesh%pw(0) + (1d0-Swe)*pg - mesh%chig(0)*mesh%pg(0))
    return
  end function getPhiG

  subroutine initializeMeshes(mesh,mediumProps,fluidProps,ctrl,inputParams)
    implicit none
    integer :: i,n
    type(meshObj) :: mesh(:)
    type(inputMediumPropertyObj) :: mediumProps
    type(inputFluidPropertyObj) :: fluidProps
    type(inputCtrlObj) :: ctrl
    type(inputMeshParam) :: inputParams(:)

    n = size(mesh)
    do i = 1,n
      mesh(i)%x=inputParams(i)%coordinate
      mesh(i)%dx=inputParams(i)%meshLength
      mesh(i)%pw(0)=inputParams(i)%waterPressure
      mesh(i)%pg(0)=inputParams(i)%gasPressure
      mesh(i)%C(0)=inputParams(i)%saltConcentration
      mesh(i)%pw(1)=mesh(i)%pw(0)
      mesh(i)%pg(1)=mesh(i)%pg(0)
      mesh(i)%C(1) = mesh(i)%C(0)

      mesh(i)%rhow(0) = getRhoW(mesh(i),mediumProps,fluidProps,ctrl)
      mesh(i)%rhog(0) = getRhoG(mesh(i),mediumProps,fluidProps,ctrl)
      mesh(i)%rhow(1)=mesh(i)%rhow(0)
      mesh(i)%rhog(1)=mesh(i)%rhog(0)
      mesh(i)%chiw(0) = getSwe(mesh(i),mediumProps,fluidProps,ctrl)
      mesh(i)%chig(0) = 1d0 - mesh(i)%chiw(0)
      mesh(i)%chiw(1)=mesh(i)%chiw(0)
      mesh(i)%chig(1)=mesh(i)%chig(0)
      mesh(i)%Sw(0) = getSw(mesh(i),mediumProps,fluidProps,ctrl)
      mesh(i)%Sg(0) = 1d0 - mesh(i)%Sw(0)
      mesh(i)%Sw(1)=mesh(i)%Sw(0)
      mesh(i)%Sg(1)=mesh(i)%Sg(0)
      mesh(i)%krw(0) = getKrW(mesh(i),mediumProps,fluidProps,ctrl)
      mesh(i)%krg(0) = getKrG(mesh(i),mediumProps,fluidProps,ctrl)
      mesh(i)%krw(1)=mesh(i)%krw(0)
      mesh(i)%krg(1)=mesh(i)%krg(0)
      mesh(i)%phiw(0) = mediumProps%initPhiW
      mesh(i)%phig(0) = mediumProps%initPhiG
      mesh(i)%phiw(1)=mesh(i)%phiw(0)
      mesh(i)%phig(1)=mesh(i)%phig(0)
    enddo
    return
  end subroutine initializeMeshes

  subroutine updateMeshesPressure(mesh,array)
    implicit none
    integer :: n
    real(8) :: array(:,:,:)
    type(meshObj) :: mesh(:)
    n = size(mesh)
    mesh(1:n)%pw(1) = array(1,1,1:n) + mesh(1:n)%pw(1)
    mesh(1:n)%pg(1) = array(2,1,1:n) + mesh(1:n)%pg(1)
    mesh(1:n)%C(1) = array(3,1,1:n) + mesh(1:n)%C(1)
  end subroutine updateMeshesPressure

  subroutine updateMeshes(mesh,mediumProps,fluidProps,ctrl)
    implicit none
    integer :: i,n
    type(meshObj) :: mesh(:)
    type(inputMediumPropertyObj) :: mediumProps
    type(inputFluidPropertyObj) :: fluidProps
    type(inputCtrlObj) :: ctrl

    n = size(mesh)
    do i = 1,n
      mesh(i)%pw(0)=mesh(i)%pw(1)
      mesh(i)%pg(0)=mesh(i)%pg(1)
      mesh(i)%C(0)=mesh(i)%C(1)

      mesh(i)%rhow(0)=mesh(i)%rhow(1)
      mesh(i)%rhog(0)=mesh(i)%rhog(1)
      mesh(i)%chiw(0)=mesh(i)%chiw(1)
      mesh(i)%chig(0)=mesh(i)%chig(1)
      mesh(i)%Sw(0)=mesh(i)%Sw(1)
      mesh(i)%Sg(0)=mesh(i)%Sg(1)
      mesh(i)%krw(0)=mesh(i)%krw(1)
      mesh(i)%krg(0)=mesh(i)%krg(1)
      mesh(i)%phiw(0)=mesh(i)%phiw(1)
      mesh(i)%phig(0)=mesh(i)%phig(1)

      mesh(i)%rhow(1) = getRhoW(mesh(i),mediumProps,fluidProps,ctrl)
      mesh(i)%rhog(1) = getRhoG(mesh(i),mediumProps,fluidProps,ctrl)
      mesh(i)%chiw(1) = getSwe(mesh(i),mediumProps,fluidProps,ctrl)
      mesh(i)%chig(1) = 1d0 - mesh(i)%chiw(1)
      mesh(i)%Sw(1) = getSw(mesh(i),mediumProps,fluidProps,ctrl)
      mesh(i)%Sg(1) = 1d0 - mesh(i)%Sw(1)
      mesh(i)%krw(1) = getKrW(mesh(i),mediumProps,fluidProps,ctrl)
      mesh(i)%krg(1) = getKrG(mesh(i),mediumProps,fluidProps,ctrl)
      mesh(i)%phiw(1) = getPhiW(mesh(i),mediumProps,fluidProps,ctrl)
      mesh(i)%phig(1) = getPhiG(mesh(i),mediumProps,fluidProps,ctrl)
    enddo
    return
  end subroutine updateMeshes
end module m_meshObj

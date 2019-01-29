module m_residual
  use m_meshObj
  implicit none
  double precision,parameter :: pi = atan(1d0)*4d0
contains
  complex(8) function boundaryWaterResidual(centerMesh,rightMesh,mediumProps,fluidProps,ctrl,index,kind)
    implicit none
    type(meshObj) :: centerMesh,rightMesh
    type(inputMediumPropertyObj) :: mediumProps
    type(inputFluidPropertyObj) :: fluidProps
    type(inputCtrlObj) :: ctrl
    integer,optional :: index,kind
    complex(8) :: rhoPhiW,massfluxW,flux
    rhoPhiW = differentialRhoPhiW(centerMesh,mediumProps,fluidProps,ctrl,index,kind)
    flux = getMassFluxW(centerMesh,rightMesh,mediumProps,fluidProps,ctrl,index,kind,1)/centerMesh%dx
    massfluxW = 2*pi*100*((centerMesh%x+centerMesh%dx)*flux - centerMesh%x*0)
    boundaryWaterResidual = pi*100*((centerMesh%x+centerMesh%dx)**2 - centerMesh%x**2)*rhoPhiW + massfluxW
    return
  end function boundaryWaterResidual

  complex(8) function boundaryGasResidual(injection,centerMesh,rightMesh,mediumProps,fluidProps,ctrl,index,kind)
    implicit none
    type(meshObj) :: centerMesh,rightMesh
    type(inputMediumPropertyObj) :: mediumProps
    type(inputFluidPropertyObj) :: fluidProps
    type(inputCtrlObj) :: ctrl
    integer,optional :: index,kind
    !index -1:left,0:center,1:right
    !kind 0:water,1:gas,2:salt
    double precision :: injection
    complex(8) :: rhoPhiG,massfluxG,flux
    rhoPhiG = differentialRhoPhiG(centerMesh,mediumProps,fluidProps,ctrl,index,kind)
    flux = getMassFluxG(centerMesh,rightMesh,mediumProps,fluidProps,ctrl,index,kind,1)
    massfluxG = 2*pi*100*((centerMesh%x+centerMesh%dx)*flux) - injection
    boundaryGasResidual = pi*100*((centerMesh%x+centerMesh%dx)**2 - centerMesh%x**2)*rhoPhiG + massfluxG
    return
  end function boundaryGasResidual

  complex(8) function boundarySaltResidual(centerMesh,rightMesh,mediumProps,fluidProps,ctrl,index,kind)
    implicit none
    type(meshObj) :: centerMesh,rightMesh
    type(inputMediumPropertyObj) :: mediumProps
    type(inputFluidPropertyObj) :: fluidProps
    type(inputCtrlObj) :: ctrl
    integer,optional :: index,kind
    !index -1:left,0:center,1:right
    !kind 0:water,1:gas,2:salt
    complex(8) :: differentialC,massfluxC,C,rho,flux
    if(present(index))then
      if(index == 0)then
        rho = getRhoW(centerMesh,mediumProps,fluidProps,ctrl,kind)
        C = getC(centerMesh,mediumProps,fluidProps,ctrl,kind)
      else
        rho = getRhoW(centerMesh,mediumProps,fluidProps,ctrl)
        C = getC(centerMesh,mediumProps,fluidProps,ctrl)
      endif
    else
      rho = getRhoW(centerMesh,mediumProps,fluidProps,ctrl)
      C = getC(centerMesh,mediumProps,fluidProps,ctrl)
    endif
    differentialC = (C-centerMesh%c(0))/ctrl%interval
    flux = getMassFluxC(centerMesh,rightMesh,mediumProps,fluidProps,ctrl,index,kind,1)
    massfluxC = 2*pi*100*((centerMesh%x+centerMesh%dx)*flux)
    boundarySaltResidual = pi*100*((centerMesh%x+centerMesh%dx)**2 - centerMesh%x**2)*rho*differentialC + massfluxC
    return
  end function boundarySaltResidual

  complex(8) function waterResidual(leftMesh,centerMesh,rightMesh,mediumProps,fluidProps,ctrl,index,kind)
    implicit none
    type(meshObj) :: leftMesh,centerMesh,rightMesh
    type(inputMediumPropertyObj) :: mediumProps
    type(inputFluidPropertyObj) :: fluidProps
    type(inputCtrlObj) :: ctrl
    integer,optional :: index,kind
    !index -1:left,0:center,1:right
    !kind 0:water,1:gas,2:salt
    complex(8) :: rhoPhiW,massfluxW
    rhoPhiW = differentialRhoPhiW(centerMesh,mediumProps,fluidProps,ctrl,index,kind)
    massfluxW = differentialMassFluxW(leftMesh,centerMesh,rightMesh,mediumProps,fluidProps,ctrl,index,kind)
    waterResidual = pi*100*((centerMesh%x+centerMesh%dx)**2 - centerMesh%x**2)*rhoPhiW + massfluxW
    return
  end function waterResidual

  complex(8) function gasResidual(leftMesh,centerMesh,rightMesh,mediumProps,fluidProps,ctrl,index,kind)
    implicit none
    type(meshObj) :: leftMesh,centerMesh,rightMesh
    type(inputMediumPropertyObj) :: mediumProps
    type(inputFluidPropertyObj) :: fluidProps
    type(inputCtrlObj) :: ctrl
    integer,optional :: index,kind
    !index -1:left,0:center,1:right
    !kind 0:water,1:gas,2:salt
    complex(8) :: rhoPhiG,massfluxG
    rhoPhiG = differentialRhoPhiG(centerMesh,mediumProps,fluidProps,ctrl,index,kind)
    massfluxG = differentialMassFluxG(leftMesh,centerMesh,rightMesh,mediumProps,fluidProps,ctrl,index,kind)
    gasResidual = pi*100*((centerMesh%x+centerMesh%dx)**2 - centerMesh%x**2)*rhoPhiG + massfluxG
    return
  end function gasResidual

  complex(8) function saltResidual(leftMesh,centerMesh,rightMesh,mediumProps,fluidProps,ctrl,index,kind)
    implicit none
    type(meshObj) :: leftMesh,centerMesh,rightMesh
    type(inputMediumPropertyObj) :: mediumProps
    type(inputFluidPropertyObj) :: fluidProps
    type(inputCtrlObj) :: ctrl
    integer,optional :: index,kind
    !index -1:left,0:center,1:right
    !kind 0:water,1:gas,2:salt
    complex(8) :: differentialC,massfluxC,C,rho
    if(present(index))then
      if(index == 0)then
        rho = getRhoW(centerMesh,mediumProps,fluidProps,ctrl,kind)
        C = getC(centerMesh,mediumProps,fluidProps,ctrl,kind)
      else
        rho = getRhoW(centerMesh,mediumProps,fluidProps,ctrl)
        C = getC(centerMesh,mediumProps,fluidProps,ctrl)
      endif
    else
      rho = getRhoW(centerMesh,mediumProps,fluidProps,ctrl)
      C = getC(centerMesh,mediumProps,fluidProps,ctrl)
    endif
    differentialC = (C-centerMesh%c(0))/ctrl%interval
    massfluxC = differentialMassFluxC(leftMesh,centerMesh,rightMesh,mediumProps,fluidProps,ctrl,index,kind)
    saltResidual = pi*100*((centerMesh%x+centerMesh%dx)**2 - centerMesh%x**2)*rho*differentialC + massfluxC
    return
  end function saltResidual

  complex(8) function differentialRhoPhiW(mesh,mediumProps,fluidProps,ctrl,index,kind)
    !residualの第１項、水用
    implicit none
    integer,optional,intent(in) :: index,kind
    type(meshObj),intent(in) :: mesh
    type(inputMediumPropertyObj) :: mediumProps
    type(inputFluidPropertyObj) :: fluidProps
    type(inputCtrlObj) :: ctrl
    complex(8) :: rhow,rhoPhi
    !新しいタイムステップでの値,complexであることに注意
    if(present(index))then
      if(index == 0)then
        rhow = getRhoW(mesh,mediumProps,fluidProps,ctrl,kind)
        rhoPhi = rhow*getPhiW(mesh,mediumProps,fluidProps,ctrl,kind)
      else
        rhow = getRhoW(mesh,mediumProps,fluidProps,ctrl)
        rhoPhi = rhow*getPhiW(mesh,mediumProps,fluidProps,ctrl)
      endif
    else
      rhow = getRhoW(mesh,mediumProps,fluidProps,ctrl)
      rhoPhi = rhow*getPhiW(mesh,mediumProps,fluidProps,ctrl)
    endif
    differentialRhoPhiW = (rhoPhi - mesh%rhow(0)*mesh%phiw(0))/ctrl%interval
    return
  end function differentialRhoPhiW

  complex(8) function differentialRhoPhiG(mesh,mediumProps,fluidProps,ctrl,index,kind)
    implicit none
    integer,optional,intent(in) :: index,kind
    type(meshObj),intent(in) :: mesh
    type(inputMediumPropertyObj) :: mediumProps
    type(inputFluidPropertyObj) :: fluidProps
    type(inputCtrlObj) :: ctrl
    complex(8) :: rhog,rhoPhi
    !新しいタイムステップでの値,complexであることに注意
    if(present(index))then
      if(index == 0)then
        rhog = getRhoG(mesh,mediumProps,fluidProps,ctrl,kind)
        rhophi = rhog*getPhiG(mesh,mediumProps,fluidProps,ctrl,kind)
      else
        rhog = getRhoG(mesh,mediumProps,fluidProps,ctrl)
        rhoPhi = rhog*getPhiG(mesh,mediumProps,fluidProps,ctrl)
      endif
    else
      rhog = getRhoG(mesh,mediumProps,fluidProps,ctrl)
      rhoPhi = rhog*getPhiG(mesh,mediumProps,fluidProps,ctrl)
    endif
    differentialRhoPhiG = (rhoPhi - mesh%rhog(0)*mesh%phig(0))/ctrl%interval
    return
  end function differentialRhoPhiG

  complex(8) function differentialMassFluxW(leftMesh,centerMesh,rightMesh,mediumProps,fluidProps,ctrl,index,kind)
    implicit none
    integer,optional,intent(in) :: index,kind
    type(meshObj),intent(in) :: leftMesh,centerMesh,rightMesh
    type(inputMediumPropertyObj) :: mediumProps
    type(inputFluidPropertyObj) :: fluidProps
    type(inputCtrlObj) :: ctrl
    complex(8),allocatable :: massflux(:)
    allocate(massflux(0:1))
    massflux(0) = getMassFluxW(leftMesh,centerMesh,mediumProps,fluidProps,ctrl,index,kind,0)
    massflux(1) = getMassFluxW(centerMesh,rightMesh,mediumProps,fluidProps,ctrl,index,kind,1)
    differentialMassFluxW = (massflux(1)*(centerMesh%x+centerMesh%dx) - massflux(0)*centerMesh%x)*2*pi*100
    return
  end function differentialMassFluxW

  complex(8) function differentialMassFluxG(leftMesh,centerMesh,rightMesh,mediumProps,fluidProps,ctrl,index,kind)
    implicit none
    integer,optional,intent(in) :: index,kind
    type(meshObj),intent(in) :: leftMesh,centerMesh,rightMesh
    type(inputMediumPropertyObj) :: mediumProps
    type(inputFluidPropertyObj) :: fluidProps
    type(inputCtrlObj) :: ctrl
    complex(8),allocatable :: massflux(:)
    allocate(massflux(0:1))
    massflux(0) = getMassFluxG(leftMesh,centerMesh,mediumProps,fluidProps,ctrl,index,kind,0)
    massflux(1) = getMassFluxG(centerMesh,rightMesh,mediumProps,fluidProps,ctrl,index,kind,1)
    differentialMassFluxG = (massflux(1)*(centerMesh%x+centerMesh%dx) - massflux(0)*centerMesh%x)*2*pi*100
    deallocate(massflux)
    return
  end function differentialMassFluxG

  complex(8) function differentialMassFluxC(leftMesh,centerMesh,rightMesh,mediumProps,fluidProps,ctrl,index,kind)
    implicit none
    integer,optional,intent(in) :: index,kind
    type(meshObj),intent(in) :: leftMesh,centerMesh,rightMesh
    type(inputMediumPropertyObj) :: mediumProps
    type(inputFluidPropertyObj) :: fluidProps
    type(inputCtrlObj) :: ctrl
    complex(8),allocatable :: massflux(:)
    allocate(massflux(0:1))
    massflux(0) = getMassFluxC(leftMesh,centerMesh,mediumProps,fluidProps,ctrl,index,kind,0)
    massflux(1) = getMassFluxC(centerMesh,rightMesh,mediumProps,fluidProps,ctrl,index,kind,1)
    differentialMassFluxC = (massflux(1)*(centerMesh%x+centerMesh%dx) - massflux(0)*centerMesh%x)*2*pi*100
    deallocate(massflux)
    return
  end function differentialMassFluxC

  complex(8) function getMassFluxW(meshA,meshB,mediumProps,fluidProps,ctrl,index,kind,meshBIndex)
    use m_physicalConstant
    implicit none
    integer,optional,intent(in) :: index,kind,meshBIndex
    type(meshObj),intent(in) :: meshA,meshB
    type(inputMediumPropertyObj) :: mediumProps
    type(inputFluidPropertyObj) :: fluidProps
    type(inputCtrlObj) :: ctrl
    complex(8) :: krw,rhow,pwB,pwA,poA,poB
    real(8) :: innerBracket
    pwA = meshA%pw(1)
    pwB = meshB%pw(1)
    poA = getOsmoticPressure(meshA,mediumProps,fluidProps,ctrl)
    poB = getOsmoticPressure(meshB,mediumProps,fluidProps,ctrl)
    innerBracket = ((pwB-pwA) - mediumProps%reflectionCoef*(poB - poA))/meshA%dx - meshA%rhow(1)*g
    if(innerBracket > 0)then
      !meshBが上流
      if(present(index))then
        if(index == meshBIndex)then
          pwB = getPw(meshB,mediumProps,fluidProps,ctrl,kind)
          poB = getOsmoticPressure(meshB,mediumProps,fluidProps,ctrl,kind)
          rhow = getRhoW(meshB,mediumProps,fluidProps,ctrl,kind)
          krw = getKrW(meshB,mediumProps,fluidProps,ctrl,kind)
        else if(index == meshBIndex-1)then
          pwA = getPw(meshA,mediumProps,fluidProps,ctrl,kind)
          poA = getOsmoticPressure(meshA,mediumProps,fluidProps,ctrl,kind)
          rhow = getRhoW(meshB,mediumProps,fluidProps,ctrl)
          krw = getKrW(meshB,mediumProps,fluidProps,ctrl)
        else
          krw = getKrW(meshB,mediumProps,fluidProps,ctrl)
          rhow = getRhoW(meshB,mediumProps,fluidProps,ctrl)
        endif
      else
        krw = getKrW(meshB,mediumProps,fluidProps,ctrl)
        rhow = getRhoW(meshB,mediumProps,fluidProps,ctrl)
      endif
    else
      if(present(index))then
        if(index == meshBIndex-1)then
          pwA = getPw(meshA,mediumProps,fluidProps,ctrl,kind)
          poA = getOsmoticPressure(meshA,mediumProps,fluidProps,ctrl,kind)
          rhow = getRhoW(meshA,mediumProps,fluidProps,ctrl,kind)
          krw = getKrW(meshA,mediumProps,fluidProps,ctrl,kind)
        else if(index == meshBIndex)then
          pwB = getPw(meshB,mediumProps,fluidProps,ctrl,kind)
          poB = getOsmoticPressure(meshB,mediumProps,fluidProps,ctrl,kind)
          rhow = getRhoW(meshA,mediumProps,fluidProps,ctrl)
          krw = getKrW(meshA,mediumProps,fluidProps,ctrl)
        else
          rhow = getRhoW(meshA,mediumProps,fluidProps,ctrl)
          krw = getKrW(meshA,mediumProps,fluidProps,ctrl)
        endif
      else
        krw = getKrW(meshA,mediumProps,fluidProps,ctrl)
        rhow = getRhoW(meshA,mediumProps,fluidProps,ctrl)
      endif
    endif
    getMassFluxW = - mediumProps%permeability*rhow*krw/fluidProps%waterViscosity* &
    & (((pwB-pwA)- mediumProps%reflectionCoef*(poB - poA))/meshA%dx - rhow*g)
    return
  end function getMassFluxW

  complex(8) function getMassFluxG(meshA,meshB,mediumProps,fluidProps,ctrl,index,kind,meshBIndex)
    use m_physicalConstant
    implicit none
    integer,optional,intent(in) :: index,kind,meshBIndex
    type(meshObj),intent(in) :: meshA,meshB
    type(inputMediumPropertyObj) :: mediumProps
    type(inputFluidPropertyObj) :: fluidProps
    type(inputCtrlObj) :: ctrl
    complex(8) :: krg,rhog,pgA,pgB
    real(8) :: innerBracket
    innerBracket = (meshB%pg(1)-meshA%pg(1))/meshA%dx - meshA%rhog(1)*g
    pgA = meshA%pg(1)
    pgB = meshB%pg(1)
    if(innerBracket > 0)then
      !meshBが上流
      if(present(index))then
        if(index == meshBIndex)then
          pgB = getPg(meshB,mediumProps,fluidProps,ctrl,kind)
          rhog = getRhoG(meshB,mediumProps,fluidProps,ctrl,kind)
          krg = getKrG(meshB,mediumProps,fluidProps,ctrl,kind)
        else if(index == meshBIndex -1)then
          pgA = getPg(meshA,mediumProps,fluidProps,ctrl,kind)
          rhog = getRhoG(meshB,mediumProps,fluidProps,ctrl)
          krg = getKrG(meshB,mediumProps,fluidProps,ctrl)
        else
          rhog = getRhoG(meshB,mediumProps,fluidProps,ctrl)
          krg = getKrG(meshB,mediumProps,fluidProps,ctrl)
        endif
      else
        rhog = getRhoG(meshB,mediumProps,fluidProps,ctrl)
        krg = getKrG(meshB,mediumProps,fluidProps,ctrl)
      endif
    else
      !meshAが上流
      if(present(index))then
        if(index == meshBIndex-1)then
          pgA = getPg(meshA,mediumProps,fluidProps,ctrl,kind)
          rhog = getRhoG(meshA,mediumProps,fluidProps,ctrl,kind)
          krg = getKrG(meshA,mediumProps,fluidProps,ctrl,kind)
        else if(index == meshBIndex)then
          pgB = getPg(meshB,mediumProps,fluidProps,ctrl,kind)
          rhog = getRhoG(meshA,mediumProps,fluidProps,ctrl)
          krg = getKrG(meshA,mediumProps,fluidProps,ctrl)
        else
          rhog = getRhoG(meshA,mediumProps,fluidProps,ctrl)
          krg = getKrG(meshA,mediumProps,fluidProps,ctrl)
        endif
      else
        rhog = getRhoG(meshA,mediumProps,fluidProps,ctrl)
        krg = getKrG(meshA,mediumProps,fluidProps,ctrl)
      endif
    endif
    getMassFluxG = - mediumProps%permeability*rhog*krg/fluidProps%gasViscosity*((pgB-pgA)/meshA%dx - rhog*g)
    return
  end function getMassFluxG

  complex(8) function getMassFluxC(meshA,meshB,mediumProps,fluidProps,ctrl,index,kind,meshBIndex)
    implicit none
    type(meshObj),intent(in) :: meshA,meshB
    type(inputMediumPropertyObj) :: mediumProps
    type(inputFluidPropertyObj) :: fluidProps
    type(inputCtrlObj) :: ctrl
    integer,optional,intent(in) :: index,kind,meshBIndex
    complex(8) :: qw,C
    complex(8),allocatable :: phiW(:),rhoC(:)
    allocate(phiW(0:1),rhoC(0:1))

    qw = getMassFluxW(meshA,meshB,mediumProps,fluidProps,ctrl)
    if(real(qw) > 0)then
      !meshBが上流
      if(present(index))then
        if(index == meshBIndex)then
          C = getC(meshB,mediumProps,fluidProps,ctrl,kind)
        else
          C = getC(meshB,mediumProps,fluidProps,ctrl)
        endif
      else
        C = getC(meshB,mediumProps,fluidProps,ctrl)
      endif
    else
      !meshAが上流
      if(present(index))then
        if(index == meshBIndex-1)then
          C = getC(meshA,mediumProps,fluidProps,ctrl,kind)
        else
          C = getC(meshA,mediumProps,fluidProps,ctrl)
        endif
      else
        C = getC(meshA,mediumProps,fluidProps,ctrl)
      endif
    endif

    if(present(index))then
      if(index == meshBIndex)then
        phiW(0) = getPhiW(meshA,mediumProps,fluidProps,ctrl)
        phiW(1) = getPhiW(meshB,mediumProps,fluidProps,ctrl,kind)
        rhoC(0) = getRhoW(meshA,mediumProps,fluidProps,ctrl)*getC(meshA,mediumProps,fluidProps,ctrl)
        rhoC(1) = getRhoW(meshB,mediumProps,fluidProps,ctrl,kind)*getC(meshB,mediumProps,fluidProps,ctrl,kind)
      else if(index == meshBIndex-1)then
        phiW(0) = getPhiW(meshA,mediumProps,fluidProps,ctrl,kind)
        phiW(1) = getPhiW(meshB,mediumProps,fluidProps,ctrl)
        rhoC(0) = getRhoW(meshA,mediumProps,fluidProps,ctrl,kind)*getC(meshA,mediumProps,fluidProps,ctrl,kind)
        rhoC(1) = getRhoW(meshB,mediumProps,fluidProps,ctrl)*getC(meshB,mediumProps,fluidProps,ctrl)
      else
        phiW(0) = getPhiW(meshA,mediumProps,fluidProps,ctrl)
        phiW(1) = getPhiW(meshB,mediumProps,fluidProps,ctrl)
        rhoC(0) = getRhoW(meshA,mediumProps,fluidProps,ctrl)*getC(meshA,mediumProps,fluidProps,ctrl)
        rhoC(1) = getRhoW(meshB,mediumProps,fluidProps,ctrl)*getC(meshB,mediumProps,fluidProps,ctrl)
      endif
    else
      phiW(0) = getPhiW(meshA,mediumProps,fluidProps,ctrl)
      phiW(1) = getPhiW(meshB,mediumProps,fluidProps,ctrl)
      rhoC(0) = getRhoW(meshA,mediumProps,fluidProps,ctrl)*getC(meshA,mediumProps,fluidProps,ctrl)
      rhoC(1) = getRhoW(meshB,mediumProps,fluidProps,ctrl)*getC(meshB,mediumProps,fluidProps,ctrl)
    endif
    getMassFluxC = - (1-mediumProps%reflectionCoef)*phiw(1-meshBIndex)*mediumProps%effectiveDiffusionCoef &
    & *(rhoC(1)-rhoC(0))/meshA%dx + C*(1-mediumProps%reflectionCoef)*qw
  end function getMassFluxC

  subroutine fillMatrix(array,mesh,mediumProps,fluidProps,ctrl,index)
    implicit none
    integer :: i,j,n
    integer,optional :: index
    real(8) :: eps,array(:,:,:)
    type(meshObj),allocatable :: mesh(:)
    type(inputMediumPropertyObj) :: mediumProps
    type(inputFluidPropertyObj) :: fluidProps
    type(inputCtrlObj) :: ctrl

    n = size(mesh)
    eps = ctrl%epsilon
    if(present(index))then
      do i = 2,n-1
        do j = 1,3
          array(1,j,i) = imag(waterResidual(mesh(i-1),mesh(i),mesh(i+1),mediumProps,fluidProps,ctrl,index,j-1))/eps
          array(2,j,i) = imag(gasResidual(mesh(i-1),mesh(i),mesh(i+1),mediumProps,fluidProps,ctrl,index,j-1))/eps
          array(3,j,i) = imag(saltResidual(mesh(i-1),mesh(i),mesh(i+1),mediumProps,fluidProps,ctrl,index,j-1))/eps
        enddo
      enddo
    else
      do i = 2,n-1
        array(1,1,i) = -waterResidual(mesh(i-1),mesh(i),mesh(i+1),mediumProps,fluidProps,ctrl)
        array(2,1,i) = -gasResidual(mesh(i-1),mesh(i),mesh(i+1),mediumProps,fluidProps,ctrl)
        array(3,1,i) = -saltResidual(mesh(i-1),mesh(i),mesh(i+1),mediumProps,fluidProps,ctrl)
      enddo
    endif
  end subroutine fillMatrix

  subroutine fillMatrixBoundary(a,b,c,r,mesh,mediumProps,fluidProps,ctrl,injection)
    !境界条件の設定
    implicit none
    integer :: i,j,n
    real(8) eps,a(:,:,:),b(:,:,:),c(:,:,:),r(:,:,:)
    type(meshObj),allocatable :: mesh(:)
    type(inputMediumPropertyObj) :: mediumProps
    type(inputFluidPropertyObj) :: fluidProps
    type(inputCtrlObj) :: ctrl
    double precision :: injection

    eps = ctrl%epsilon
    n = size(mesh)
    do i = 1,3
      a(1,i,1) = imag(boundaryWaterResidual(mesh(1),mesh(2),mediumProps,fluidProps,ctrl,0,i-1))/eps
      b(1,i,1) = imag(boundaryWaterResidual(mesh(1),mesh(2),mediumProps,fluidProps,ctrl,1,i-1))/eps
      a(2,i,1) = imag(boundaryGasResidual(injection,mesh(1),mesh(2),mediumProps,fluidProps,ctrl,0,i-1))/eps
      b(2,i,1) = imag(boundaryGasResidual(injection,mesh(1),mesh(2),mediumProps,fluidProps,ctrl,1,i-1))/eps
      a(3,i,1) = imag(boundarySaltResidual(mesh(1),mesh(2),mediumProps,fluidProps,ctrl,0,i-1))/eps
      b(3,i,1) = imag(boundarySaltResidual(mesh(1),mesh(2),mediumProps,fluidProps,ctrl,1,i-1))/eps
    enddo
    r(1,1,1) = - boundaryWaterResidual(mesh(1),mesh(2),mediumProps,fluidProps,ctrl)
    r(2,1,1) = - boundaryGasResidual(injection,mesh(1),mesh(2),mediumProps,fluidProps,ctrl)
    r(3,1,1) = - boundarySaltResidual(mesh(1),mesh(2),mediumProps,fluidProps,ctrl)
    do i=1,3
      do j=1,3
        a(i,j,n) = 0
        c(i,j,n) = 0
      enddo
      a(i,i,n) = 1
    enddo
    r(1,1,n) = ctrl%initPw - mesh(n)%pw(1)
    r(2,1,n) = 0; r(3,1,n) = 0
  end subroutine
end module m_residual

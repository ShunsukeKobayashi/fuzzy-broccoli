!this is the sample of ctrl
module m_ctrl
	implicit none; private

	type, public :: ctrl_obj
		integer :: newtonmax,neq,ndof
		real(8) :: residual_ratio,correction_ratio,correction_limit,correction_limit_min,correction_limit_ratio
		real(8), allocatable :: max_residual(:),max_correction(:),zero_residual(:),zero_correction(:)
		integer, allocatable :: imax_res_loc(:,:),imax_cor_loc(:,:)
		real(8), allocatable :: dof_change(:,:,:,:),mean_flow(:,:,:,:)
		logical, allocatable :: residual_flag(:,:,:,:),correction_flag(:,:,:,:)
		logical :: total_res_flag,total_cor_flag
		real(8), allocatable :: delta(:),delta_PEtrans(:)
	contains
		procedure :: InitialSet
		procedure :: allocateLoc
		procedure :: ResidualCheck
		procedure :: CorrectionCheck
		procedure :: dealloc
	end type

	contains
!---------------------------------------------------------------------------------------------------------------------------
	subroutine CorrectionCheck(ctrl,r)
	implicit none
	class(ctrl_obj), intent(inout) :: ctrl
	real(8), intent(in) :: r(:,:,:,:)
	integer :: idof,ix,iy,iz,ndof,nx,ny,nz

	ndof=size(r,1)
	nx=size(r,2)
	ny=size(r,3)
	nz=size(r,4)

	ctrl%correction_flag=.false.

	ctrl%imax_cor_loc=0
	ctrl%max_correction=0

	do iz=1,nz
		do iy=1,ny
			do ix=1,nx
				do idof=1,ndof
					if(abs(r(idof,ix,iy,iz))>abs(ctrl%max_correction(idof))) then
						ctrl%imax_cor_loc(idof,1)=ix
						ctrl%imax_cor_loc(idof,2)=iy
						ctrl%imax_cor_loc(idof,3)=iz
						ctrl%max_correction(idof)=r(idof,ix,iy,iz)
					endif
					if(abs(r(idof,ix,iy,iz))<ctrl%correction_ratio*abs(ctrl%dof_change(idof,ix,iy,iz)) &
						&.or. abs(r(idof,ix,iy,iz))<ctrl%zero_correction(idof)) then
						ctrl%correction_flag(idof,ix,iy,iz)=.true.
					endif
				enddo
			enddo
		enddo
	enddo

	return
	end subroutine
!---------------------------------------------------------------------------------------------------------------------------
	subroutine ResidualCheck(ctrl,r)
	implicit none
	class(ctrl_obj), intent(inout) :: ctrl
	real(8), intent(in) :: r(:,:,:,:)
	integer :: ieq,ix,iy,iz,neq,nx,ny,nz

	neq=size(r,1)
	nx=size(r,2)
	ny=size(r,3)
	nz=size(r,4)

	ctrl%residual_flag=.false.

	ctrl%imax_res_loc=0
	ctrl%max_residual=0

	do iz=1,nz
		do iy=1,ny
			do ix=1,nx
				do ieq=1,neq
					if(abs(r(ieq,ix,iy,iz))>abs(ctrl%max_residual(ieq))) then
						ctrl%imax_res_loc(ieq,1)=ix
						ctrl%imax_res_loc(ieq,2)=iy
						ctrl%imax_res_loc(ieq,3)=iz
						ctrl%max_residual=r(ieq,ix,iy,iz)
					endif
					if(abs(r(ieq,ix,iy,iz))<ctrl%residual_ratio*abs(ctrl%mean_flow(ieq,ix,iy,iz)) .or. abs(r(ieq,ix,iy,iz))<ctrl%zero_residual(ieq)) then
						ctrl%residual_flag(ieq,ix,iy,iz)=.true.
					endif
				enddo
			enddo
		enddo
	enddo

	return
	end subroutine
!---------------------------------------------------------------------------------------------------------------------------
	subroutine InitialSet(ctrl,number)
	use stdio
	implicit none
	class(ctrl_obj), intent(inout) :: ctrl
	integer, intent(in), optional :: number

	type(file_obj) :: InputFile
	character(ncard_width) :: reading_card

	InputFile%filename="./input/ctrl.env"
	if(present(number)) then
		InputFile%nunit=number*10
		open(unit=InputFile%nunit,file=InputFile%filename,status="old",action="read")
	else
		open(newunit=InputFile%nunit,file=InputFile%filename,status="old",action="read")
	endif

	do
		call InputFile%readline(reading_card)
		if(InputFile%nerr/=0) exit
		if(trim(reading_card)=="#NEWTONMAX") then
			call InputFile%readline(reading_card)
			read(reading_card,*) ctrl%newtonmax
		endif
		if(trim(reading_card)=="#NEQUATION") then
			call InputFile%readline(reading_card)
			read(reading_card,*) ctrl%neq
		endif
	enddo
	ctrl%ndof=ctrl%neq
	close(InputFile%nunit)
	open(unit=InputFile%nunit,file=InputFile%filename,status="old",action="read")

	if(allocated(ctrl%delta)) deallocate(ctrl%delta)
	if(allocated(ctrl%delta_PEtrans)) deallocate(ctrl%delta_PEtrans)
	if(allocated(ctrl%zero_residual)) deallocate(ctrl%zero_residual)
	if(allocated(ctrl%zero_correction)) deallocate(ctrl%zero_correction)
	allocate(ctrl%delta(ctrl%ndof),ctrl%delta_PEtrans(ctrl%ndof))
	allocate(ctrl%zero_residual(ctrl%ndof),ctrl%zero_correction(ctrl%ndof))

	do
		call InputFile%readline(reading_card)
		if(InputFile%nerr/=0) exit
		if(trim(reading_card)=="#DELTA") then
			call InputFile%readline(reading_card)
			read(reading_card,*) ctrl%delta
		endif
		if(trim(reading_card)=="#TRANS_DAMPING") then
			call InputFile%readline(reading_card)
			read(reading_card,*) ctrl%delta_PEtrans
		endif
		if(trim(reading_card)=="#ZERO_CORRECTION") then
			call InputFile%readline(reading_card)
			read(reading_card,*) ctrl%zero_correction
		endif
		if(trim(reading_card)=="#CORRECTION_RATIO") then
			call InputFile%readline(reading_card)
			read(reading_card,*) ctrl%correction_ratio
		endif
		if(trim(reading_card)=="#CORRECTION_LIMIT") then
			call InputFile%readline(reading_card)
			read(reading_card,*) ctrl%correction_limit,ctrl%correction_limit_min,ctrl%correction_limit_ratio
		endif
		if(trim(reading_card)=="#ZERO_RESIDUAL") then
			call InputFile%readline(reading_card)
			read(reading_card,*) ctrl%zero_residual
		endif
		if(trim(reading_card)=="#RESIDUAL_RATIO") then
			call InputFile%readline(reading_card)
			read(reading_card,*) ctrl%residual_ratio
		endif
	enddo

	close(InputFile%nunit)



	return
	end subroutine
!---------------------------------------------------------------------------------------------------------------------------
	subroutine allocateLoc(ctrl,nx,ny,nz)
	implicit none
	class(ctrl_obj), intent(inout) :: ctrl
	integer, intent(in) :: nx,ny,nz

	if(allocated(ctrl%residual_flag)) deallocate(ctrl%residual_flag)
	if(allocated(ctrl%correction_flag)) deallocate(ctrl%correction_flag)
	if(allocated(ctrl%dof_change)) deallocate(ctrl%dof_change)
	if(allocated(ctrl%mean_flow)) deallocate(ctrl%mean_flow)
	if(allocated(ctrl%imax_res_loc)) deallocate(ctrl%imax_res_loc)
	if(allocated(ctrl%imax_cor_loc)) deallocate(ctrl%imax_cor_loc)
	if(allocated(ctrl%max_residual)) deallocate(ctrl%max_residual)
	if(allocated(ctrl%max_correction)) deallocate(ctrl%max_correction)
	allocate(ctrl%residual_flag(ctrl%neq,nx,ny,nz),ctrl%correction_flag(ctrl%ndof,nx,ny,nz))
	allocate(ctrl%dof_change(ctrl%ndof,nx,ny,nz),ctrl%mean_flow(ctrl%neq,nx,ny,nz))
	allocate(ctrl%imax_res_loc(ctrl%neq,3),ctrl%imax_cor_loc(ctrl%ndof,3))
	allocate(ctrl%max_residual(ctrl%neq),ctrl%max_correction(ctrl%ndof))

	return
	end subroutine
!---------------------------------------------------------------------------------------------------------------------------
	subroutine dealloc(ctrl)
	implicit none
	class(ctrl_obj), intent(inout) :: ctrl

	if(allocated(ctrl%delta)) deallocate(ctrl%delta)
	if(allocated(ctrl%delta_PEtrans)) deallocate(ctrl%delta_PEtrans)
	if(allocated(ctrl%zero_residual)) deallocate(ctrl%zero_residual)
	if(allocated(ctrl%zero_correction)) deallocate(ctrl%zero_correction)
	if(allocated(ctrl%residual_flag)) deallocate(ctrl%residual_flag)
	if(allocated(ctrl%correction_flag)) deallocate(ctrl%correction_flag)
	if(allocated(ctrl%dof_change)) deallocate(ctrl%dof_change)
	if(allocated(ctrl%mean_flow)) deallocate(ctrl%mean_flow)
	if(allocated(ctrl%imax_res_loc)) deallocate(ctrl%imax_res_loc)
	if(allocated(ctrl%imax_cor_loc)) deallocate(ctrl%imax_cor_loc)
	if(allocated(ctrl%max_residual)) deallocate(ctrl%max_residual)
	if(allocated(ctrl%max_correction)) deallocate(ctrl%max_correction)

	return
	end subroutine
!---------------------------------------------------------------------------------------------------------------------------
end module

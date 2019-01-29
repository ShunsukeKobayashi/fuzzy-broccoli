module m_inputFluidProperty
	implicit none; private

	type, public :: inputFluidPropertyObj
    real(8) :: molecularWeight,temperature,gamma,molarConcentration,waterViscosity,gasViscosity,waterSaturateRate,gasSaturateRate, &
		& rhow0,p0,kw
	contains
		procedure :: InitialSet
	end type

	contains
	subroutine InitialSet(inputFluidProperty,number)
	use stdio
	implicit none
	class(inputFluidPropertyObj), intent(inout) :: inputFluidProperty
	integer, intent(in), optional :: number
  integer :: i

	type(fileObj) :: InputFile
	character(ncardWidth) :: readingCard

	inputFluidProperty%gamma = 0.8
	!per unit salt concentration Viscosity increase rate
	inputFluidProperty%molarConcentration = 58.44d-3
	inputFluidProperty%molecularWeight = 44.01d-3
	inputFluidProperty%rhow0 = 1d3
	inputFluidProperty%p0 = 101325
	inputFluidProperty%kw = 2.2d9

	InputFile%filename="./input/inputFluidProperty.env"
	if(present(number)) then
		InputFile%nunit=number*10
		open(unit=InputFile%nunit,file=InputFile%filename,status="old",action="read")
	else
		open(newunit=InputFile%nunit,file=InputFile%filename,status="old",action="read")
	endif

	do
		call InputFile%readline(readingCard)
		if(InputFile%nerr/=0) exit
		if(trim(readingCard)=="#temperature") then
			call InputFile%readline(readingCard)
			read(readingCard,*) inputFluidProperty%temperature
		endif
		if(trim(readingCard)=="#waterViscosity") then
			call InputFile%readline(readingCard)
			read(readingCard,*) inputFluidProperty%waterViscosity
		endif
		if(trim(readingCard)=="#gasViscosity") then
			call InputFile%readline(readingCard)
			read(readingCard,*) inputFluidProperty%gasViscosity
		endif
		if(trim(readingCard)=="#waterSaturateRate") then
			call InputFile%readline(readingCard)
			read(readingCard,*) inputFluidProperty%waterSaturateRate
		endif
		if(trim(readingCard)=="#gasSaturateRate") then
			call InputFile%readline(readingCard)
			read(readingCard,*) inputFluidProperty%gasSaturateRate
		endif
	enddo
	close(InputFile%nunit)

	return
	end subroutine
end module

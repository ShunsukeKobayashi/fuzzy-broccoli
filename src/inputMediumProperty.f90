module m_inputMediumProperty
	implicit none; private

	type, public :: inputMediumPropertyObj
    real(8) :: compressibility,permeability,arufa,n,L,reflectionCoef,effectiveDiffusionCoef,initPhiW,initPhiG
	contains
		procedure :: InitialSet
	end type

	contains
	subroutine InitialSet(inputMediumProperty,number)
	use stdio
	implicit none
	class(inputMediumPropertyObj), intent(inout) :: inputMediumProperty
	integer, intent(in), optional :: number
  integer :: i

	type(fileObj) :: InputFile
	character(ncardWidth) :: readingCard

	InputFile%filename="./input/inputMediumProperty.env"
	if(present(number)) then
		InputFile%nunit=number*10
		open(unit=InputFile%nunit,file=InputFile%filename,status="old",action="read")
	else
		open(newunit=InputFile%nunit,file=InputFile%filename,status="old",action="read")
	endif

	do
		call InputFile%readline(readingCard)
		if(InputFile%nerr/=0) exit
		if(trim(readingCard)=="#compressibility") then
			call InputFile%readline(readingCard)
			read(readingCard,*) inputMediumProperty%compressibility
		endif
		if(trim(readingCard)=="#permeability") then
			call InputFile%readline(readingCard)
			read(readingCard,*) inputMediumProperty%permeability
		endif
		if(trim(readingCard)=="#arufa") then
			call InputFile%readline(readingCard)
			read(readingCard,*) inputMediumProperty%arufa
		endif
		if(trim(readingCard)=="#n") then
			call InputFile%readline(readingCard)
			read(readingCard,*) inputMediumProperty%n
		endif
		if(trim(readingCard)=="#L") then
			call InputFile%readline(readingCard)
			read(readingCard,*) inputMediumProperty%L
		endif
		if(trim(readingCard)=="#reflectionCoef") then
			call InputFile%readline(readingCard)
			read(readingCard,*) inputMediumProperty%reflectionCoef
		endif
		if(trim(readingCard)=="#effectiveDiffusionCoef") then
			call InputFile%readline(readingCard)
			read(readingCard,*) inputMediumProperty%effectiveDiffusionCoef
		endif
		if(trim(readingCard)=="#INITPHIW") then
			call InputFile%readline(readingCard)
			read(readingCard,*) inputMediumProperty%initPhiW
		endif
		if(trim(readingCard)=="#INITPHIG") then
			call InputFile%readline(readingCard)
			read(readingCard,*) inputMediumProperty%initPhiG
		endif
	enddo
	close(InputFile%nunit)

	return
	end subroutine
end module

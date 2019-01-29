module m_inputCtrl
	implicit none; private

	type, public :: inputCtrlObj
		real(8) :: period,interval,epsilon,maxNewton,reductionThreshold,initPw,initPg,initC
  contains
		procedure :: InitialSet
	end type

	contains
	subroutine InitialSet(inputCtrl,number)
	use stdio
	implicit none
	class(inputCtrlObj), intent(inout) :: inputCtrl
	integer, intent(in), optional :: number

	type(fileObj) :: InputFile
	character(ncardWidth) :: readingCard

	inputCtrl%epsilon = 1d-10
	inputCtrl%reductionThreshold = 1d5
	inputCtrl%maxNewton = 2000
	InputFile%filename="./input/inputCtrl.env"
	if(present(number)) then
		InputFile%nunit=number*10
		open(unit=InputFile%nunit,file=InputFile%filename,status="old",action="read")
	else
		open(newunit=InputFile%nunit,file=InputFile%filename,status="old",action="read")
	endif

	do
		call InputFile%readline(readingCard)
		if(InputFile%nerr/=0) exit
		if(trim(readingCard)=="#PERIOD") then
			call InputFile%readline(readingCard)
			read(readingCard,*) inputCtrl%period
		endif
		if(trim(readingCard)=="#INTERVAL") then
			call InputFile%readline(readingCard)
			read(readingCard,*) inputCtrl%interval
		endif
	enddo
	close(InputFile%nunit)

	return
	end subroutine

end module

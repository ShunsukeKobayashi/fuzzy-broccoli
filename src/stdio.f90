module stdio
	implicit none; private
	integer, parameter, public :: ncardWidth=511
	type, public :: fileObj
		integer :: ncard=ncardWidth
		character(13) :: fcard
		character(1) :: commark="!",keymark="#"
		integer :: nunit,nerr,iline=0
		character(:), allocatable :: filename
	contains
		procedure :: Initialize		!! ������
		procedure :: CountLine		!! ([EndKey]) �w�肵���L�[���[�h�������̓t�@�C���I���L�^�܂ł̍s���𐔂���
		procedure :: ReadError		!! (readingCard) �ǂݍ��݃G���[���N�������t�@�C�����A�s���A���̂Ƃ��̓ǂݍ��݃f�[�^���\�����Ē��~
		procedure :: OpenError		!! () �t�@�C�����J���Ƃ��ɃG���[���N�������t�@�C�������\�����Ē��~
		procedure :: ReadLine		!! (readingCard) �R�����g�s�������͋��s�łȂ��s�𕶎����Ƃ��ēǂݍ���
		procedure, private :: commentout	!! (readingCard) �R�����g���������菜��(�B��)
	end type

	contains
!-----------------------------------------------------------------------------------------------------------------------------------
	subroutine Initialize(file,ncard,commark,keymark)
	implicit none
	class(fileObj), intent(inout) :: file
	integer, intent(in), optional :: ncard
	character(1), intent(in), optional :: commark,keymark
	character(10) :: num_card

	if(present(ncard)) then
		file%ncard=ncard
		write(num_card,'(I0)') ncard
		file%fcard="(a"//trim(num_card)//")"
	endif

	if(present(commark)) file%commark=commark
	if(present(keymark)) file%keymark=keymark

	return
	end subroutine initialize
!-----------------------------------------------------------------------------------------------------------------------------------
	integer function CountLine(file,EndKey)
	implicit none
	class(fileObj), intent(inout) :: file
	character(len=*), intent(in), optional :: EndKey

	character(file%ncard) :: readingCard
	integer :: i,line0

	CountLine=0
	line0=file%iline
	if(present(EndKey)) then
		do
			call file%ReadLine(readingCard)
			if(trim(readingCard(1:len(EndKey)))==EndKey .or. file%nerr/=0) then
				exit
			else
				CountLine=CountLine+1
			endif
		enddo
	else
		do
			call file%ReadLine(readingCard)
			if(file%nerr/=0) then
				exit
			else
				CountLine=CountLine+1
			endif
		enddo
	endif

	!! �ǂݏo���J�n�ʒu�ɖ߂�
	do i=1,file%iline-line0
		backspace(file%nunit)
	enddo

	file%iline=line0

	return
	end function CountLine
!-----------------------------------------------------------------------------------------------------------------------------------
	subroutine ReadError(file,readingCard)
	implicit none
	class(fileObj), intent(in) :: file
	character(len=*) :: readingCard

	write(*,*) "ERROR!: Reading: ",trim(file%filename)
	write(*,*) "LINE:",file%iline
	write(*,*) "IMAGE: ", trim(readingCard)

	stop
	end subroutine ReadError
!-----------------------------------------------------------------------------------------------------------------------------------
	subroutine OpenError(file)
	implicit none
	class(fileObj), intent(in) :: file

	write(*,*) "Input file does not exist!"
	write(*,*) "File name: ",trim(file%filename)

	stop
	end subroutine OpenError
!-----------------------------------------------------------------------------------------------------------------------------------
	subroutine ReadLine(file,readingCard)
	implicit none
	class(fileObj), intent(inout) :: file
	character(len=*), intent(out) :: readingCard

	call file%Initialize(ncard=len(readingCard))
	do
		file%iline=file%iline+1
		read(file%nunit,'(a511)',iostat=file%nerr) readingCard
		if(file%nerr==0) then
			call file%commentout(readingCard)
			if(readingCard(1:1)/="") exit
		else
			exit
		endif
	enddo

	return
	end subroutine
!-----------------------------------------------------------------------------------------------------------------------------------
	pure subroutine commentout(file,readingCard)
	implicit none
	class(fileObj), intent(in) :: file
	character(len=*), intent(inout) :: readingCard
	integer :: icom

	readingCard=adjustl(readingCard)
	icom=index(readingCard,file%commark)
	if(icom>0) readingCard(icom:)=""

	return
	end subroutine commentout
!-----------------------------------------------------------------------------------------------------------------------------------
end module stdio

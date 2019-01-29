module m_inputMesh
	implicit none; private

	type, public :: inputMeshParam
		real(8) :: coordinate,meshLength,waterPressure,gasPressure,saltConcentration
	end type inputMeshParam

	type, public :: inputMeshObj
		integer :: numberOfMeshes
		type(inputMeshParam),allocatable :: mesh(:)
		contains
			procedure :: InitialSet
	end type inputMeshObj

	contains
		subroutine InitialSet(inputMesh,number)
		use stdio
		implicit none
		class(inputMeshObj), intent(inout) :: inputMesh
		integer, intent(in), optional :: number
		integer :: i,j=1,paramHead=0,paramLength=1,numberOfMeshParams=5
		character(20),allocatable :: meshParams(:)

		type(fileObj) :: InputFile
		character(ncardWidth) :: readingCard

		InputFile%filename="./input/inputMesh.env"
		if(present(number)) then
			InputFile%nunit=number*10
			open(unit=InputFile%nunit,file=InputFile%filename,status="old",action="read")
		else
			open(newunit=InputFile%nunit,file=InputFile%filename,status="old",action="read")
		endif

		do
			call InputFile%readline(readingCard)
			if(InputFile%nerr/=0) exit
				if(trim(readingCard)=="#NUMBEROFMESHES") then
				call InputFile%readline(readingCard)
				read(readingCard,*) inputMesh%numberOfMeshes
			endif
		enddo
		close(InputFile%nunit)
		open(unit=InputFile%nunit,file=InputFile%filename,status="old",action="read")

		allocate(inputMesh%mesh(inputMesh%numberOfMeshes))
		do
			call InputFile%readline(readingCard)
			if(InputFile%nerr/=0) exit
			if(trim(readingCard)=="#MESHESDATA") then
				do i=1,inputMesh%numberOfMeshes
					allocate(meshParams(numberOfMeshParams))
					paramHead=0;paramLength=1
					call InputFile%readline(readingCard)
					readingCard=adjustl(readingCard)
					do j=1,numberOfMeshParams
						paramLength=index(readingCard(paramHead:len(readingCard)),',')
						if(paramLength==0)then
							meshParams(j)=adjustl(readingCard(paramHead:len(readingCard)))
						else
							meshParams(j)=adjustl(readingCard(paramHead:(paramHead+paramLength-2)))
						endif
						paramHead=paramHead+paramLength
					enddo
					read(meshParams(1),*) inputMesh%mesh(i)%coordinate
					read(meshParams(2),*) inputMesh%mesh(i)%meshLength
					read(meshParams(3),*) inputMesh%mesh(i)%waterPressure
					read(meshParams(4),*) inputMesh%mesh(i)%gasPressure
					read(meshParams(5),*) inputMesh%mesh(i)%saltConcentration
					deallocate(meshParams)
				enddo
			endif
		enddo
		close(InputFile%nunit)
		return
	end subroutine
end module

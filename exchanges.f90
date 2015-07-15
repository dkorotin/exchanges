program exchanges

	use general
	use parameters
	use iomodule

  implicit none

  integer :: i,j

!  call read_hamilt()
  call read_crystal()

  write(stdout,*)
  write(stdout,'(a14,f12.9,a5)') 'Cell constant:', alat, 'Bohr'
  write(stdout,'(20a)') 'Cell vectors (rows):'
  do i=1,3
  	write(stdout,'(3f9.5)') cell(i,:)
  end do
  write(stdout,'(a26)') 'Atoms of the initial cell:'
  do i=1, natoms
  	write(stdout,'(a3,x,3f9.5)') atomlabel(i), tau(i,:)
  end do
  write(stdout,*)
  


  call clear()

end program exchanges

subroutine read_crystal()

	use iomodule
	use parameters, only : dp
	use general
	
	implicit none
	integer :: i,j
	character(40) :: line

	call open_input_file(iunsystem,'system.am')

	call find_section(iunsystem,'&cell')
	read(iunsystem,*) alat
	
	do i = 1,3
		read(iunsystem,*) cell(:,i)
	end do
	
	read(iunsystem,*) natoms
	allocate( tau(natoms,3) )
	allocate( atomlabel(natoms) )

	do i=1, natoms
		read(iunsystem,*) atomlabel(i), tau(i,:)
	end do

	close(iunsystem)

end subroutine read_crystal

subroutine clear()

	use general

	if(allocated(tau)) deallocate(tau)
	if(allocated(atomlabel)) deallocate(atomlabel)

end subroutine clear
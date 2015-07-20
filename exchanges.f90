program exchanges

	use general
	use parameters
	use iomodule

  implicit none

  integer :: i,j
  character(len=3) :: fmt='   ' ! is used for pretty output only

  real(dp) :: &
  						distance ! readius of the sphere for nearest neighbours search

  call read_hamilt()
  call read_crystal()
  
  ! Output of the readed values
  write(stdout,'(/,a14,f12.9,a5)') 'Cell constant:', alat, 'Bohr'
  write(stdout,'(20a)') 'Cell vectors (rows):'
  do i = 1, 3
  	write(stdout,'(3f9.5)') cell(i,:)
  end do
  write(stdout,'(a26)') 'Atoms of the initial cell:'
  do i = 1, natoms
  	write(stdout,'(a3,x,3f9.5)') atomlabel(i), tau(i,:)
  end do
  
  write(stdout,'(/,a71)') 'We have the following basis for the hamiltonian and the green function:'
  write(stdout,'(3x,i3,a20,i2,a8)') hdim, ' orbitals grouped in', nblocks, ' blocks:'
  do i=1, nblocks
  	write(fmt,'(i1)') block_dim(i)
  	write(stdout,'(5x,a3,a7,i2,a2,3x,i1,x,a,a11,'//adjustl(fmt)//'a11,a1)') &
  			adjustr(atomlabel(block_atom(i))), ' (atom ', block_atom(i), '):', &
  			block_dim(i),block_l(i), '-orbitals (', (orbitals(block_orbitals(i,j)), j=1, block_dim(i)), ')'
  end do
  
  call find_nnbrhds(natoms,tau,cell,block_atom(1),distance,taunew,parents)


  call clear()

end program exchanges

subroutine read_crystal()

	use iomodule
	use parameters, only : dp
	use general
	
	implicit none
	integer :: i,j
	character(len=3) :: dummy
	character :: l_symbol

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

	call find_section(iunsystem,'&basis')
	read(iunsystem,*) i, nblocks
	if(i .ne. hdim) stop 'Basis dimention in system.am and hamilt.am are unequal'

	allocate( block_atom(nblocks) )
	allocate( block_l(nblocks) )
	allocate( block_dim(nblocks) )
	allocate( block_orbitals(nblocks,7) )
	block_orbitals = 0

	do i = 1, nblocks
		read(iunsystem,*) dummy, block_atom(i), block_l(i), block_dim(i), block_orbitals(i,1:block_dim(i))
	end do

	! let's find the first orbital of each block - for simplicity
	allocate( block_start(nblocks) )
	block_start(1) = 1
	do i = 2, nblocks
		block_start(i) = block_start(i-1) + block_dim(i-1)
	end do

	call find_section(iunsystem,'&efermi')
	read(iunsystem,*) efermi

	close(iunsystem)

end subroutine read_crystal

subroutine read_hamilt()

	use iomodule
	use parameters, only : dp
	use general
	
	implicit none
	integer :: i,j,ispin, ik
	real(dp) :: Hre, Him

	call open_input_file(iunhamilt,'hamilt.am')

	call find_section(iunhamilt,'&hamilt')
	read(iunhamilt,*) nkp, hdim, nspin
	allocate( h(hdim,hdim,nkp,nspin) )
	allocate( wk(nkp,nspin) )
	allocate( xk(3,nkp,nspin) )

	h = cmplx(0.0,0.0)
	wk = 0.0
	xk = 0.0

	do ispin = 1, nspin
		do ik = 1, nkp
			read(iunhamilt,*) wk(ik,ispin), xk(:,ik,ispin)
			do i=1, hdim
				do j=1, hdim
					read(iunhamilt,*) Hre, Him
					h(i,j,ik,ispin) = cmplx(Hre,Him)
				end do
			end do
		end do
	end do
	
	close(iunhamilt)

end subroutine read_hamilt

subroutine clear()

	use general

	if( allocated(tau) ) deallocate(tau)
	if( allocated(atomlabel) ) deallocate(atomlabel)
	if( allocated(h) ) deallocate(h)
	if( allocated(wk) ) deallocate(wk)
	if( allocated(xk) ) deallocate(xk)
	if( allocated(block_atom) ) deallocate( block_atom )
	if( allocated(block_l) ) deallocate( block_l )
	if( allocated(block_dim) ) deallocate( block_dim )
	if( allocated(block_orbitals) ) deallocate( block_orbitals )
	if( allocated(block_start) ) deallocate(block_start)

end subroutine clear
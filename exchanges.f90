! Copyright (C) Dmitry Korotin dmitry@korotin.name

program exchange_parameters

	use general
	use parameters
	use iomodule

  implicit none

  integer :: i, j, time_start, time_end, count_rate, nz
  character(len=3) :: fmt='   ' ! is used for pretty output only
  real(dp) :: pos_delta
  complex(dp), allocatable :: z(:), G(:,:,:,:,:,:)
  complex(dp) :: zstep

  ! Input parameters are below:

  integer :: iverbosity ! verbosity of output
  
  ! Intergation mesh:
  integer :: nz1, nz2, nz3
  real(dp) :: emin, emax, height

  real(dp) :: distance ! distance for the nearest neighbours search
  										 ! could be used as integer do define a coordination sphere number
  										 ! if mode='csphere'
  character(len=10) :: mode ! mode for the nearest neighbours search

  namelist /exchanges/ nz1, nz2, nz3, height, emin, emax, distance, iverbosity, mode

  ! default values
	nz1 = 150
	nz2 = 350
	nz3 = 150
	height = 0.5
	emin = -30
	emax = 0.05
	distance = 8.d-1
	iverbosity = 0
	mode = 'distance'
	
  call system_clock(time_start,count_rate)

  read(stdin, exchanges, iostat=ios)
  if( ios .ne. 0 ) stop "Can't read input"
  
  call read_hamilt()
  call read_crystal()
  
  ! Output of the readed values
  write(stdout,'(/,a14,f12.9,a5)') 'Cell constant:', alat, 'Bohr'
  write(stdout,'(20a)') 'Cell vectors (rows):'
  do i = 1, 3
  	write(stdout,'(3f9.5)') cell(:,i)
  end do
  write(stdout,'(a26)') 'Atoms of the initial cell:'
  do i = 1, natoms
  	write(stdout,'(a3,x,3f9.5)') atomlabel(i), tau(:,i)
  end do
  
  write(stdout,'(/,a71)') 'We have the following basis for the hamiltonian and the green function:'
  write(stdout,'(3x,i3,a20,i2,a8)') hdim, ' orbitals grouped in', nblocks, ' blocks:'
  do i=1, nblocks
  	write(fmt,'(i1)') block_dim(i)
  	write(stdout,'(5x,a3,a7,i2,a2,3x,i1,x,a,a11,'//adjustl(fmt)//'a11,a1)') &
  			adjustr(atomlabel(block_atom(i))), ' (atom ', block_atom(i), '):', &
  			block_dim(i),block_l(i), '-orbitals (', (orbitals(block_orbitals(i,j)), j=1, block_dim(i)), ')'
  end do
  
  call atoms_list(mode,distance,1)

  ! Pretty output
  write(stdout,'(/a17,i4,a8)') 'We will consider ', nnnbrs,' atoms: '

  do i = 1, nnnbrs
  	pos_delta = sqrt( (taunew(1,i)-tau(1,block_atom(1)))**2 + &
  								(taunew(2,i)-tau(2,block_atom(1)))**2 + &
									(taunew(3,i)-tau(3,block_atom(1)))**2 )
  	write(stdout,'(a3,x,3f9.5,3x,a,f9.5,a)') atomlabel(parents(i)), taunew(:,i), '(', pos_delta, ')'
  end do


  ! define the z mesh
  ! i will generate nz+1 mesh points as a trick for zstep calculation below
  ! only nz points will be used in fact.

  nz = nz1 + nz2 + nz3
  allocate(z(nz+1))
  z = cmplx(0.0,0.0,dp)

  z(1) = cmplx(emin,0.d0,dp)

  zstep = cmplx(0.d0,(height/nz1),dp)
  do i = 2, nz1+1
    z(i) = z(i-1) + zstep
  end do

  zstep = cmplx((emax-emin)/nz2,0.d0,dp)
  do i = nz1+2, nz1+nz2
    z(i) = z(i-1) + zstep
  end do

  zstep = -1.d0*cmplx(0.d0,(height/nz3),dp)
  do i = nz1+nz2+1, nz+1
    z(i) = z(i-1) + zstep
  end do

  !Now compute full Green function for all atoms
  allocate(G(nz,nnnbrs,nnnbrs,MAXVAL(block_dim),MAXVAL(block_dim),nspin))

  call compute_g(nz,nnnbrs,nblocks,MAXVAL(block_dim),G,H,z(1:nz),parents(1:nnnbrs), &
                        taunew(:,1:nnnbrs),block_start(1:nblocks),block_dim(1:nblocks))

  if( allocated(z) ) deallocate(z)
  if( allocated(G) ) deallocate(G)
  call clear()

  call system_clock(time_end,count_rate)

  write(stdout,'(/5x,a15,i5,a8)') 'Execution time:', (time_end - time_start)/count_rate, 'seconds'

end program exchange_parameters

subroutine atoms_list(mode,distance,block_num)

	use iomodule
	use parameters, only : dp, maxnnbrs
	use general
	
	implicit none
	character(len=10), intent(in) :: mode
	real(dp), intent(in) :: distance
	integer, intent(in) :: block_num

	! temp storage
	real(dp) :: taunew_(3,maxnnbrs)
	integer :: parents_(maxnnbrs), nnnbrs_

	integer :: i,j, iblock
	logical :: have_atom_already

	call find_nnbrs(natoms,tau,cell,block_atom(block_num),distance,nnnbrs_,taunew_,parents_)

	! We want to consider all hamilt atoms anyway
	do i = 1, natoms
		call haa(taunew_, tau(:,i), have_atom_already)
		if( .not. have_atom_already ) then
			nnnbrs_ = nnnbrs_+1
			taunew_(:,nnnbrs_) = tau(:,i)
			parents_(nnnbrs_) = i
		end if
	end do

	nnnbrs = 0
	parents = -1
	taunew = -1000

	!filter atoms by l
	do i = 1, nnnbrs_
		do iblock = 1, nblocks
			if( block_atom(iblock) .eq. parents_(i) .and. block_l(iblock) .eq. block_l(block_num) ) then
				nnnbrs = nnnbrs + 1
				parents(nnnbrs) = parents_(i)
				taunew(:,nnnbrs) = taunew_(:,i)
			end if
		end do 
	end do

end subroutine atoms_list

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
	allocate( tau(3,natoms) )
	allocate( atomlabel(natoms) )

	do i=1, natoms
		read(iunsystem,*) atomlabel(i), tau(:,i)
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

subroutine inverse_complex_matrix(dim,a)
  
  use parameters, only : dp
  implicit none

  integer :: dim, info, ipiv(dim)
  complex(dp) :: a(dim,dim), work(dim)

  call ZGETRF(dim,dim,a,dim,ipiv,info)
  if(info /= 0) stop "inverse_complex_matrix Error in ZGETRF"
  call ZGETRI(dim,a,dim,ipiv,work,dim,info)
  if(info /= 0) stop"inverse_complex_matrix Error in ZGETRI"

end subroutine
! Copyright (C) Dmitry Korotin dmitry@korotin.name

program exchange_parameters

	use general
	use parameters
	use iomodule

  implicit none

  integer :: i, j, time_start, time_end, count_rate, nz, ia, ja, idim, jdim, iz, istart, jstart, iend, jend
  character(len=3) :: fmt='   ' ! is used for pretty output only
  real(dp) :: pos_delta
  complex(dp), allocatable :: z(:), G(:,:,:,:,:,:), delta(:,:), tmp1(:,:), hksum(:,:,:)
  complex(dp) :: zstep, tmp2
  real(dp), allocatable :: Jexc(:,:), Jexc0(:), Jorb(:,:,:,:)
  

  ! Input parameters are below:
  
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
	mode = 'distance'
	
  call system_clock(time_start,count_rate)

  write(stdout,'(/,5x,a66)') '------------------------------------------------------------------'
  write(stdout,'(5x,a66)')   '                      Program EXCHANGES                           '
  write(stdout,'(5x,a66)')   ' for calculation of exchange parameters of the Heisenberg model.  '
	write(stdout,'(/,5x,a66)') 'Please cite "D. M. Korotin et al., Phys. Rev. B 91, 224405 (2015)"'
	write(stdout,'(5x,a66)')   '    in publications or presentations arising from this work.      '
	write(stdout,'(5x,a66,/)') '------------------------------------------------------------------'

  read(stdin, exchanges, iostat=ios)
  if( ios .ne. 0 ) stop "Can't read input"
  
  call read_hamilt()
  call read_crystal()

  ! debug
 	if( iverbosity .ge. 3) then
		! output H(0)
		allocate( hksum(hdim,hdim,nspin) )
		hksum = cmplx(0.0,0.0,dp)

		do i=1, nkp
			hksum(:,:,:) = hksum(:,:,:) + h(:,:,i,:)*wk(i)
		end do

		do i=1, nspin
			write(stdout,'(/,5x,a13,i2,a1)') 'H(0) for spin', i, ':'
			call output_matrix_by_blocks(hdim,dreal(hksum(:,:,i)))
		end do

		deallocate( hksum )	
	end if
	! end of debug
  
  ! Output of the readed values
  write(stdout,'(/,5x,a14,f12.9,a5)') 'Cell constant:', alat, 'Bohr'
  write(stdout,'(5x,20a)') 'Cell vectors (rows):'
  do i = 1, 3
  	write(stdout,'(7x,3f9.5)') cell(:,i)
  end do
  write(stdout,'(5x,a26)') 'Atoms of the initial cell:'
  do i = 1, natoms
  	write(stdout,'(7x,a3,x,3f9.5)') atomlabel(i), tau(:,i)
  end do
  
  write(stdout,'(/,5x,a71)') 'We have the following basis for the hamiltonian and the green function:'
  write(stdout,'(5x,i3,a20,i2,a8)') hdim, ' orbitals grouped in', nblocks, ' blocks:'
  do i=1, nblocks
  	write(fmt,'(i1)') block_dim(i)
  	write(stdout,'(5x,a3,a7,i2,a2,3x,i1,x,a,a11,'//adjustl(fmt)//'a11,a1)') &
  			adjustr(atomlabel(block_atom(i))), ' (atom ', block_atom(i), '):', &
  			block_dim(i),block_l(i), '-orbitals (', (orbitals(block_orbitals(i,j)), j=1, block_dim(i)), ')'
  end do
  
  call atoms_list(mode,distance,1)

  ! Pretty output
  write(stdout,'(/5x,a17,i4,a8)') 'We will consider ', nnnbrs,' atoms: '

  do i = 1, nnnbrs
  	pos_delta = sqrt( (taunew(1,i)-tau(1,block_atom(1)))**2 + &
  								(taunew(2,i)-tau(2,block_atom(1)))**2 + &
									(taunew(3,i)-tau(3,block_atom(1)))**2 )
  	write(stdout,'(5x,i2,a2,a3,x,3f9.5,3x,a12,f9.5,a)') i, ': ', atomlabel(parent(i)), taunew(:,i), '( distance =', pos_delta, ')'
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

  call compute_g(nz,nnnbrs,nblocks,MAXVAL(block_dim),G,H,z(1:nz),parent(1:nnnbrs), &
                        taunew(:,1:nnnbrs),block_start(1:nblocks),block_dim(1:nblocks))

  allocate(delta(hdim,hdim))
  call compute_delta(H,delta)

  ! debug
  if( iverbosity .ge. 3) then
  	write(stdout,'(/,5x,a15)') 'Computed delta:'
  	call output_matrix_by_blocks(hdim,dreal(delta))
  end if
  ! end of debug

  allocate(Jexc(nnnbrs,nnnbrs))
  Jexc = cmplx(0.0,0.0,dp)

  allocate( Jorb(nnnbrs,nnnbrs,MAXVAL(block_dim),MAXVAL(block_dim)))
  Jorb = cmplx(0.0,0.0,dp)

  allocate( tmp1(MAXVAL(block_dim),MAXVAL(block_dim)) )

  DO ia=1,nnnbrs
    DO ja=ia+1,nnnbrs
        ! arrays indexes
        istart = block_start(parent(ia))
        idim = block_dim(parent(ia))        
        iend = block_start(parent(ia)) + block_dim(parent(ia)) - 1

        jstart = block_start(parent(ja))
        jdim = block_dim(parent(ja))
        jend = block_start(parent(ja)) + block_dim(parent(ja)) - 1
        
        IF(idim .NE. jdim ) THEN
          write(stdout,*) ia,ja,idim,jdim
          stop 'Not equal subblocks size'
        END IF        

        DO iz=1,nz
          zstep = z(iz+1) - z(iz)

          tmp1 = cmplx(0.0,0.0,dp)

          tmp1(1:idim,1:idim) = MATMUL( &
                            MATMUL(delta(istart:iend,istart:iend),G(iz,ia,ja,1:idim,1:jdim,2)), &
                            MATMUL(delta(jstart:jend,jstart:jend),G(iz,ja,ia,1:jdim,1:idim,1)) &
                            )

          Jorb(ia,ja,1:idim,1:idim) = &
              Jorb(ia,ja,1:idim,1:idim) + DIMAG(tmp1(1:idim,1:idim)*zstep)

        END DO

        ! Trace
        DO i=1,idim
          Jexc(ia,ja) = Jexc(ia,ja) + Jorb(ia,ja,i,i)
        END DO

    END DO
  END DO
  deallocate(tmp1)

  Jexc = (-1.d0/tpi)*Jexc
  Jorb = (-1.d0/tpi)*Jorb

  DO ia = 1, nnnbrs
    DO ja = ia+1, nnnbrs
        ! Compute distance between atoms for pretty output
        pos_delta= SQRT((taunew(1,ia) - taunew(1,ja))**2+ &
                        (taunew(2,ia) - taunew(2,ja))**2+ &
                        (taunew(3,ia) - taunew(3,ja))**2 )
        !
        write(stdout,'(/5x,a35,i3,a4,i3,a2,f7.3,a7,i5,a13,f6.3,a)') &
              'Exchange interaction between atoms ', ia, 'and', ja, ': ', &
              Jexc(ia,ja)*1.0d3, ' meV = ', INT(Jexc(ia,ja)/kb_ev), ' K (distance:',pos_delta,')'
        
        write(stdout,'(/7x,a62)') "Orbital exchange interaction matrix J_{i,j,m,n} (in K and meV):"
        DO i=1,block_dim(parent(ia))
          write(stdout, '(7x,5i5,8x,5f9.5)' ) INT(Jorb(ia,ja,i,:)/kb_ev), Jorb(ia,ja,i,:)*1.d3
        END DO

    END DO
  END DO

  write(stdout,*)
  write(stdout,*) '    Computed orbitals occupations should coincide with your DFT results'
  write(stdout,*) '    If they differ significantly - check your integration contour.'
  DO ia=1,nnnbrs    
    DO j=1,nspin
      write(stdout,'(/5x,a8,i3,a5,i2)') 'For atom', ia, 'spin', j
      DO i = 1, block_dim(parent(ia))
          tmp2 = cmplx(0.0,0.0,dp)
          DO iz=1,nz
            zstep = z(iz+1) - z(iz)
            tmp2 = tmp2 + ((-1.d0/pi)*DIMAG(G(iz,ia,ia,i,i,j)*zstep))
          END DO
        write(stdout,'(7x,a8,i2,a13,f6.3)') 'Orbital', i, ' occupation: ', DREAL(tmp2)
      END DO
    END DO
  END DO


  if( allocated(z) ) deallocate(z)
  if( allocated(G) ) deallocate(G)
  if( allocated(delta) ) deallocate(delta)
  if( allocated(Jorb) ) deallocate(Jorb)
  if( allocated(Jexc) ) deallocate(Jexc)
  call clear()

  call system_clock(time_end,count_rate)

  write(stdout,'(/5x,a15,i5,a8)') 'Execution time:', (time_end - time_start)/count_rate, 'seconds'

end program exchange_parameters

subroutine compute_delta(h,delta)
  use parameters, only : dp
  use general, only : hdim, nspin, nkp, wk
  
  implicit none
  complex(dp) :: h(hdim,hdim,nkp,nspin), delta(hdim,hdim)
  integer :: ik

  delta = cmplx(0.0,0.0,dp)
  
  ! delta = h(spin_up) - h(spin_down)
  do ik=1, nkp
    delta(:,:) = delta(:,:) + wk(ik)*( h(:,:,ik,1) - h(:,:,ik,2) )
  end do

  delta = dreal(delta) 

END SUBROUTINE compute_delta

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
	integer :: parent_(maxnnbrs), nnnbrs_

	integer :: i,j, iblock
	logical :: have_atom_already

	call find_nnbrs(natoms,tau,cell,block_atom(block_num),distance,nnnbrs_,taunew_,parent_)

	! We want to consider all hamilt atoms anyway
	do i = 1, natoms
		call haa(taunew_, tau(:,i), have_atom_already)
		if( .not. have_atom_already ) then
			nnnbrs_ = nnnbrs_+1
			taunew_(:,nnnbrs_) = tau(:,i)
			parent_(nnnbrs_) = i
		end if
	end do

	nnnbrs = 0
	parent = -1
	taunew = -1000

	!filter atoms by l
	do i = 1, nnnbrs_
		do iblock = 1, nblocks
			if( block_atom(iblock) .eq. parent_(i) .and. block_l(iblock) .eq. block_l(block_num) ) then
				nnnbrs = nnnbrs + 1
				parent(nnnbrs) = parent_(i)
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
	allocate( block_start(nblocks) )
	block_orbitals = 0

	do i = 1, nblocks
		read(iunsystem,*) dummy, block_atom(i), block_l(i), block_dim(i), block_start(i), block_orbitals(i,1:block_dim(i))
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

	call find_section(iunhamilt,'&nspin')
	read(iunhamilt,*) nspin

	call find_section(iunhamilt,'&nkp')
	read(iunhamilt,*) nkp

	call find_section(iunhamilt,'&dim')
	read(iunhamilt,*) hdim

	allocate( h(hdim,hdim,nkp,nspin) )
	allocate( wk(nkp) )
	allocate( xk(3,nkp) )

	h = cmplx(0.0,0.0)
	wk = 0.0
	xk = 0.0

	call find_section(iunhamilt,'&kpoints')
	do ik=1, nkp
		read(iunhamilt,*) wk(ik), xk(:,ik)
	end do
			
	call find_section(iunhamilt,'&hamiltonian')
	
	do ispin = 1, nspin
		do ik = 1, nkp
			do i=1, hdim
				do j=i, hdim
					read(iunhamilt,*) Hre, Him
					h(i,j,ik,ispin) = cmplx(Hre,Him)
					h(j,i,ik,ispin) = dconjg( h(i,j,ik,ispin) )
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
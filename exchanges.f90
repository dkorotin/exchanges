! Copyright (C) Dmitry Korotin dmitry@korotin.name

program exchange_parameters

  use general
  use parameters
  use iomodule
  use omp_lib

  implicit none

  integer :: i, j, time_start, time_end, count_rate, nz, ia, ja, idim, jdim, iz, istart, jstart, iend, jend
  character(len=3) :: fmt='   ' ! is used for pretty output only
  real(dp) :: pos_delta(3)
  complex(dp), allocatable :: z(:), G(:,:,:,:,:,:), delta(:,:), tmp1(:,:), tmp2(:,:), hksum(:,:,:)
  complex(dp) :: zstep, tmp3
  real(dp), allocatable :: Jexc(:,:), Jexc0(:), Jorb(:,:,:,:)
  

  ! Input parameters are below:
  
  ! Intergation mesh:
  integer :: nz1, nz2, nz3
  real(dp) :: emin, emax, height
  
  integer :: atom_of_interest
  character(len=1) :: l
  
  real(dp) :: distance ! distance for the nearest neighbours search
                       ! could be used as integer do define a coordination sphere number
                       ! if mode='csphere'
  character(len=10) :: mode ! mode for the nearest neighbours search

  namelist /exchanges/ nz1, nz2, nz3, height, emin, emax, distance, iverbosity, mode, l, atom_of_interest

  ! default values
  nz1 = 150
  nz2 = 3000
  nz3 = 150
  height = 0.01
  emin = -100.0
  emax = 0.0
  distance = 8.d-1
  mode = 'distance'
  l = 'd'
  atom_of_interest = -100
  
  call system_clock(time_start,count_rate)

  write(stdout,'(/,5x,a66)') '------------------------------------------------------------------'
  write(stdout,'(5x,a66)')   '                      Program EXCHANGES                           '
  write(stdout,'(5x,a66)')   ' for calculation of exchange parameters of the Heisenberg model.  '
  write(stdout,'(/,5x,a66)') 'Please cite "D. M. Korotin et al., Phys. Rev. B 91, 224405 (2015)"'
  write(stdout,'(5x,a66)')   '    in publications or presentations arising from this work.      '
  write(stdout,'(5x,a66,/)') '------------------------------------------------------------------'
  write(stdout,'(5x,a66,/)') '   We are using the model with the exchange term defined as:      '
  write(stdout,'(5x,a66,/)') '                   H = \sum_ij J_{ij} e_i e_j,                    '
  write(stdout,'(5x,a66,/)') ' where e_i,j are unit vectors and sum runs once over ions pairs   '
  write(stdout,'(5x,a66,/)') '------------------------------------------------------------------'

  write(stdout,'(/,5x,a11,i3,a8,/)') 'Running in ', OMP_get_max_threads(), ' threads'

  read(stdin, exchanges, iostat=ios)
  if( ios .ne. 0 ) stop "Can't read input"

  call read_hamilt()
  call read_crystal()

  if(atom_of_interest == -100) atom_of_interest = block_atom(1)

  ! let's find emin as the lowest hamiltonian eigenvalue if emin is not set
  if( emin == -100.0 ) call find_emin(emin)
  
  ! Output of the readed values

  write(stdout,'(5x,a34)') 'Parameters of integration contour:'
  write(stdout,'(5x,a7,f7.3,a9,f6.3,a11,f5.3,a5)') 'emin = ', emin, '  emax = ', emax, '  height = ', height,' (eV)'
  write(stdout,'(5x,a6,i4,a8,i4,a8,i4)') 'nz1 = ', nz1, '  nz2 = ', nz2, '  nz3 = ', nz3

  ! check integration countur
  if ( (emax - emin) / nz2 > height ) then
    write(stdout, '(/,5x,a55)' ) '*******************************************************'
    write(stdout, '(5x,a55)' )   'WARNING: increase nz2 or check the integration contour!'
    write(stdout, '(5x,a55)' )   '*******************************************************'
  end if
  
  write(stdout,'(/,5x,a20,f13.9,a5)') 'Cell constant (alat):', alat, 'Bohr'
  write(stdout,'(5x,20a)') 'Cell vectors (rows):'
  do i = 1, 3
    write(stdout,'(7x,3f9.5)') cell(:,i)
  end do
  write(stdout,'(5x,a26)') 'Atoms of the initial cell:'
  do i = 1, natoms
    write(stdout,'(7x,a3,x,3f9.5)') atomlabel(i), tau(:,i)
  end do
  
  write(stdout,'(/,5x,a71)') 'We have the following basis for the hamiltonian and the green function:'
  write(stdout,'(5x,i3,a20,i3,a8)') hdim, ' orbitals grouped in', nblocks, ' blocks:'
  do i=1, nblocks
    write(fmt,'(i1)') block_dim(i)
    write(stdout,'(5x,a3,a7,i2,a2,3x,i1,x,a,a11,'//adjustl(fmt)//'a11,a1)') &
        adjustr(atomlabel(block_atom(i))), ' (atom ', block_atom(i), '):', &
        block_dim(i),block_l(i), '-orbitals (', (orbitals(block_orbitals(i,j)), j=1, block_dim(i)), ')'
  end do
  
  call atoms_list(mode,distance,atom_of_interest,l)

  ! Pretty output
  write(stdout,'(/5x,a17,i4,a8)') 'We will consider ', nnnbrs,' atoms: '

  do i = 1, nnnbrs
    pos_delta =   taunew(:,i)-tau(:,atom_of_interest)
    write(stdout,'(5x,i2,a2,a3,x,3f9.5,3x,a12,f9.5,a7)') i, ': ', &
      atomlabel( block_atom( parent(i) ) ), taunew(:,i), '( distance =', norm2(pos_delta), ' alat )'
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

  allocate( tmp1(MAXVAL(block_dim),MAXVAL(block_dim)), tmp2(MAXVAL(block_dim),MAXVAL(block_dim)) )

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

        DO iz=1,nz-1
          zstep = z(iz+1) - z(iz)

          tmp1 = cmplx(0.0,0.0,dp)

          if (iz > 1) then
            tmp1 = tmp2
          else
            tmp1(1:idim,1:idim) = MATMUL( &
                              MATMUL(delta(istart:iend,istart:iend),G(iz,ia,ja,1:idim,1:jdim,2)), &
                              MATMUL(delta(jstart:jend,jstart:jend),G(iz,ja,ia,1:jdim,1:idim,1)) &
                              )
          end if

          tmp2(1:idim,1:idim) = MATMUL( &
                            MATMUL(delta(istart:iend,istart:iend),G(iz+1,ia,ja,1:idim,1:jdim,2)), &
                            MATMUL(delta(jstart:jend,jstart:jend),G(iz+1,ja,ia,1:jdim,1:idim,1)) &
                            )

          ! integration using trapezoidal rule
          Jorb(ia,ja,1:idim,1:idim) = &
              Jorb(ia,ja,1:idim,1:idim) + DIMAG( ( tmp1(1:idim,1:idim)+tmp2(1:idim,1:idim) )/2.0*zstep)

        END DO

        ! Trace
        DO i=1,idim
          Jexc(ia,ja) = Jexc(ia,ja) + Jorb(ia,ja,i,i)
        END DO

    END DO
  END DO
  deallocate(tmp1, tmp2)

  Jexc = (-1.d0/tpi)*Jexc
  Jorb = (-1.d0/tpi)*Jorb

  DO ia = 1, nnnbrs
    DO ja = ia+1, nnnbrs
        ! Compute distance between atoms for pretty output
        pos_delta= taunew(:,ia) - taunew(:,ja)
        !
        write(stdout,'(/5x,a35,a3,a,i3,a6,a3,a,i3,a4,f7.3,a7,i5,a13,f6.3,a)') &
              'Exchange interaction between atoms ', &
              atomlabel( block_atom( parent(ia) ) ), '(', ia, ') and ', &
              atomlabel( block_atom( parent(ja) ) ), '(', ja, ') : ', &
              Jexc(ia,ja)*1.0d3, ' meV = ', INT(Jexc(ia,ja)/kb_ev), ' K (distance:',norm2(pos_delta),')'

        write(stdout,'(5x,a35,3f8.4)') 'Connecting vector in cell vectors: ', pos_delta
        write(stdout,'(5x,a35,3f8.4)') 'Connecting vector in Cart. coords: ', MATMUL( pos_delta, TRANSPOSE(cell) )
        
        write(stdout,'(/7x,a62)') "Orbital exchange interaction matrix J_{i,j,m,n} (in K and meV):"
        DO i=1,block_dim(parent(ia))
          write(stdout, '(7x,5i5,8x,5f9.5)' ) INT(Jorb(ia,ja,i,1:block_dim(parent(ia)))/kb_ev), &
                                              Jorb(ia,ja,i,1:block_dim(parent(ia)))*1.d3
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
          tmp3 = cmplx(0.0,0.0,dp)
          DO iz=1,nz-1
            zstep = z(iz+1) - z(iz)
            tmp3 = tmp3 + ((-1.d0/pi)*( DIMAG(G(iz,ia,ia,i,i,j)) + DIMAG(G(iz+1,ia,ia,i,i,j)) )/2.0*zstep)
          END DO
        write(stdout,'(7x,a8,i2,a13,f6.3)') 'Orbital', i, ' occupation: ', DREAL(tmp3)
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

subroutine find_emin(emin)
  use parameters, only : dp
  use general, only : h, nkp, wk, nspin, hdim, efermi
  use iomodule, only: stdout, iverbosity, output_matrix_by_blocks

  implicit none

  complex(dp), allocatable :: hksum(:,:,:), work(:)
  real(dp), allocatable :: ev(:), rwork(:)
  integer :: i, j, info, lwork
  real(dp), intent(inout) :: emin

  allocate( hksum(hdim,hdim,nspin) )
  hksum = cmplx(0.0,0.0,dp)

  do i=1, nkp
    hksum(:,:,:) = hksum(:,:,:) + h(:,:,i,:)*wk(i)
  end do

  allocate( ev(hdim) )
  emin = 0.0

  allocate(rwork(3*hdim-2))

  ! obtain optimal lwork parameter
  allocate( work(1) )
  
  call zheev('N', 'U', hdim, hksum(:,:,1), hdim, ev, work, -1, rwork, info)
  lwork = work(1)
  deallocate(work)

  allocate( work(lwork) )
  do i = 1, nspin
    call zheev('N', 'U', hdim, hksum(:,:,i), hdim, ev, work, lwork, rwork, info)
    emin = min( ev(1),emin )
  end do
  deallocate( work, rwork, ev )

  emin = floor( emin - efermi -2.0)

  ! debug
  if( iverbosity .ge. 3) then
    do i=1, nspin
      write(stdout,'(/,5x,a13,i2,a1)') 'H(0) for spin', i, ':'
      call output_matrix_by_blocks(hdim,dreal(hksum(:,:,i)))
    end do
  ! end of debug
  end if

  deallocate( hksum ) 
end subroutine

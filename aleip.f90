PROGRAM ALEIP

USE kinds, ONLY: DP, nwan, nks, natoms, nspin, ZERO, ZEROC, pi, tpi, &
                 rytoev, K_BOLTZMANN_RY, iout, igeometry, nmaxatoms

IMPLICIT NONE

INTEGER                  :: i, j, iz, nz, nz1, nz2, nz3, nwand
!INTEGER                 ::  nblocks
INTEGER                  :: ia,ja,istart,idim,iend,jstart,jdim,jend
!INTEGER, allocatable     :: parent(:)       !tmp
!INTEGER, allocatable     :: block_start(:)  !tmp
!INTEGER, allocatable     :: block_dim(:)    !tmp
!REAL(DP), allocatable    :: dist(:,:)       !tmp
INTEGER                  :: parent(nmaxatoms)       !tmp
INTEGER                  :: block_start(nmaxatoms)  !tmp
INTEGER                  :: block_dim(nmaxatoms)    !tmp
REAL(DP)                 :: dist(3,nmaxatoms)       !tmp

REAL(DP)                 :: emin, emax, height, ef
COMPLEX(DP)              :: zstep, tmp2
COMPLEX(DP), allocatable :: H(:,:,:,:)       ! hamiltonian 
COMPLEX(DP), allocatable :: G(:,:,:,:,:,:) !  
COMPLEX(DP), allocatable :: z(:)           ! mesh for the integration, while calculating J
COMPLEX(DP), allocatable :: delta(:,:)     ! spin difference 
COMPLEX(DP), allocatable :: tmp1(:,:)
REAL(DP), allocatable    :: wk(:)          ! weights of k-points
REAL(DP), allocatable    :: xk(:,:)        ! mesh of k-points PWSCF: 2Pi/cell(1), cell(1) is in Au
REAL(DP), allocatable    :: Jexc(:,:), Jorb(:,:,:,:)

NAMELIST /geometry/  parent, block_start, block_dim, dist

! temprorary define some variables 
OPEN (10, file = 'hamilt', status = 'old', form = 'formatted')
read(10,*) nks, nwan
CLOSE (10)

nks = nks / 2  !FOR PWSCF, where hamiltonian is written for 2 spins
nspin  = 2 

call read_pwscf(nz1, nz2, nz3, height, emin, emax, ef, nwand, natoms)

if (natoms>nmaxatoms) stop 'increase nmaxatoms'

OPEN (iout, file = 'aleip.out', status = 'unknown', form = 'formatted')
write(iout,*) '---> INPUT PARAMETERS <---'
write(iout,*) ''
write(iout,'(a34,i3)') 'Number of Wannier functions:', nwan
write(iout,'(a34,i5)') 'Number of k-points:', nks
write(iout,'(a34,i4,i4,i4)') 'Number of points for integration:',nz1,nz2,nz3
write(iout,'(a34,i3)') 'Number of orbitals on imp. site:',  nwand 
write(iout,'(a34,i3)') 'Number of impurity sites:', natoms
write(iout,'(a34,3f8.5)') 'Emin, Emax, height:', emin, emax, height 
write(iout,'(a34,f8.5)') 'Fermi energy:', ef 

allocate(H(nwan, nwan, nks, nspin))
allocate(wk(nks))
allocate(xk(3,nks))
allocate(delta(nwan,nwan))

!allocate(block_start(nblocks))
!allocate(block_dim(nblocks))
!allocate(dist(3,natoms))    !tmp
!allocate(parent(natoms))    
!allocate(block_start(natoms))
!allocate(block_dim(natoms))

!default values
dist(1,1) = 0.d0
dist(2,1) = 0.d0
dist(3,1) = 0.d0
dist(1,2) = 0.d0
dist(2,2) = 0.d0
dist(3,2) = 0.d0
parent(1) = 1
parent(2) = 2
block_start(1) = 1
block_start(2) = 6
block_dim(1) = 5
block_dim(2) = 5

open (igeometry, file ='geometry.in')
read (igeometry, nml = geometry)
close (igeometry)
do i=1, natoms
  write(iout,'(a16,i2,a3,3f8.5)') 'Position of ion',i,':',(dist(j,i),j=1,3)  
enddo
do i=1, natoms
  write(iout,'(a29,i2,a3,i3)') 'Index of 1st orbital for ion',i,':',block_start(i)  
  write(iout,'(a29,i2,a3,i3)') 'Dimension of loc ham for ion',i,':',block_dim(i)  
enddo

delta = ZEROC

write(iout,*) ''
write(iout,*) '---> OUTPUT <---'
CALL read_hamilt(H,xk,wk)
CALL compute_delta(H,wk,delta)

!make a mesh for the integration
nz = nz1 + nz2 + nz3
allocate(G(nz,natoms,natoms,nwand,nwand,nspin))
allocate(z(nz+1))
z = ZEROC
z(1) = dcmplx(emin,0.d0)
zstep = dcmplx(0.d0,(height/nz1))
DO i = 2, nz1+1
   z(i) = z(i-1) + zstep
END DO
zstep = dcmplx((emax-emin)/nz2,0.d0)
DO i = nz1+2, nz1+nz2
   z(i) = z(i-1) + zstep
END DO
zstep = -1.d0*dcmplx(0.d0,(height/nz3))
DO i = nz1+nz2+1, nz+1
   z(i) = z(i-1) + zstep
END DO


!CALL compute_g(nz,nblocks,nwand,G,H,xk,wk,z,parent(1:natoms), &
!               dist(:,1:natoms),block_start(1:nblocks),block_dim(1:nblocks),ef)
CALL compute_g(nz,natoms,nwand,G,H,xk,wk,z,parent(1:natoms), &
               dist(:,1:natoms),block_start(1:natoms),block_dim(1:natoms),ef)

allocate(Jexc(natoms,natoms))
allocate(Jorb(natoms,natoms,MAXVAL(block_dim),MAXVAL(block_dim)))
allocate(tmp1(MAXVAL(block_dim),MAXVAL(block_dim)))
Jexc = ZERO
Jorb = ZERO

do ia=1,natoms
   do ja=ia+1,natoms
        ! arrays indexes
        istart = block_start(parent(ia))
        idim   = block_dim(parent(ia))
        iend   = block_start(parent(ia)) + block_dim(parent(ia)) - 1
        jstart = block_start(parent(ja))
        jdim   = block_dim(parent(ja))
        jend   = block_start(parent(ja)) + block_dim(parent(ja)) - 1

        IF(idim .NE. jdim ) THEN
          write(*,*) ia,ja,idim,jdim
          write(*,*) 'leip: Not equal subblocks size'
        END IF

        do iz=1,nz
          zstep = z(iz+1) - z(iz)
          tmp1 = ZEROC
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
  Jexc = (-1.d0/tpi)*Jexc
  Jorb = (-1.d0/tpi)*Jorb

  DO ia = 1, natoms
    DO ja = ia+1, natoms
        ! Compute distance between atoms for pretty output
        !distance = SQRT((tau(1,parent(ia))+newdist(1,ia) - tau(1,parent(ja))-newdist(1,ja))**2+ &
        !                (tau(2,parent(ia))+newdist(2,ia) - tau(2,parent(ja))-newdist(2,ja))**2+ &
        !                (tau(3,parent(ia))+newdist(3,ia) - tau(3,parent(ja))-newdist(3,ja))**2 )
        !
        write(iout,'(/5x,a35,i3,a4,i3,a2,f7.3,a7,i5)') &
              'Exchange interaction between atoms ', ia, 'and', ja, ': ', &
              Jexc(ia,ja)*rytoev*1.0d3, ' meV = ', INT(Jexc(ia,ja)/K_BOLTZMANN_RY)
        !      Jexc(ia,ja)*rytoev*1.0d3, ' meV = ', INT(Jexc(ia,ja)/K_BOLTZMANN_RY), ' K (distance:',distance,')'

        write(iout,'(/7x,a62)') "Orbital exchange interaction matrix J_{i,j,m,n} (in K and meV):"
        DO i=1,block_dim(parent(ia))
          write(iout, '(7x,5i5,8x,5f9.5)' ) INT(Jorb(ia,ja,i,:)/K_BOLTZMANN_RY), Jorb(ia,ja,i,:)*rytoev*1.d3
        END DO

    END DO
  END DO

  write(iout,*)
  write(iout,*) '    Computed orbitals occupations should coincide with your DFT results'
  write(iout,*) '    If they differ significanlty - check your integration contour.'
  DO ia=1,natoms
    DO j=1,nspin
      write(iout,'(/5x,a8,i3,a5,i2)') 'For atom', ia, 'spin', j
      DO i = 1, block_dim(parent(ia))
          tmp2 = ZEROC
          DO iz=1,nz
            zstep = z(iz+1) - z(iz)
            tmp2 = tmp2 + ((-1.d0/pi)*DIMAG(G(iz,ia,ia,i,i,j)*zstep))
          END DO
        write(iout,'(7x,a8,i2,a13,f6.3)') 'Orbital', i, ' occupation: ', DREAL(tmp2)
      END DO
    END DO
  END DO




!clean up
deallocate(H,G,xk,wk,z,delta,Jexc,Jorb,tmp1)
CLOSE(iout)
END

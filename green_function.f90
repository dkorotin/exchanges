SUBROUTINE compute_g(nz,natoms,nblocks,gdim,G,H,z,parent,taunew,block_start,block_dim)

  use parameters, only : dp, tpi
  use general, only : hdim, nkp, nspin, xk, wk, tau, block_atom

  implicit none

  INTEGER, INTENT(IN) :: nz, natoms, nblocks, gdim, parent(natoms), block_start(nblocks), block_dim(nblocks)
  COMPLEX(DP), INTENT(IN) :: H(hdim,hdim,nkp,nspin), z(nz)
  REAL(DP), INTENT(IN) :: taunew(3,natoms)
  COMPLEX(DP), INTENT(OUT) :: G(nz,natoms,natoms,gdim,gdim,nspin)

  INTEGER :: i,j,k, ia, ja, ik, ispin, iz, istart, jstart
  COMPLEX(DP) :: kphase
  COMPLEX(DP), ALLOCATABLE :: Gloc(:,:)

  G = cmplx(0.0,0.0,dp)

  allocate(Gloc(hdim,hdim))
 
  !$OMP PARALLEL PRIVATE(Gloc,kphase,istart,jstart)
  do ispin = 1, nspin
    DO ik = 1, nkp

      !$OMP DO
      DO iz = 1, nz

        CALL compute_gloc(Gloc,H(:,:,ik,ispin),z(iz))

        DO ia = 1, natoms
          DO ja = 1, natoms
            ! exp(i*k*(Ri-Rj))
            kphase = cdexp( 1.d0*DCMPLX(0.d0,1.d0)*tpi*&
                DOT_PRODUCT( xk(:,ik), &
                  ((taunew(:,ia)-tau(:,block_atom(parent(ia))))-(taunew(:,ja)-tau(:,block_atom(parent(ja)))) ) ))

            istart = block_start(parent(ia))-1
            jstart = block_start(parent(ja))-1
                          
            DO i = 1, block_dim(parent(ia))
              DO j = 1, block_dim(parent(ja))
                  G(iz,ia,ja,i,j,ispin) = G(iz,ia,ja,i,j,ispin) + wk(ik)*Gloc(istart+i,jstart+j)*kphase
              END DO
            END DO

          END DO
        END DO

      END DO
      !$OMP END DO

    END DO
  end do
  !$OMP END PARALLEL 

  deallocate(Gloc)

END SUBROUTINE compute_g

subroutine compute_gloc(Gloc,H,z)
  use parameters, only : dp
  use general, only: hdim, efermi, nkp

  implicit none
  
  complex(dp), intent(in) :: H(hdim,hdim), z
  complex(dp), intent(out) :: Gloc(hdim,hdim)
  
  complex(dp) ::  tmp(hdim,hdim)
  integer :: i, j

  Gloc = cmplx(0.0,0.0,dp)

  tmp(:,:) = -1.d0*H(:,:)

  do i=1, hdim
    tmp(i,i) = z + tmp(i,i) + cmplx(efermi,0.d0)
  end do 

  call inverse_complex_matrix(hdim,tmp)
    
  Gloc(:,:) = tmp(:,:)

end subroutine compute_gloc
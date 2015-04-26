SUBROUTINE compute_g(nz,nblocks,gdim,G,H,xk,wk,z,parent,dist,block_start,block_dim,ef)

USE kinds, only : DP,nwan,natoms,nspin,nks,ZEROC,ZERO,tpi
implicit none

INTEGER            :: nz, nblocks, gdim, parent(natoms), block_start(nblocks), block_dim(nblocks)
COMPLEX(DP)        :: H(nwan,nwan,nks,nspin), z(nz+1)
COMPLEX(DP)        :: G(nz,natoms,natoms,gdim,gdim,nspin)
REAL(DP), INTENT(IN)     :: dist(3,natoms), xk(3,nks), wk(nks)

INTEGER                  :: ik, iz, istart, jstart, ia, ja, ispin, i, j
REAL(DP)                 :: v, ek, ef
COMPLEX(DP)              :: kphase

COMPLEX(DP), ALLOCATABLE :: Gloc(:,:)


G = ZEROC
allocate(Gloc(nwan,nwan))

do ispin = 1, nspin
   write(*,*) ispin
 do ik = 1, nks
   ek = ZERO !not needed
   v  = ZERO
   do ia = 1, natoms
     do ja = 1, natoms
        ! exp(i*k*(Ri-Rj))
        ! arg = -1.d0*tpi*DOT_PRODUCT( xk(:,ik), (dist(:,ia)-dist(:,ja)) )
        ! kphase = DCMPLX (COS (arg), -1.d0*SIN (arg) )
        kphase = cdexp( 1.d0*DCMPLX(0.d0,1.d0)*tpi*DOT_PRODUCT( xk(:,ik), (dist(:,ia)-dist(:,ja)) ) )
        istart = block_start(parent(ia))-1
        jstart = block_start(parent(ja))-1

        do iz = 1, nz
          call compute_gloc(Gloc,H(:,:,:,ispin),ik,z(iz),ef)
          do i = 1, block_dim(parent(ia))
            do j = 1, block_dim(parent(ja))
              G(iz,ia,ja,i,j,ispin) = G(iz,ia,ja,i,j,ispin) + wk(ik)*Gloc(istart+i,jstart+j)*kphase
            END DO
          END DO
         END DO !iz
    END DO !ya
  END DO !ia
 END DO !ik
END DO !ispin

END

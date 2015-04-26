SUBROUTINE compute_gloc(Gloc,H,ik,z,ef)

USE kinds, only : DP, ZEROC, nwan, nks
implicit none

INTEGER     :: ik,i,j
COMPLEX(DP) :: H(nwan,nwan,nks), Gloc(nwan,nwan), z, tmp(nwan,nwan)
REAL(DP)    :: ef

  Gloc = ZEROC
  tmp(:,:) = -1.d0*H(:,:,ik)
  DO i=1, nwan
    tmp(i,i) = z + tmp(i,i) + DCMPLX(ef,0.d0)
  END DO

             !write(20,*) 'compute_gloc: H + z' 
             !DO i=1, nwan
             ! write(20,'(16f8.5)') (DREAL(tmp(i,j)),j=1,nwan)
             !END DO
             !write(20,*) 'Diagonal:'
             !write(20,*) (DREAL(tmp(j,j)),j=1,nwan)
             !write(20,*) '---'
             !DO i=1, nwan
             ! write(20,'(16f8.5)') (DIMAG(tmp(i,j)),j=1,nwan)
             !END DO
             !write(20,*) '  '

  CALL inverse_complex_matrix(nwan,tmp)

  Gloc(:,:) = tmp(:,:)
END

SUBROUTINE inverse_complex_matrix(dim,a)

  USE kinds, ONLY : DP
  implicit none

  integer     :: dim, info, ipiv(dim)
  complex(DP) :: a(dim,dim), work(dim)

  call ZGETRF(dim,dim,a,dim,ipiv,info)
  call ZGETRI(dim,a,dim,ipiv,work,dim,info)

END SUBROUTINE

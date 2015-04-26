SUBROUTINE compute_delta(H,wk,delta)

USE kinds, only : DP, nwan, nks, nspin, iout
implicit none
COMPLEX(DP) :: H(nwan,nwan,nks,nspin), delta(nwan,nwan)
REAL(DP)    :: wk(nks)          ! weights of k-points
INTEGER     :: i,ik


    do ik=1, nks  !spin up
      delta(:,:) = delta(:,:) + H(:,:,ik,1)*wk(ik)
    end do
    do ik=1, nks  !spin down
      delta(:,:) = delta(:,:) - H(:,:,ik,2)*wk(ik)
    end do

    !Dima makes Im(delta)=0 by hand
    delta = DREAL(delta)


    write(iout,'(/5x,a36)') 'DEBUG: Real part of Delta (Ry):'
    DO i=1, nwan
      write(iout,'(7x,30f9.5)') DREAL(delta(i,:))
    END DO
    write(iout,'(/5x,a36)') 'DEBUG: Imaginary part of Delta (Ry):'
    DO i=1, nwan
      write(iout,'(7x,30f9.5)') AIMAG(delta(i,:))
    END DO

    delta = DREAL(delta)

END

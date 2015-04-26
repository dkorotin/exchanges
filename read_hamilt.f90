SUBROUTINE read_hamilt(H,xk,wk)

!USE kinds, ONLY: DP, nwan, nks, ZERO, ZEROC, rytoev
USE kinds 

implicit none
REAL(DP)    :: wk(nks)
REAL(DP)    :: xk(3,nks)
COMPLEX(DP) :: Hsum(nwan,nwan,nspin)
COMPLEX(DP) :: H(nwan,nwan,nks,nspin)


!tmp variable
REAL        :: dummy
INTEGER     :: nkp
INTEGER     :: ik,i,j,ispin
REAL(DP)    :: HRe, HIm


!read k-mesh
OPEN (iunhamilt, file = 'klist.inp', status = 'old', form = 'formatted')
   do ik = 1, nks
     read(iunhamilt,*) xk(:,ik)
     read(iunhamilt,*) wk(ik)
     !write(10,*) ik,xk(:,ik)
   end do
close (iunhamilt)

!read hamilt
OPEN (iunhamilt, file = 'hamilt', status = 'old', form = 'formatted')
read (iunhamilt,*) nkp, dummy
if(dummy /= nwan) STOP 'nwan and hamilt dimentions are not equal'
!if(nkp /= nks)    STOP 'nks and number of k-points in hamilt are not equal'

  HRe  = ZERO
  HIm  = ZERO
  H    = ZEROC
  Hsum = ZEROC


  do ispin = 1, nspin
   do ik = 1, nks
    read(iunhamilt,*) dummy
    do i = 1, nwan
      do j=1, nwan
        read(iunhamilt,*) HRe, HIm
        H(i,j,ik,ispin) = dcmplx(HRe,HIm)/rytoev
      end do
    end do
    !wk(ik) = 1.0/nks
    Hsum(:,:,ispin) = Hsum(:,:,ispin) + wk(ik)*H(:,:,ik,ispin)
   end do !ik
  end do !ispin

  write(iout,*)
  write(iout,*) 'WARNING: Note that wk(k) are defined ad hoc'
  write(iout,*) 'WARNING: We suppose that first half of k-points is for spin 1'
  write(iout,*) '         and second half is for spin 2'
  write(iout,*)

    do ispin = 1, nspin
      write(iout,*) 'Real part of H(0) for spin', ispin
      do i=1, nwan
        write(iout,'(5x,80f9.5)') (DREAL(Hsum(i,j,ispin)),j=1,nwan )
      end do
      write(iout,*) 'Imagianry part of H(0) for spin', ispin
      do i=1, nwan
        write(iout,'(5x,80f9.5)') (DIMAG(Hsum(i,j,ispin)),j=1,nwan )
      end do
    end do

close (iunhamilt)

END

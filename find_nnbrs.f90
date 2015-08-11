! Copyright (C) Dmitry Korotin dmitry@korotin.name

subroutine find_nnbrs(natoms,tau,cell,aoi,maxdistance,nnnbrs,taunew,parent)
	
	use parameters, only : dp, maxnnbrs
	implicit none

	integer, intent(in) :: natoms, & ! Number of atoms in initial cell
						   					 aoi ! atom of interest. We will search neighbours of this atom
	real(dp), intent(in) :: tau(3,natoms), & ! positions of atoms in the initial cell
													cell(3,3), & ! cell vectors row-wise (in units of alat)
													maxdistance ! max distance betwee aoi and neighbouring atom

	real(dp), intent(out) :: taunew(3,maxnnbrs) ! nnbrhds positions
	integer, intent(out) :: nnnbrs, & ! total number of nearest neighbours
													parent(maxnnbrs) ! arent atom for each neighbour

  integer :: i, j, k, ntransl, iatom, new_a, parenttmp
  real(dp) v(3), new_pos(3), distance, dist(maxnnbrs), distt, tautmp(3)
  logical :: have_atom_already

  taunew = -1000.0
  parent = -1
  new_a = 0

  ntransl = ceiling(distance) + 2
  do i = -ntransl, ntransl
    do j = -ntransl, ntransl
      do k = -ntransl, ntransl
        ! translation vector:
        v = real(i,dp)*cell(:,1) + real(j,dp)*cell(:,2) + real(k,dp)*cell(:,3)

        do iatom = 1, natoms
          new_pos = tau(:,iatom) + v
          distance = SQRT((new_pos(1)-tau(1,aoi))**2 + &
                          (new_pos(2)-tau(2,aoi))**2 + &
                          (new_pos(3)-tau(3,aoi))**2 )
          
          if( distance .le. maxdistance ) then
            call haa(taunew(:,:),new_pos,have_atom_already)
            if( .not. have_atom_already ) then
              new_a=new_a+1
              taunew(:,new_a) = new_pos
              parent(new_a) = iatom
              dist(new_a) = distance
              IF (new_a .eq. maxnnbrs) &
                stop 'Error find_nnbrs: max number of nearest neighbours reached'
            end if
          end if

        end do
      end do  
    end do
  end do

  nnnbrs = new_a

  ! Now let's sort atoms depending on distance
  ! Trivial buuble algorithm
  do j=1, nnnbrs-1
    do i=1, nnnbrs-j
      if(dist(i) > dist(i+1)) then
        ! swap
        distt = dist(i)
        dist(i) = dist(i+1)
        dist(i+1) = distt

        tautmp = taunew(:,i)
        taunew(:,i) = taunew(:,i+1)
        taunew(:,i+1) = tautmp

        parenttmp = parent(i)
        parent(i) = parent(i+1)
        parent(i+1) = parenttmp
      end if
    end do
  end do

end subroutine find_nnbrs

subroutine haa(storage, new_pos,h)

! Determines if storage contains new_pos already or not
! If yes: h = .true.

  use parameters
  
  implicit none
  
  real(dp), intent(in) :: storage(3,maxnnbrs), new_pos(3)
  logical, intent(out) :: h
  
  integer :: i
  real(dp),parameter :: eps=1.d-5 
  
  h = .false.
  do i=1, maxnnbrs
    if((ABS(storage(1,i)-new_pos(1)) .LE. eps) &
      .AND. (ABS(storage(2,i)-new_pos(2)) .LE. eps) &
      .AND. (ABS(storage(3,i)-new_pos(3)) .LE. eps)) h=.true.
  end do
  
end subroutine
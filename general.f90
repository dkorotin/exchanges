module general
 
  use parameters, only : dp 	 
  implicit none
  save

  !Crystal cell
  real(dp) :: alat, cell(3,3)
  real(dp), allocatable :: tau(:,:)
  integer :: natoms
  character(len=3), allocatable :: atomlabel(:)

end module general
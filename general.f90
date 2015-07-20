! Copyright (C) Dmitry Korotin dmitry@korotin.name

module general
 
  use parameters, only : dp 	 
  implicit none
  save

  !Crystal cell
  real(dp) :: alat, cell(3,3)
  real(dp), allocatable :: tau(:,:)
  integer :: natoms
  character(len=3), allocatable :: atomlabel(:)

  !Hamiltonian
  integer :: &
  			nkp, &  ! number of k-points
  			hdim, & ! hamiltonian dimention
  			nspin	!

  real(dp), allocatable :: &
  			wk(:,:), & ! k-points weight (nspin,nkp)
  			xk(:,:,:) ! k-points coordinates in 2pi/a units (3,nkp,nspin)
  complex(dp), allocatable :: h(:,:,:,:) ! Hamiltonian matrix (hdim,hdim,nspin)

  
  !Basis
  integer :: nblocks
  integer, allocatable :: block_atom(:), block_dim(:), block_orbitals(:,:)
  character, allocatable :: block_l(:)

  integer, allocatable :: block_start(:)

  !Other
  real(dp) :: efermi


  ! That is how the atomic orbitals are encoded in the system.am file
  character(len=11) :: orbitals(16)
  data orbitals/'s', 'y', 'z', 'x', 'xy', 'yz', '3z^2-1', 'xz', 'x^2-y^2', 'y(3x^2-y^2)', 'xyz', &
  				'y(5z^2-1)', 'z(5z^2-3)', 'x(5z^2-1)', 'z(x^2-y^2)', 'x(3y^2-x^2)'/


end module general
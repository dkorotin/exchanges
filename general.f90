! Copyright (C) Dmitry Korotin dmitry@korotin.name

module general
 
  use parameters, only : dp, maxnnbrs 	 
  implicit none
  save

  !Crystal cell
  real(dp) :: alat, cell(3,3), inv_cell(3,3)
  real(dp), allocatable :: tau(:,:)
  integer :: natoms
  character(len=3), allocatable :: atomlabel(:)

  !Hamiltonian
  integer :: &
  			nkp, &  ! number of k-points
  			hdim, & ! hamiltonian dimention
  			nspin	!

  real(dp), allocatable :: &
  			wk(:), & ! k-points weight (nkp)
  			xk(:,:) ! k-points coordinates in 2pi/a units (3,nkp)
  complex(dp), allocatable :: h(:,:,:,:) ! Hamiltonian matrix (hdim,hdim,nkp,nspin)

  
  !Basis
  integer :: nblocks
  integer, allocatable :: block_atom(:), block_dim(:), block_orbitals(:,:)
  character, allocatable :: block_l(:)

  integer, allocatable :: block_start(:)

  !Other
  real(dp) :: efermi

  ! Nearest neighbours
  real(dp) :: &
             taunew(3,maxnnbrs) ! positions of nearest neighbours
  integer :: nnnbrs, & ! total number of nearest neighbours
             parent(maxnnbrs) ! parent atom for every neighbour


  ! That is how the atomic orbitals are encoded in the system.am file
  character(len=11) :: orbitals(16)
  data orbitals/'s', 'y', 'z', 'x', 'xy', 'yz', '3z^2-1', 'xz', 'x^2-y^2', 'y(3x^2-y^2)', 'xyz', &
  				'y(5z^2-1)', 'z(5z^2-3)', 'x(5z^2-1)', 'z(x^2-y^2)', 'x(3y^2-x^2)'/

  contains
  subroutine inverse_complex_matrix(dim,a)
  
    use parameters, only : dp
    implicit none

    integer :: dim, info, ipiv(dim)
    complex(dp) :: a(dim,dim), work(dim)

    call ZGETRF(dim,dim,a,dim,ipiv,info)
    if(info /= 0) stop "inverse_complex_matrix Error in ZGETRF"
    call ZGETRI(dim,a,dim,ipiv,work,dim,info)
    if(info /= 0) stop "inverse_complex_matrix Error in ZGETRI"

  end subroutine
  
  subroutine read_crystal()

    use iomodule
    use parameters, only : dp
    
    implicit none
    integer :: i,j
    character(len=3) :: dummy
    character :: l_symbol

    call open_input_file(iunsystem,'system.am')

    call find_section(iunsystem,'&cell')
    read(iunsystem,*) alat
    
    do i = 1,3
      read(iunsystem,*) cell(:,i)
    end do

    inv_cell = matinv3(cell)
    
    call find_section(iunsystem,'&atoms')
    read(iunsystem,*) natoms
    allocate( tau(3,natoms) )
    allocate( atomlabel(natoms) )

    do i=1, natoms
      read(iunsystem,*) atomlabel(i), tau(:,i)
    end do

    call find_section(iunsystem,'&basis')
    read(iunsystem,*) i, nblocks
    if(i .ne. hdim) stop 'Basis dimention in system.am and hamilt.am are unequal'

    allocate( block_atom(nblocks) )
    allocate( block_l(nblocks) )
    allocate( block_dim(nblocks) )
    allocate( block_orbitals(nblocks,7) )
    allocate( block_start(nblocks) )
    block_orbitals = 0

    do i = 1, nblocks
      read(iunsystem,*) dummy, block_atom(i), block_l(i), block_dim(i), block_start(i), block_orbitals(i,1:block_dim(i))
    end do

    call find_section(iunsystem,'&efermi')
    read(iunsystem,*) efermi

    close(iunsystem)

  end subroutine read_crystal

  subroutine read_hamilt()

  use iomodule
  use parameters, only : dp
  
  implicit none
  integer :: i,j,ispin, ik
  real(dp) :: Hre, Him

  call open_input_file(iunhamilt,'hamilt.am')

  call find_section(iunhamilt,'&nspin')
  read(iunhamilt,*) nspin
  if (nspin .ne. 2) stop 'A spin-polarized hamiltonian is neccessary'

  call find_section(iunhamilt,'&nkp')
  read(iunhamilt,*) nkp

  call find_section(iunhamilt,'&dim')
  read(iunhamilt,*) hdim

  allocate( h(hdim,hdim,nkp,nspin) )
  allocate( wk(nkp) )
  allocate( xk(3,nkp) )

  h = cmplx(0.0,0.0)
  wk = 0.0
  xk = 0.0

  call find_section(iunhamilt,'&kpoints')
  do ik=1, nkp
    read(iunhamilt,*) wk(ik), xk(:,ik)
  end do
      
  call find_section(iunhamilt,'&hamiltonian')
  
  do ispin = 1, nspin
    do ik = 1, nkp
      do i=1, hdim
        do j=i, hdim
          read(iunhamilt,*) Hre, Him
          h(i,j,ik,ispin) = cmplx(Hre,Him)
          h(j,i,ik,ispin) = dconjg( h(i,j,ik,ispin) )
        end do
      end do
    end do
  end do
  
  close(iunhamilt)

end subroutine read_hamilt

subroutine atoms_list(mode,distance,atom_of_interest,l_of_interest)

  use iomodule
  use parameters, only : dp, maxnnbrs
  
  implicit none
  character(len=10), intent(in) :: mode
  real(dp), intent(in) :: distance
  character(len=1), intent(in) :: l_of_interest
  integer, intent(in) :: atom_of_interest

  ! temp storage
  real(dp) :: taunew_(3,maxnnbrs), vect(3), d
  integer :: parent_(maxnnbrs), nnnbrs_

  integer :: i,j, iblock, nvect, index
  logical :: have_atom_already

  character(len=40) :: currentline

  nnnbrs = 0
  parent = -1
  taunew = -1000

  !call find_nnbrs(natoms,tau,cell,block_atom(1),distance,nnnbrs_,taunew_,parent_)
  call find_nnbrs(natoms,tau,cell,atom_of_interest,distance,nnnbrs_,taunew_,parent_)

  select case ( trim(mode) )
    case ('list') ! list mode
      ! a primitive way to find ATOMS_LIST string in stdin
      currentline = ''
      do while( ios .ge. 0 .AND. trim(currentline) .ne. 'ATOMS_LIST' )
        read(stdin, fmt="(a40)", iostat=ios ) currentline
        if ( ios .gt. 0 ) stop "Can not find ATOMS_LIST in stdin"
      end do
      !call find_section(stdin,'ATOMS_LIST')
      read(stdin,*) nvect

      do i = 1, nvect
        read(stdin,*) vect(1:3)
        do j = 1, nnnbrs_
          call haa(taunew_, vect, have_atom_already, index)
          if( .not. have_atom_already ) then
            write(stdout,'(/,"Error: Can not find atom in position", 3f9.5, " within sphere of ",f9.5," with ",a2," orbitals"/)') &
              vect,distance,l_of_interest
            stop ' '
          else
            nnnbrs = nnnbrs + 1
            parent(nnnbrs) = parent_(index)
            taunew(:,nnnbrs) = vect
            exit
          end if
        end do
      end do
      ! reset arrays
      nnnbrs_ = nnnbrs
      parent_ = parent
      taunew_ = taunew

    case default ! 'distance mode'

      ! We want to consider all hamilt atoms anyway
      do i = 1, natoms
        call haa(taunew_, tau(:,i), have_atom_already, index)
        if( .not. have_atom_already ) then
          nnnbrs_ = nnnbrs_+1
          taunew_(:,nnnbrs_) = tau(:,i)
          parent_(nnnbrs_) = i
        end if
      end do

  end select

  nnnbrs = 0
  parent = -1
  taunew = -1000

  !filter atoms by l
  do i = 1, nnnbrs_
    do iblock = 1, nblocks
      if( block_atom(iblock) .eq. parent_(i) .and. block_l(iblock) .eq. l_of_interest ) then
        nnnbrs = nnnbrs + 1
        parent(nnnbrs) = iblock
        taunew(:,nnnbrs) = taunew_(:,i)
      end if
    end do 
  end do

end subroutine atoms_list

pure function matinv3(A) result(B)
  !! Performs a direct calculation of the inverse of a 3×3 matrix.
  real(dp), intent(in) :: A(3,3)   !! Matrix
  real(dp)             :: B(3,3)   !! Inverse matrix
  real(dp)             :: detinv

  ! Calculate the inverse determinant of the matrix
  detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
            - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
            + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

  ! Calculate the inverse of the matrix
  B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
  B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
  B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
  B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
  B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
  B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
  B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
  B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
  B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
end function

end module general
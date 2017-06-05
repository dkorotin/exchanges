! Copyright (C) Dmitry Korotin dmitry@korotin.name

module iomodule

  implicit none
  save

  integer, parameter :: &
    stdin = 5, &
    stdout = 6, &
    iunhamilt = 31, & ! hamiltonian
    iunsystem = 32    ! data about crystal

  integer :: globalhash = 0
  integer :: ios = 0
  integer :: iverbosity = 0 ! verbosity of the output

  contains
  
  subroutine find_section(unit, sect)
    
    implicit none
    integer, intent(in) :: unit
    character(len=*), intent(in) :: sect

    character(40) :: currentline = '', sec

    sec = trim(sect)
    rewind(unit)
    
    do while( ios .ge. 0 .AND. trim(currentline) .ne. sec )
      read(unit=unit, fmt="(a40)", iostat=ios ) currentline
      if ( ios .gt. 0 ) stop "Read error in file unit"

    end do

  end subroutine find_section

  subroutine open_input_file(unit,filename)

    integer, intent(in) :: unit
    character(len=*), intent(in) :: filename

    character(40) :: hashline
    integer :: hash

    open(unit=unit, file=trim(filename), iostat=ios, status="old", action="read")
    if ( ios /= 0 ) then
      write(stdout,*) "Error opening file "//filename
      stop 
    end if

    !here I want to check the hash
    call find_section(unit,'&hash')
    read(unit,*) hashline
    read(hashline,*) hash
    if (globalhash .eq. 0) then
      globalhash = hash
    else
      if ( globalhash .ne. hash ) then
        write(stdout,*) 'Error: Input files have different hashes', globalhash, hash
        stop
      end if
    end if
    
  end subroutine open_input_file

  subroutine output_matrix_by_blocks(dim,matrix)
    use parameters, only : dp
    use general, only : nblocks, block_start, block_dim
    implicit none

    integer, intent(in) :: dim
    real(dp), intent(in) :: matrix(dim,dim)
    integer :: iblock, i, j, block_end

    character :: fmt

    do iblock = 1, nblocks

      block_end = block_start(iblock) + block_dim(iblock) - 1
      write(stdout,'(7x,a5,i3,a1)') 'Block', iblock, ':'
      write(fmt,'(i1)') block_dim(iblock)
      do i=block_start(iblock), block_end
        write(stdout,'(7x,'//fmt//'f9.5)') (matrix(i,j), j=block_start(iblock), block_end)
      end do
    
    end do

  end subroutine output_matrix_by_blocks

end module iomodule

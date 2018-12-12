! Copyright (C) Dmitry Korotin dmitry@korotin.name

program atoms_list_prog

  use general
  use parameters
  use iomodule

  implicit none

  integer :: i, j, time_start, time_end, count_rate, nz, ia, ja, idim, jdim, iz, istart, jstart, iend, jend
  character(len=3) :: fmt='   ' ! is used for pretty output only
  character(len=10) :: arg ! cli argument
  real(dp) :: pos_delta(3)
  integer :: atom_of_interest
  character(len=1) :: l
  
  real(dp) :: distance ! distance for the nearest neighbours search

  call get_command_argument(1, arg)
  if (len_trim(arg) == 0) then
    distance = 8.d-1
  else
    read(arg, *) distance
  end if
  
  l = 'd'
  atom_of_interest = -100
  
  call system_clock(time_start,count_rate)

  call read_hamilt()
  call read_crystal()

  if(atom_of_interest == -100) atom_of_interest = block_atom(1)
 
  write(stdout,'(/,5x,a20,f13.9,a5)') 'Cell constant (alat):', alat, 'Bohr'
  write(stdout,'(5x,20a)') 'Cell vectors (rows):'
  do i = 1, 3
    write(stdout,'(7x,3f9.5)') cell(:,i)
  end do
  write(stdout,'(5x,a26)') 'Atoms of the initial cell:'
  do i = 1, natoms
    write(stdout,'(7x,a3,x,3f9.5)') atomlabel(i), tau(:,i)
  end do
    
  call atoms_list('distance  ',distance,atom_of_interest,l)

  ! Pretty output
  write(stdout,'(/5x,a17,i4,a8)') 'We will consider ', nnnbrs,' atoms: '

  do i = 1, nnnbrs
    pos_delta =   taunew(:,i)-tau(:,atom_of_interest)
    write(stdout,'(5x,i3,a2,a3,x,3f9.5,3x,a12,f9.5,a7)') i, ': ', &
      atomlabel( block_atom( parent(i) ) ), taunew(:,i), '( distance =', norm2(pos_delta), ' alat )'
    write(stdout,'(10x,a37,3f8.4,a3,3f8.4)') 'Connecting vector (Cart. / crystal): ', pos_delta, &
                                             ' / ', MATMUL( pos_delta, TRANSPOSE(inv_cell)  )
  end do
        
  call clear()

  call system_clock(time_end,count_rate)

  write(stdout,'(/5x,a15,i5,a8)') 'Execution time:', (time_end - time_start)/count_rate, 'seconds'

end program atoms_list_prog

subroutine clear()

  use general

  if( allocated(tau) ) deallocate(tau)
  if( allocated(atomlabel) ) deallocate(atomlabel)
  if( allocated(block_atom) ) deallocate( block_atom )
  if( allocated(block_l) ) deallocate( block_l )
  if( allocated(block_dim) ) deallocate( block_dim )
  if( allocated(block_orbitals) ) deallocate( block_orbitals )
  if( allocated(block_start) ) deallocate(block_start)

end subroutine clear

! Copyright (C) Dmitry Korotin dmitry@korotin.name

module parameters
  
  implicit none
  save

  integer, parameter :: &
    dp = selected_real_kind(14,200), &
    maxnnbrs = 100 ! maximum number of nearest neighbours
  

end module parameters
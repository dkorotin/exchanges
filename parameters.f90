! Copyright (C) Dmitry Korotin dmitry@korotin.name

module parameters
  
  implicit none
  save

  integer, parameter :: &
    dp = selected_real_kind(14,200), &
    maxnnbrs = 300 ! maximum number of nearest neighbours
  
  real(dp), parameter :: pi               = 3.14159265358979323846_dp
  real(dp), parameter :: tpi              = 2.0_dp * pi
  real(dp), parameter :: kb_ev = 8.6173324d-5
  

end module parameters
!------------------------------------------------------------------------------!
    MODULE kinds
!------------------------------------------------------------------------------!

      IMPLICIT NONE
      SAVE
! ... kind definitions
      INTEGER, PARAMETER     :: DP               = selected_real_kind(14,200)
      INTEGER, PARAMETER     :: nmaxatoms        = 20
      REAL(DP), PARAMETER    :: ELECTRONVOLT_SI  = 1.602176487E-19_DP    
      REAL(DP), PARAMETER    :: HARTREE_SI       = 4.35974394E-18_DP   
      REAL(DP), PARAMETER    :: AUTOEV           = HARTREE_SI / ELECTRONVOLT_SI
      REAL(DP), PARAMETER    :: RYTOEV           = AUTOEV / 2.0_DP
      REAL(DP), PARAMETER    :: pi               = 3.14159265358979323846_DP
      REAL(DP), PARAMETER    :: tpi              = 2.0_DP * pi
      REAL(DP), PARAMETER    :: ZERO             = 0.0_DP
      REAL(DP), PARAMETER    :: RYDBERG_SI       = HARTREE_SI/2.0_DP   ! J
      REAL(DP), PARAMETER    :: K_BOLTZMANN_SI   = 1.3806504E-23_DP    ! J K^-1 
      REAL(DP), PARAMETER    :: K_BOLTZMANN_RY   = K_BOLTZMANN_SI / RYDBERG_SI
      COMPLEX(DP), PARAMETER :: ZEROC            = (0.0_DP, 0.0_DP)

      !files
      INTEGER, PARAMETER     :: iout             = 20
      INTEGER, PARAMETER     :: igeometry        = 30
      INTEGER, PARAMETER     :: iunhamilt        = 40

      PUBLIC ::  DP

 
      INTEGER  :: &
                 nwan,             &! number of wannier functions
                 natoms,           &! number of atoms
                 nspin,            &! number of spin
                 nks                ! number of k-points
!
!------------------------------------------------------------------------------!
!
    CONTAINS
!
!------------------------------------------------------------------------------!
!
!
!------------------------------------------------------------------------------!
    END MODULE kinds
!------------------------------------------------------------------------------!

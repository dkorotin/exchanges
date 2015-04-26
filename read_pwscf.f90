SUBROUTINE read_pwscf(nz1,nz2,nz3,height,emin,emax,ef,nwand,natoms)

USE kinds, only :   DP, nwan, rytoev          !nwan - total number of Wannier functions (not only d)

implicit none
INTEGER         ::   nwand                     !number of localized orbitals per site
INTEGER         ::   natoms                    !number of atoms with localized orbitals
INTEGER         ::   nz1, nz2, nz3             !number of points in contour for energy integration  
REAL(DP)        ::   emin, emax, height        !energy points for the integration
REAL(DP)        ::   ef                        !fermi energy

NAMELIST /leip/   nwan, nz1, nz2, nz3, height, emin, emax, ef, nwand, natoms

!default values
nz1    = 250
nz2    = 450
nz3    = 250
emin   = -15.0           !in Ry
emax   = 0.1             !in Ry
height = 0.2             !in Ry
ef     = 0.865150783431619 
nwand  = 5
natoms = 2

READ (5, leip)

!transform to eV
emin = emin/rytoev
emax = emax/rytoev
height = height/rytoev

END

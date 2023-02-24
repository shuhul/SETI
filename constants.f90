MODULE constants      
  implicit none
  integer, parameter :: DP = SELECTED_REAL_KIND(14)
  real(KIND=DP), parameter :: pi  = 3.14159265358979324
  real(KIND=DP), parameter :: kB  = 1.3807d-16 ! Boltz. constant [erg/K]
  real(KIND=DP), parameter :: mp  = 1.6726d-24 ! Mass of proton [g]
  real(KIND=DP), parameter :: GG  = 6.672d-8 ! Gravitational constant
  real(KIND=DP), parameter :: ss = 5.6703d-5 ! Stefan-Boltz const  [erg/cm^2/K^4/s]
  real(KIND=DP), parameter :: au=1.496d13    ! Astronomical unit [cm]
  real(KIND=DP), parameter :: solm=1.99d33   ! Solar mass [g]
  real(KIND=DP), parameter :: solr=6.99d10  ! Solar radius [cm]
END MODULE constants

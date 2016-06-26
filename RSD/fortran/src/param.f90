module param
  implicit none

  real*8,parameter :: om_m0=0.26
  real*8,parameter :: om_v0=0.74
  real*8,parameter :: h0=0.70

  real*8,parameter :: redshift=3.0
  real*8,parameter :: Hubble=100.0*sqrt(om_v0+om_m0*(1.0+redshift)**3) ! km/s/(Mpc/h)

  real*8,parameter :: pi=3.14159265359
  real*8,parameter :: c=2.998e5 ! km/s

  ! model parameters
  real*8 :: r0, slope_xi                   ! [cMpc/h], [-] 
  real*8 :: v_inflow, r_inflow, slope_v12  ! [km/s], [cMpc/h], [-] 
  real*8 :: v_outflow, r_outflow           ! [km/s], [cMpc/h]
  real*8 :: s12_0                          ! [km/s]
  real*8 :: r_eq, betaN                    ! [cMpc/h], [-]

  ! astrophysical model parameters
  real*8,parameter :: T_abs=60000.d0                     ! [K] temperature of absorber
  real*8,parameter :: a=0.0472d0/sqrt(T_abs)            ! a-parameter of Voigt profile
  real*8,parameter :: doppler_width=1.057d9*sqrt(T_abs) ! [Hz]
  real*8,parameter :: cross_sec=0.011d0                  ! [cm2*Hz]
  real*8,parameter :: freq_lya=2.466d15                  ! [Hz] lya resonance frequency


end module param

module array
  implicit none
  integer,parameter :: N=100                           ! length of array
  real*8, dimension(:), allocatable :: xi_array, v_array
  real*8,parameter :: vmin=-1000.d0              ! km/s, min/max of velocity-space
  real*8,parameter :: vmax=+1000.d0              !       correlation function
  real*8,parameter :: dv=(vmax-vmin)/N

  integer,parameter :: M=100
  real*8,parameter ::  NHI_min=10.0**12.0         ! [1/cm2] log(NHI) minimum
  real*8,parameter ::  NHI_max=10.0**21.55        ! [1/cm2] log(NHI) maximum
  real*8,parameter :: dlogNHI=(log10(NHI_max)-log10(NHI_min))/M
  real*8, dimension(:), allocatable :: dtau_array, NHI_array


end module array

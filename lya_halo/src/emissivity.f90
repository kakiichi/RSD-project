function emissivity(r)
  use param
  use array
  implicit none
  real*8 :: emissivity, r
  real*8 :: lnNHI_min,lnNHI_max
  integer :: i
  real*8 :: Lv, xi
  ! quadpack parameters
  real*8, parameter :: epsabs=0.0d0
  real*8, parameter :: epsrel=1.d-6
  real*8 :: result,abserr
  integer :: neval,ier

  external integrand_emissivity

  ! pre-compute the velocity-integrated luminosity
  allocate(Lv_array(M))
  allocate(NHI_array(M))
  do i=1,M
     NHI_array(i)=10.0**( log10(NHI_min)+(i-0.5)*dlogNHI )
     Lv_array(i)=Lv(r,NHI_array(i))
  end do
 
  lnNHI_min=log(NHI_array(2))
  lnNHI_max=log(NHI_array(M-1))

  call qags(integrand_emissivity, lnNHI_min, lnNHI_max, epsabs, epsrel, result, abserr, neval, ier) 
  emissivity=(1.d0+xi(r))/(4.d0*pi*r*r)*result ! erg/s/(cMpc/h)^3 comoving unit

  ! covert it to proper unit [erg/s/cm^3]
  emissivity=emissivity*(1.d0+redshift)**3*(h0/Mpc2cm)**3

  deallocate(Lv_array)
  deallocate(NHI_array)

  return
end function emissivity

function integrand_emissivity(lnNHI)
  use param
  use array
  implicit none
  real*8 :: integrand_emissivity, lnNHI, NHI
  real*8 :: Lv_interp, drdz
  real*8 :: dNdNHIdz
  integer,parameter :: dim=1

  NHI=exp(lnNHI)
  call interp_linear(dim,M,NHI_array,Lv_array,1,NHI,Lv_interp)
  drdz=c/Hubble ! cMpc/h
  integrand_emissivity=NHI*dNdNHIdz(NHI)* Lv_interp / drdz ! erg/s/(cMpc/h)

  return
end function integrand_emissivity


function Lv(r,NHI)
  implicit none
  real*8 :: Lv, r, NHI
  real*8 :: vmin,vmax
  ! quadpack parameters
  real*8, parameter :: epsabs=0.0d0
  real*8, parameter :: epsrel=1.d-6
  real*8 :: result,abserr
  integer :: neval,ier

  real*8 :: r_share, NHI_share
  common /share_Lv/ r_share, NHI_share
  external integrand_Lv

  r_share = r
  NHI_share = NHI
  vmin=-1000.d0  ! km/s
  vmax=+1000.d0  ! km/s
  call qags(integrand_Lv, vmin, vmax, epsabs, epsrel, result, abserr, neval, ier) 
  Lv=result
  
  return
end function Lv

function integrand_Lv(v)
  use param
  implicit none
  real*8 :: integrand_Lv, v
  real*8 :: vabs, Labs
  real*8 :: v12, s12
 
  real*8 :: r_share, NHI_share
  common /share_Lv/ r_share, NHI_share

  vabs=Hubble*r_share/(1.d0+redshift) + v
  integrand_Lv=Labs(vabs,NHI_share) & 
       *(1.d0/(dsqrt(2.d0*pi)*s12(r_share)))*exp(-(v-v12(r_share))**2/(2.d0*(s12(r_share))**2))

  return
end function integrand_Lv



! absorber's Lya luminosity [erg/s]
function Labs(vabs,NHI)
  implicit none
  real*8 :: Labs, vabs, NHI
  real*8 :: vmin, vmax
  ! quadpack parameters
  real*8, parameter :: epsabs=0.0d0
  real*8, parameter :: epsrel=1.d-6
  real*8 :: result,abserr
  integer :: neval,ier

  real*8 :: vabs_share, NHI_share
  common /share_Labs/ vabs_share, NHI_share
  external integrand_La

  vabs_share = vabs
  NHI_share = NHI
  vmin=-1000.d0  ! km/s
  vmax=+1000.d0  ! km/s
  call qags(integrand_La, vmin, vmax, epsabs, epsrel, result, abserr, neval, ier) 
  Labs=result

  return
end function Labs

function integrand_La(v_e)
  use param
  implicit none
  real*8 :: integrand_La, v_e
  real*8 :: x_inj
  real*8 :: optdpt, intrinsic_profile

  real*8 :: vabs_share, NHI_share
  common /share_Labs/ vabs_share, NHI_share

  x_inj = -freq_lya/doppler_width*( (vabs_share+v_e)/c )
  integrand_La = ( 1.d0-exp(-optdpt(x_inj,NHI_share)) ) * intrinsic_profile(v_e)

  return
end function integrand_La

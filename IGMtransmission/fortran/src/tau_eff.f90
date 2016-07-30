function tau_eff(vel)
  use array
  implicit none
  real*8 :: tau_eff, vel
  real*8 :: xi_v
  integer :: i,j
  ! quadpack parameters
  real*8, parameter :: epsabs=0.0d0
  real*8, parameter :: epsrel=1.d-6
  real*8 :: result,abserr
  integer :: neval,ier
  real*8 :: lnNHI_min,lnNHI_max

  real*8 :: dNdNHIdz ! test, can remove after test

  real*8 :: dtau
  external tau_eff_integrand

  ! pre-compute the velocity-space correlation function
  ! module array holds: xi(:), v_para(:)
  allocate(xi_array(N))
  allocate(v_array(N))
  do i=1,N
     v_array(i)=vmin+(i-0.5)*dv
     xi_array(i)=xi_v(v_array(i))
  end do

  ! compute the integration over the velocity-space correlation function
  allocate(dtau_array(M))
  allocate(NHI_array(M))
  do j=1,M
     NHI_array(j)=10.0**( log10(NHI_min)+(j-0.5)*dlogNHI )
     dtau_array(j)=dtau(vel,NHI_array(j))
  end do

  lnNHI_min=log(NHI_array(2))
  lnNHI_max=log(NHI_array(M-1))


  call qags(tau_eff_integrand, lnNHI_min, lnNHI_max, epsabs, epsrel, tau_eff, abserr, neval, ier) 

  deallocate(xi_array,v_array)
  deallocate(dtau_array,NHI_array)

  return
end function tau_eff


function tau_eff_integrand(lnNHI)
  use param
  use array
  implicit none

  real*8 :: tau_eff_integrand, lnNHI, NHI
  real*8 :: dNdNHIdz
  real*8 :: dtau_int
  integer,parameter :: dim=1

  NHI=exp(lnNHI)

  call interp_linear(dim,M,NHI_array,dtau_array,1,NHI,dtau_int)

  tau_eff_integrand=NHI*dNdNHIdz(NHI)*dtau_int*Hubble/c

  return
end function tau_eff_integrand

function dtau(vel,NHI)
  use param
  use array
  implicit none
  real*8 :: dtau, vel, NHI
  ! quadpack parameters
  real*8, parameter :: epsabs=0.0d0
  real*8, parameter :: epsrel=1.d-6
  real*8 :: result,abserr
  integer :: neval,ier

  real*8 :: vel_share, NHI_share
  common /share_dtau/ vel_share, NHI_share
  external dtau_integrand

  vel_share=vel
  NHI_share=NHI
  call qags(dtau_integrand, vmin+dv, vmax-dv, epsabs, epsrel, dtau, abserr, neval, ier)
  
  return
end function dtau

function dtau_integrand(v_para)
  use param
  use array
  implicit none
  integer,parameter :: dim=1
  real*8 :: dtau_integrand, v_para
  real*8 :: x
  real*8 :: xi_velocity,optdpt

  real*8 :: vel, NHI
  real*8 :: vel_share, NHI_share
  common /share_dtau/ vel_share, NHI_share

  vel=vel_share
  NHI=NHI_share

  x=-freq_lya/doppler_width*(vel+v_para)/c

  call interp_linear(dim,N,v_array,xi_array,1,v_para,xi_velocity)
  
  dtau_integrand=(1.d0+xi_velocity)*(1.d0-exp(-optdpt(x,NHI))) * (1.d0+redshift)/Hubble
  
  return
end function dtau_integrand


function xi_v(v_para,r_perp)
  implicit none
  real*8 :: xi_v,v_para,r_perp
  real*8 :: r_min,r_max ! [cMpc/h] min & max limits of the integral
  ! quadpack parameters
  real*8, parameter :: epsabs=0.0d0
  real*8, parameter :: epsrel=1.d-6
  real*8 :: result,abserr
  integer :: neval,ier

  real*8 :: v_para_share, r_perp_share
  common /share_xi/ v_para_share, r_perp_share
  external integrand_xi

  v_para_share=v_para
  r_perp_share=r_perp

  r_min=-1000.d0
  r_max=+1000.d0
  call qags(integrand_xi, r_min, r_max, epsabs, epsrel, result, abserr, neval, ier)
  xi_v=result-1.d0

  return
end function xi_v


function integrand_xi(r_para)
  use param
  implicit none
  real*8 :: integrand_xi, r_para
  real*8 :: r, mu, f_v12
  real*8 :: xi, v12, s12, phr

  real*8 :: v_para, r_perp
  real*8 :: v_para_share, r_perp_share
  common /share_xi/ v_para_share, r_perp_share

  v_para=v_para_share
  r_perp=r_perp_share

  r=sqrt(r_para**2+r_perp**2)
  mu=r_para/r

  ! gaussian streaming model
  f_v12=1.d0/(dsqrt(2.d0*pi)*s12(r)) * &
        exp( -(v_para-Hubble/(1.d0+redshift)*r_para-mu*v12(r))**2/(2.d0*(s12(r))**2) )

  integrand_xi=phr(r)*(1.d0+xi(r))*f_v12 * Hubble/(1.d0+redshift)

  return
end function integrand_xi

function phr(r)
  use param
  implicit none
  real*8 :: phr,r
  if (r_eq .eq. 0.d0) then
     phr=1.d0
  else
     phr=((r/r_eq)**(-2)+1.d0)**(-betaN+1.)
  endif
  return
end function phr
  
! 2-point galaxy-absorber correlation function 
function xi(r) 
  use param
  implicit none
  real*8 :: xi,r
  xi=(r/r0)**(-slope_xi) 
  return
end function xi

! pairwise velocity field
function v12(r)
  use param
  implicit none
  real*8 :: v12,r
  real*8 :: v12_inflow, v12_outflow, cover_fraction

  ! only with inflow -----------
  !v12=-v_inflow*(r/r_inflow)/(1.0+(r/r_inflow)**slope_v12) ! self-similar clustering
  !-----------------------------

  ! with outflow ---------------
  v12_inflow=-v_inflow*(r/r_inflow)/(1.0+(r/r_inflow)**slope_v12) ! self-similar clustering

  ! if (v_outflow > 0.0) then
  !    v_esc=155.0 ! km/s escape velocity from isothermal halo of 10^11 Msun at z=3
  !    if (1.d0-(v_esc/v_outflow)**2*log(r/r_outflow) >= 0.0001) then
  !       v12_outflow=v_outflow*dsqrt( 1.d0-(v_esc/v_outflow)**2*log(r/r_outflow) )
  !    else
  !       v12_outflow=0.d0
  !    end if
  ! end if
  v12_outflow=v_outflow
  cover_fraction=exp(-r/r_outflow) ! covering fraction of outflowing gas

  !v12=v12_inflow+v12_outflow
  v12=(1.d0-cover_fraction)*v12_inflow+cover_fraction*v12_outflow
  !-----------------------------

  return
end function v12

! pairwise velocity dispersion
function s12(r)
  use param
  implicit none
  real*8 :: s12,r
  s12=s12_0
  return
end function s12


function optdpt(x,NHI)
  use param
  implicit none
  real*8 :: optdpt,NHI
  real*8 :: x,voigt
  optdpt=cross_sec*NHI*voigt(a,x)/doppler_width
  return
end function optdpt


! column density distribution function from Haardt&Madau 2014 [cm2]
function dNdNHIdz(NHI)
  use param, only : redshift
  implicit none
  real*8 :: dNdNHIdz,NHI
  real*8 :: A,beta,gamma,C
  real*8 :: A1,beta1,gamma1,C1
  real*8 :: A2,beta2,gamma2,C2
  real*8 :: logy1,logy2

  if ( (log10(NHI) >= 11.0) .and. (log10(NHI) < 15.0) ) then
     A=1.2e7
     beta=1.5
     gamma=3.0
     dNdNHIdz=A*NHI**(-beta)*(1.0+redshift)**gamma
  else if ( (log10(NHI) >= 15.0) .and. (log10(NHI) < 17.5) ) then
     A=3.8e14
     beta=2.0
     gamma=3.0
     dNdNHIdz=A*NHI**(-beta)*(1.0+redshift)**gamma
  else if ( (log10(NHI) >= 19.0) .and. (log10(NHI) < 20.3) ) then
     A=0.45
     beta=1.05
     gamma=1.27
     dNdNHIdz=A*NHI**(-beta)*(1.0+redshift)**gamma
  else if ( (log10(NHI) >= 20.3) .and. (log10(NHI) < 21.55) ) then
     A=8.7e18
     beta=2.0
     gamma=1.27
     dNdNHIdz=A*NHI**(-beta)*(1.0+redshift)**gamma
  else if ( (log10(NHI) >= 17.5) .and. (log10(NHI) < 19.0) ) then
     A1=3.8e14
     beta1=2.0
     gamma1=3.0
     C1=A1*(1.0+redshift)**gamma1

     A2=0.45
     beta2=1.05
     gamma2=1.27
     C2=A2*(1.0+redshift)**gamma2
            
     logy1=log10(C1)-17.5*beta1
     logy2=log10(C2)-19.0*beta2
        
     beta=(logy1-logy2)/(19.0-17.5)
     C=(10.**logy1)*(10.**17.5)**beta

     dNdNHIdz=C*NHI**(-beta)
  else
     dNdNHIdz=0.0
  end if

  return
end function dNdNHIdz


! Voigt-Hjerting function for Voigt profile
! Analytic form from Argyro Tasitsiomi (2006)
! dimensionless Voigt profile
function voigt(a,x)
  implicit none
  real*8 :: voigt, a, x
  ! internal variables used for analytic fit
  real*8 :: xx, y, q, pi
  pi=3.14159265359d0
  xx=x*x
  y=(xx-0.855d0)/(xx+3.42d0)
  if (y .le. 0.d0) then
     q=0.d0
  else if (y .gt. 0.d0) then
     q=y*(1.d0+21.d0/xx)* a/(pi*(xx+1.d0)) * ( 0.1117d0+y*( 4.421d0+y*(-9.207d0+5.674d0*y) ) )
  else
     print*, 'error in Voigt-Hjerting function. stop.'
     stop
  endif

  if (xx .lt. 500.d0) then  ! this prevents the underflow 
     voigt=q+exp(-xx)/1.77245385d0
  else
     voigt=q
  endif

  return
end function voigt

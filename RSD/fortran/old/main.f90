program main
  use param
  implicit none
  integer, parameter :: N=100
  real*8, parameter :: s_min=-8.d0, s_max=+8.d0 ! cMpc/h
  integer :: i,j
  real*8 :: xi_s, s_perp, s_para, ds
  real*8 :: xi, r, dr
  real*8 :: v12,s12
  character(len=128) :: infile, outfile

  call read_arg
  ! read input parameters
  open(1,file=infile)
  read(1,*) r0,slope_xi
  read(1,*) v_inflow,r_inflow,slope_v12
  read(1,*) v_outflow,r_outflow
  read(1,*) s12_0
  read(1,*) r_eq, betaN 
  close(1)

  write(*,*) 'input parameter from : ',infile
  write(*,'(A,2F8.3)') 'r0, slope_xi: ',r0,slope_xi
  write(*,'(A,3F8.3)') 'v_inflow, r_inflow, slope_v12: ',v_inflow,r_inflow,slope_v12
  write(*,'(A,2F8.3)') 'v_outflow, r_outflow: ',v_outflow,r_outflow
  write(*,'(A,F8.3)') 's12_0: ', s12_0
  write(*,'(A,2F8.3)') 'r_eq, betaN: ',r_eq, betaN 

  ! write the redshift-space correlation function
  print*,'output to : ', outfile
  open(10,file=outfile)
  ds=(s_max-s_min)/N
  do i=1,N
     do j=1,N

        s_para=s_min+(i-0.5)*ds
        s_perp=s_min+(j-0.5)*ds
        !if (abs(s_perp) <= 1.e-3) then
        !   print*, 'division by almost zero will happen'
        !end if
        xi= xi_s(s_perp,s_para)
        write(10,*) i,j,s_para,s_perp, xi
        write(*,'(2I5,3ES15.4)') i,j,s_para, s_perp,xi


     enddo
     write(10,*) ''
  enddo
  close(10)

  
contains                                           
  subroutine read_arg
    implicit none
    integer            :: i,n 
    integer            :: iargc                          
    character(len=4)   :: opt                       
    character(len=128) :: arg                       

    n = iargc()
    if (n < 3) then
       print *, 'usage: RSD  -in  param.in'
       print *, '            -out RSD.dat'
       print *, 'e.g: skewers -in param.in -out RSD.dat'
       print *, ' '
       stop
    end if
    
    do i = 1,n,2
       call getarg(i,opt)
       if (i == n) then
          print '("option ",a2," has no argument")', opt
          stop 2
       end if
       call getarg(i+1,arg)
       select case (opt)
       case('-in')
          infile=trim(arg)
       case ('-out')
          outfile=trim(arg)
       case default
          print '("unknown option ",a2," ignored")', opt
       end select
    end do
   
    return    
  end subroutine read_arg


end program main



function xi_s(s_perp,s_para)
  implicit none
  real*8 :: xi_s,s_perp,s_para
  ! quadpack parameters
  real*8, parameter :: epsabs=0.0d0
  real*8, parameter :: epsrel=1.d-6
  real*8 :: result,abserr,y_min,y_max
  integer :: neval,ier
  real*8 :: s1,s2
  common /svar/ s1,s2
  external intg

  s1=s_perp
  s2=s_para
  y_min=-1000.d0
  y_max=+1000.d0
  call qags(intg, y_min, y_max, epsabs, epsrel, result, abserr, neval, ier)
  xi_s=result-1.d0

  return
end function xi_s


! function xi_s(s_perp,s_para)
!   implicit none
!   real*8 :: xi_s,s_perp,s_para
!   ! quadpack parameters
!   real*8, parameter :: eps=1.0e-6
!   real*8 :: result,y_min,y_max
!   integer :: nint
!   real*8 :: s1,s2
!   common /svar/ s1,s2
!   external intg

!   s1=s_perp
!   s2=s_para
!   y_min=-500.0
!   y_max=+500.0
!   call gauss2(intg,y_min,y_max,eps,result,nint)
!   xi_s=result-1.d0

!   return
! end function xi_s

function intg(y)
  use param
  implicit none
  real*8 :: intg,y
  real*8 :: r,mu,f_v12
  real*8 :: xi,v12,s12,phr
  real*8 :: s1,s2
  common /svar/ s1,s2

  r=sqrt(s1**2+y**2)
  mu=y/r
  ! gaussian streaming model
  f_v12=1.d0/(dsqrt(2.d0*pi)*s12(r)/Hubble) * &
        exp( -(s2-y-mu*v12(r)/Hubble)**2/(2.d0*(s12(r)/Hubble)**2) )      

  intg=phr(r)*(1.d0+xi(r))*f_v12

  return
end function intg


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
  v12=-v_inflow/(1.0+(r/r_inflow)**slope_v12)
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

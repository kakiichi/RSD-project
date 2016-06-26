program main
  use param
  implicit none
  integer, parameter :: Ng=100
  !real*8, parameter :: vel_min=-900.d0, vel_max=+900.d0     ! km/s
  !real*8, parameter :: r_perp_min=-10.d0, r_perp_max=+10.d0 ! cMpc/h
  real*8, parameter :: vel_min=0.d0, vel_max=+900.d0     ! km/s
  real*8, parameter :: r_perp_min=0.0d0, r_perp_max=+10.d0 ! cMpc/h
  integer :: i,j
  real*8 :: s_para, s_perp
  real*8,allocatable,dimension(:) :: vel, r_perp
  real*8,allocatable,dimension(:,:) :: optdpt_eff
  real*8,parameter :: dvel=(vel_max-vel_min)/Ng
  real*8,parameter :: dr_perp=(r_perp_max-r_perp_min)/Ng
  real*8 :: tau_eff
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

  ! compute the RSD of galaxy-lya flux clustering
  allocate(optdpt_eff(Ng,Ng))
  allocate(vel(Ng))
  allocate(r_perp(Ng))

  do i=1,Ng
     do j=1,Ng
        vel(i)=vel_min+(i-0.5)*dvel
        r_perp(j)=r_perp_min+(j-0.5)*dr_perp

        optdpt_eff(i,j)=tau_eff(vel(i),r_perp(j))
        write(*,'(2I5,3ES15.4)') i,j,vel(i),r_perp(j), optdpt_eff(i,j)

     end do
  end do

  ! write the redshift-space correlation function
  print*,'output to : ', outfile
  open(10,file=outfile)
  write(10,*) '# 1) i 2) j 3) velocity shift [km/s] &
              4) comoving perpendicular distance r_perp(j) [cMpc/h] &
              5) proper z-space distance line-of-sight, s_para [pMpc] 6) s_perp, optdpt_eff(i,j)'
  do i=1,Ng
     do j=1,Ng
        s_para=vel(i)/Hubble/h0            ! pMpc
        s_perp=r_perp(j)/(1.+redshift)/h0  ! pMpc
        write(10,'(2I6, 4F12.4, ES15.4)') i,j,vel(i),r_perp(j), s_para, s_perp, optdpt_eff(i,j)

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



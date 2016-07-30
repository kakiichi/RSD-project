program main
  use param
  implicit none
  integer, parameter :: Ng=100
  real*8, parameter :: vel_min=-900.d0, vel_max=+900.d0     ! km/s
  integer :: i
  real*8 :: s_para
  real*8,allocatable,dimension(:) :: vel
  real*8,allocatable,dimension(:) :: optdpt_eff
  real*8,parameter :: dvel=(vel_max-vel_min)/Ng
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
  allocate(optdpt_eff(Ng))
  allocate(vel(Ng))

  do i=1,Ng
     vel(i)=vel_min+(i-0.5)*dvel
     optdpt_eff(i)=tau_eff(vel(i))
     write(*,'(1I5,2ES15.4)') i, vel(i), optdpt_eff(i)
  end do

  ! write the redshift-space correlation function
  print*,'output to : ', outfile
  open(10,file=outfile)
  write(10,*) '# 1) i 2) velocity shift [km/s] 3) proper line-of-sight distance, s_para [pMpc] 4) optdpt_eff(i)'
  do i=1,Ng
     s_para=vel(i)/Hubble/h0            ! pMpc
     write(10,'(1I6, 2F12.4, ES15.4)') i,vel(i), s_para, optdpt_eff(i)
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



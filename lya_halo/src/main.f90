! lyman alpha halo calculation
program main
  use param
  implicit none
  integer, parameter :: Ng=500
  real*8, parameter :: r_min=0.0d0, r_max=+10.d0 ! cMpc/h
  real*8,parameter :: dr=(r_max-r_min)/Ng
  integer :: i
  real*8,allocatable,dimension(:) :: epsilon_lya
  real*8,allocatable,dimension(:) :: r
  real*8 :: emissivity
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

  ! compute the lyman-alpha emissivity
  allocate(epsilon_lya(Ng))
  allocate(r(Ng))

  do i=1,Ng
     r(i)=r_min+(i-0.5)*dr
     epsilon_lya(i)=emissivity(r(i))

     write(*,'(I5,3ES15.4)') i,r(i),r(i)/(1.d0+redshift)/h0,epsilon_lya(i)
  end do


  ! write the redshift-space correlation function
  print*,'output to : ', outfile
  open(10,file=outfile)
  write(10,*) '# 1) i 2) r [pMpc] 3) Lya emissivity [erg/s/cm3, proper]'
  do i=1,Ng
     write(10,'(I5,2ES15.4)') i,r(i)/(1.d0+redshift)/h0,epsilon_lya(i)
  enddo
  close(10)

  deallocate(epsilon_lya)
  deallocate(r)
  
contains                                           
  subroutine read_arg
    implicit none
    integer            :: i,n 
    integer            :: iargc                          
    character(len=4)   :: opt                       
    character(len=128) :: arg                       

    n = iargc()
    if (n < 3) then
       print *, 'usage: lya_halo  -in  param.in'
       print *, '            -out lya_halo.dat'
       print *, 'e.g: lya_halo -in param.in -out lya_halo.dat'
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

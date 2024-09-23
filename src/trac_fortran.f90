PROGRAM trac_fortran
  USE iso_c_binding
  USE mptrac_struct
  USE mptrac_func
  USE, intrinsic :: iso_fortran_env
  USE, intrinsic :: iso_c_binding
  !! ToDo clean up and check for consitency
  !!What is the difference between  "USE iso_c_binding" and "USE, intrinsic :: iso_c_binding"?
  IMPLICIT NONE

  CHARACTER(len=40) :: filename_ctl, filename_atm, dirname
  INTEGER(c_int) :: argc
  CHARACTER(len=15) :: arg1, arg2, arg3, arg4, arg5
  CHARACTER(c_char), TARGET :: temp1, temp2, temp3, temp4, temp5
  TYPE(c_ptr), DIMENSION(5), TARGET :: argv_ptrs
  TYPE(ctl_t) :: ctl
  TYPE(atm_t) :: atm
  TYPE(clim_t) :: clim
  TYPE(met_t), TARGET :: met0, met1
  TYPE(met_t), POINTER :: met0p, met1p
  REAL(real64) :: t
  REAL(real64), DIMENSION(npp) :: dt
  
  filename_ctl = "data/trac.ctl"
  filename_atm = "atm_split.tab"

  !! ToDo: Create command line arguments
  argc = 5
  arg1 = "FortranWrapper" 
  arg2 = "abc2"
  arg3 = "abc3"
  arg4 = "abc4"
  arg5 = "abc5"

  ! Create 5 different pointers and store them in array
  temp1 = trim(arg1)//c_null_char
  argv_ptrs(1) = c_loc(temp1)
  temp2 = trim(arg2)//c_null_char
  argv_ptrs(2) = c_loc(temp2)
  temp3 = trim(arg3)//c_null_char
  argv_ptrs(3) = c_loc(temp3)
  temp4 = trim(arg4)//c_null_char
  argv_ptrs(4) = c_loc(temp4)
  temp5 = trim(arg5)//c_null_char
  argv_ptrs(5) = c_loc(temp5)

  print *,'Hello World!'
  
  ! ! Read control parameters... 
  CALL mptrac_read_ctl(TRIM(filename_ctl)//char(0), argc, argv_ptrs, ctl)
  print *,'Hello Control'
  
  ! ! Read climatological data... 
  CALL mptrac_read_clim(ctl, clim)
  print *,'Hello Clim'
  
  ! Read atmospheric data...
  CALL mptrac_read_atm(TRIM(filename_atm)//char(0), ctl, atm)
  print *, "Hello Atm"
  
  ! Initialize timesteps...
  CALL mptrac_module_timesteps_init(ctl, atm)
  print *, "Hello Timesteps"

  print *,c_sizeof(met1)
  print *,c_sizeof(atm)
  print *,c_sizeof(ctl)
  print *,c_sizeof(clim)
  
  
  WRITE(*,*) ctl%t_start
  ! Get meteo data...
  met0p => met0
  met1p => met1
  CALL mptrac_get_met(ctl, clim, ctl%t_start, met0p, met1p)
  print *, "Hello Met"
  
  dirname = "data"

  t = ctl%t_start

  DO WHILE (ctl%direction * (t - ctl%t_stop) < ctl%dt_mod)

     ! Adjust length of final time step... 
     IF (ctl%direction * (t - ctl%t_stop) > 0) THEN
        t = ctl%t_stop
     ENDIF

     ! Set time steps of air parcels...
     CALL mptrac_module_timesteps(ctl, met0, atm, dt, t)

     IF (t .NE. ctl%t_start) THEN
        ! Get meteo data...
        CALL mptrac_get_met(ctl, clim, t, met0p, met1p)
     ENDIF

     ! Advection...
     CALL mptrac_module_advect(ctl, met0, met1, atm, dt)

     ! Write output...
     CALL mptrac_write_output(TRIM(dirname)//char(0), ctl, met0, met1, atm, t)

     t = t + ctl%direction * ctl%dt_mod

  END DO

END PROGRAM trac_fortran

  

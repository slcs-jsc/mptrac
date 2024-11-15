PROGRAM trac_fortran
  USE mptrac_struct
  USE mptrac_func
  USE iso_fortran_env
  USE iso_c_binding

  IMPLICIT NONE

  CHARACTER(len=40) :: filename_ctl, filename_atm, dirname
  INTEGER(c_int) :: argc
  TYPE(ctl_t) :: ctl
  TYPE(atm_t) :: atm
  TYPE(clim_t) :: clim
  TYPE(met_t), TARGET :: met0, met1
  TYPE(met_t), POINTER :: met0p, met1p
  REAL(real64) :: t
  REAL(real64), DIMENSION(npp) :: dt
  CHARACTER(len=32) :: arg
  TYPE(c_ptr), ALLOCATABLE, DIMENSION(:) :: argv_ptrs
  CHARACTER(len=32), ALLOCATABLE, DIMENSION(:), TARGET :: tmp
  INTEGER :: i, stat
  
  ! Read command line arguments...
  argc = command_argument_count()

  IF (argc < 4) THEN
     WRITE(*,*) "Give parameters: <dirlist> <ctl> <atm_in>"
     CALL ABORT
  ENDIF
  
  ALLOCATE(tmp(argc), argv_ptrs(argc))
  
  DO i = 1, argc
     CALL get_command_argument(i, arg)
     IF (LEN_TRIM(arg) == 0) EXIT
     tmp(i) = TRIM(arg)//c_null_char
     argv_ptrs(i) = c_loc(tmp(i))
  ENDDO
  
  ! Open directory list...
  OPEN(10,file=tmp(1),iostat=stat)
  IF (stat .ne. 0) THEN
     WRITE(*,*) "Cannot open directory list!"
     CALL ABORT
  ENDIF
  
  DO WHILE (1 .eq. 1)
     READ(10,'(a)', END=200) dirname
 
     filename_ctl = TRIM(dirname)//"/"//tmp(2)
     filename_atm = TRIM(dirname)//"/"//tmp(3)

     ! Read control parameters...
     CALL mptrac_read_ctl(TRIM(filename_ctl)//c_null_char, argc, argv_ptrs, ctl)
  
     ! Read climatological data... 
     CALL mptrac_read_clim(ctl, clim)
  
     ! Read atmospheric data...
     CALL mptrac_read_atm(TRIM(filename_atm)//c_null_char, ctl, atm)
  
     ! Initialize timesteps...
     CALL mptrac_module_timesteps_init(ctl, atm)

     ! Get meteo data...
     met0p => met0
     met1p => met1
     CALL mptrac_get_met(ctl, clim, ctl%t_start, met0p, met1p)

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
        CALL mptrac_write_output(TRIM(dirname)//c_null_char, ctl, met0, met1, atm, t)

        t = t + ctl%direction * ctl%dt_mod

     END DO

  END DO
200  CONTINUE

END PROGRAM trac_fortran

  

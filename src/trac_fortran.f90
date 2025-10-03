! This file is part of MPTRAC.
! 
! MPTRAC is free software: you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! MPTRAC is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with MPTRAC. If not, see <http://www.gnu.org/licenses/>.
! 
! Copyright (C) 2013-2025 Forschungszentrum Juelich GmbH

PROGRAM trac_fortran
  
  USE mptrac_struct
  USE mptrac_func
  USE iso_fortran_env
  USE iso_c_binding
  
  IMPLICIT NONE
  
  CHARACTER(len=40) :: filename_ctl, filename_atm, dirname
  INTEGER(c_int) :: argc, ntask = -1
  TYPE(ctl_t), TARGET :: ctl
  TYPE(ctl_t), POINTER :: ctlp
  TYPE(cache_t), TARGET :: cache
  TYPE(cache_t), POINTER :: cachep
  TYPE(atm_t), TARGET :: atm
  TYPE(atm_t), POINTER :: atmp
  TYPE(clim_t), TARGET :: clim
  TYPE(clim_t), POINTER :: climp
  TYPE(met_t), TARGET :: met0, met1
  TYPE(met_t), POINTER :: met0p, met1p
  TYPE(dd_t), TARGET :: dd
  TYPE(dd_t), POINTER :: ddp
  REAL(real64) :: t
  CHARACTER(len=5000) :: arg
  TYPE(c_ptr), ALLOCATABLE, DIMENSION(:) :: argv_ptrs
  CHARACTER(len=5000), ALLOCATABLE, DIMENSION(:), TARGET :: tmp
  INTEGER :: i, stat
  
  ! Read command line arguments...
  argc = command_argument_count()
  
  IF (argc < 4) THEN
     WRITE(*,*) "Error: Give parameters: <dirlist> <ctl> <atm_in>"
     CALL EXIT
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
     WRITE(*,*) "Error: Cannot open directory list!"
     CALL EXIT
  ENDIF
  
  ! Set pointers...
  ctlp => ctl
  cachep => cache
  climp => clim
  met0p => met0
  met1p => met1
  atmp => atm
  ddp => dd
  
  ! Endless loop...
  DO WHILE (1 .eq. 1)
     READ(10,'(a)', END=200) dirname
     
     filename_ctl = TRIM(dirname)//"/"//tmp(2)
     filename_atm = TRIM(dirname)//"/"//tmp(3)
     
     ! Allocate memory...
     CALL mptrac_alloc(ctlp, cachep, climp, met0p, met1p, atmp, ddp)
     
     ! Read control parameters...
     CALL mptrac_read_ctl(TRIM(filename_ctl)//c_null_char, argc, argv_ptrs, ctl)
     
     ! Read climatological data... 
     CALL mptrac_read_clim(ctl, clim)
     
     ! Read atmospheric data...
     CALL mptrac_read_atm(TRIM(filename_atm)//c_null_char, ctl, atm)
     
     ! Initialize MPTRAC...
     CALL mptrac_init(ctl, cache, clim, atm, ntask)

     ! Set start time...
     t = ctl%t_start

     ! Loop over time steps...
     DO WHILE (ctl%direction * (t - ctl%t_stop) < ctl%dt_mod)
        
        ! Adjust length of final time step... 
        IF (ctl%direction * (t - ctl%t_stop) > 0) THEN
           t = ctl%t_stop
        ENDIF
        
        ! Get meteo data...
        CALL mptrac_get_met(ctl, clim, t, met0p, met1p, ddp)
        
        ! Check time step...
        IF (ctl%dt_mod > ABS(met0p%lon(2) - met0p%lon(1)) * 111132. / 150.) THEN
           WRITE(*,*) "Violation of CFL criterion! Check DT_MOD!"
        ENDIF
        
        ! Run a single time step..."
        CALL mptrac_run_timestep(ctl, cache, clim, met0p, met1p, atm, t, ddp)
        
        ! Write output...
        CALL mptrac_write_output(TRIM(dirname)//c_null_char, ctl, met0p, met1p, atm, t)

        ! Set time...
        t = t + ctl%direction * ctl%dt_mod
        
     END DO
     
     ! Free memory...
     CALL mptrac_free(ctlp, cachep, climp, met0p, met1p, atmp, ddp)
     
  END DO
200 CONTINUE
  
END PROGRAM trac_fortran

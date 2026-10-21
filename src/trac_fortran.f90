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
! Copyright (C) 2013-2026 Forschungszentrum Juelich GmbH

PROGRAM trac_fortran
  
  USE mptrac_struct
  USE mptrac_func
  USE iso_fortran_env
  USE iso_c_binding
#ifdef MPI
  USE mpi
#endif
  
  IMPLICIT NONE
  
  CHARACTER(len=length) :: filename_ctl, filename_atm, dirname
  INTEGER(c_int) :: argc, ntask = -1, read_atm_status
  TYPE(ctl_t), POINTER :: ctlp
  TYPE(cache_t), POINTER :: cachep
  TYPE(atm_t), POINTER :: atmp
  TYPE(depo_t), POINTER :: depop
  TYPE(clim_t), POINTER :: climp
  TYPE(met_t), POINTER :: met0p, met1p
  TYPE(dd_t), POINTER :: ddp
  REAL(real64) :: t
  CHARACTER(len=length) :: arg
  TYPE(c_ptr), ALLOCATABLE, DIMENSION(:) :: argv_ptrs
  CHARACTER(len=length), ALLOCATABLE, DIMENSION(:), TARGET :: tmp
  INTEGER :: i, stat, rank = 0, size = 1
#ifdef MPI
  INTEGER :: ierr
#endif

#ifdef MPI
  CALL MPI_Init(ierr)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  CALL MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
#endif
  
  ! Read command line arguments...
  argc = command_argument_count()
  
  IF (argc < 3) THEN
     WRITE(*,*) "Error: Give parameters: <dirlist> <ctl> <atm_in>"
     CALL EXIT
  ENDIF
  
  ALLOCATE(tmp(0:argc), argv_ptrs(0:argc))
  
  DO i = 0, argc
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

  NULLIFY(ctlp, cachep, climp, met0p, met1p, atmp, depop, ddp)
  
  ! Endless loop...
  DO WHILE (1 .eq. 1)
     READ(10,'(a)', END=200) dirname

     ntask = ntask + 1
     IF (MOD(ntask, size) /= rank) CYCLE
     
     filename_ctl = TRIM(dirname)//"/"//tmp(2)
     filename_atm = TRIM(dirname)//"/"//tmp(3)
     
     ! Allocate memory...
     CALL mptrac_alloc(ctlp, cachep, climp, met0p, met1p, atmp, depop, ddp)
     
     ! Read control parameters...
     CALL mptrac_read_ctl(TRIM(filename_ctl)//c_null_char, argc + 1, argv_ptrs, ctlp)
     ctlp%dd_supported = 1
     
     ! Read climatological data... 
     CALL mptrac_read_clim(ctlp, climp)
     
     ! Read atmospheric data...
     read_atm_status = mptrac_read_atm(TRIM(filename_atm)//c_null_char, ctlp, atmp)
     IF (read_atm_status == 0) THEN
        WRITE(*,*) "Error: Cannot open file!"
        CALL EXIT
     ENDIF
     
     ! Initialize MPTRAC...
     CALL mptrac_init(ctlp, cachep, climp, atmp, depop, ntask)

     ! Set start time...
     t = ctlp%t_start

     ! Loop over time steps...
     DO WHILE (ctlp%direction * (t - ctlp%t_stop) < ctlp%dt_mod)
        
        ! Adjust length of final time step... 
        IF (ctlp%direction * (t - ctlp%t_stop) > 0) THEN
           t = ctlp%t_stop
        ENDIF
        
        ! Get meteo data...
        CALL mptrac_get_met(ctlp, climp, t, met0p, met1p, ddp)
        
        ! Check time step...
        IF (ctlp%dt_mod > ABS(met0p%lon(2) - met0p%lon(1)) * 111132. / 150.) THEN
           WRITE(*,*) "Violation of CFL criterion! Check DT_MOD!"
        ENDIF

#ifdef DD
        IF (t == ctlp%t_start .AND. ctlp%dd /= 0) CALL dd_init(ctlp, ddp, atmp)
#else
        IF (ctlp%dd /= 0 .OR. ctlp%dd_subdomains_zonal * ctlp%dd_subdomains_meridional > 1) THEN
           WRITE(*,*) "Model configured for domain decomposition, but not compiled for it!"
           STOP 1
        ENDIF
#endif
        
        ! Run a single time step..."
        CALL mptrac_run_timestep(ctlp, cachep, climp, met0p, met1p, atmp, depop, t, ddp)
        
        ! Write output...
        CALL mptrac_write_output(TRIM(dirname)//c_null_char, ctlp, met0p, met1p, atmp, depop, t)

        ! Set time...
        t = t + ctlp%direction * ctlp%dt_mod
        
     END DO
     
     ! Free memory...
     CALL mptrac_free(ctlp, cachep, climp, met0p, met1p, atmp, depop, ddp)
     
  END DO
200 CONTINUE
  CLOSE(10)
#ifdef MPI
  CALL MPI_Finalize(ierr)
#endif
  
END PROGRAM trac_fortran

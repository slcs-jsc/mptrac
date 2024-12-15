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
! Copyright (C) 2013-2024 Forschungszentrum Juelich GmbH

! MPTRAC data sturctures
MODULE mptrac_struct
  
  USE iso_c_binding
  IMPLICIT NONE

  !! Values must be identical to equivalent variables in mptrac.h!
  !! The values for ex, ey, ep are suited for ERA5 data. However, they
  !! exceed some size limit of NVHPC (2GB???), where the compiler fails
  !! to build the fortran wrapper. GCC and ICX work fine.
  INTEGER, PARAMETER :: ex = 1202
  INTEGER, PARAMETER :: ey = 602
  INTEGER, PARAMETER :: ep = 140
  !! Alternative smaller values, good enough for ERA-interim.
  ! INTEGER, PARAMETER :: ex = 481
  ! INTEGER, PARAMETER :: ey = 241
  ! INTEGER, PARAMETER :: ep = 60
  INTEGER, PARAMETER :: npp = 10000000 !NP
  INTEGER, PARAMETER :: nqq = 15       !NQ
  INTEGER, PARAMETER :: length = 5000  !LEN
  INTEGER, PARAMETER :: cyy = 250      !CY
  INTEGER, PARAMETER :: co3 = 30
  INTEGER, PARAMETER :: cp = 70
  INTEGER, PARAMETER :: csza = 50
  INTEGER, PARAMETER :: ct = 12
  INTEGER, PARAMETER :: cts = 1000
  
  !! The order of the variables in each struct matters!
  TYPE, bind(c) :: atm_t
     INTEGER(c_int) :: np
     REAL(c_double), DIMENSION(npp) :: time
     REAL(c_double), DIMENSION(npp) :: p
     REAL(c_double), DIMENSION(npp) :: lon
     REAL(c_double), DIMENSION(npp) :: lat
     REAL(c_double), DIMENSION(npp,nqq) :: q
  END TYPE atm_t

  !! in 3D arrays first and third index swapped
  TYPE, bind(c) :: clim_photo_t
     INTEGER(c_int) :: np
     INTEGER(c_int) :: nsza
     INTEGER(c_int) :: no3c
     REAL(c_double), DIMENSION(cp) :: p
     REAL(c_double), DIMENSION(csza) :: sza
     REAL(c_double), DIMENSION(co3) :: o3c
     REAL(c_double), DIMENSION(co3,csza,cp) :: n2o
     REAL(c_double), DIMENSION(co3,csza,cp) :: ccl4
     REAL(c_double), DIMENSION(co3,csza,cp) :: ccl3f
     REAL(c_double), DIMENSION(co3,csza,cp) :: ccl2f2
     REAL(c_double), DIMENSION(co3,csza,cp) :: o2
     REAL(c_double), DIMENSION(co3,csza,cp) :: o3_1
     REAL(c_double), DIMENSION(co3,csza,cp) :: o3_2
     REAL(c_double), DIMENSION(co3,csza,cp) :: h2o2
     REAL(c_double), DIMENSION(co3,csza,cp) :: h2o
  END TYPE clim_photo_t

  TYPE, bind(c) :: clim_ts_t
     INTEGER(c_int) :: ntime
     REAL(c_double), DIMENSION(cts) :: time
     REAL(c_double), DIMENSION(cts) :: vmr
  END TYPE clim_ts_t

  TYPE, bind(c) :: clim_zm_t
     INTEGER(c_int) :: ntime
     INTEGER(c_int) :: nlat
     INTEGER(c_int) :: np
     REAL(c_double), DIMENSION(ct) :: time
     REAL(c_double), DIMENSION(cyy) :: lat
     REAL(c_double), DIMENSION(cp) :: p
     REAL(c_double), DIMENSION(cyy,cp,ct) :: vmr
  END TYPE clim_zm_t

  TYPE, bind(c) :: clim_t
     INTEGER(c_int) :: tropo_ntime
     INTEGER(c_int) :: tropo_nlat
     REAL(c_double), DIMENSION(12) :: tropo_time
     REAL(c_double), DIMENSION(73) :: tropo_lat
     REAL(c_double), DIMENSION(73,12) :: tropo
     TYPE(clim_photo_t) :: photo
     TYPE(clim_zm_t) :: hno3
     TYPE(clim_zm_t) :: oh
     TYPE(clim_zm_t) :: h2o2
     TYPE(clim_zm_t) :: ho2
     TYPE(clim_zm_t) :: o1d
     TYPE(clim_ts_t) :: ccl4
     TYPE(clim_ts_t) :: ccl3f
     TYPE(clim_ts_t) :: ccl2f2
     TYPE(clim_ts_t) :: n2o
     TYPE(clim_ts_t) :: sf6
  END TYPE clim_t
  
  TYPE, bind(c) :: ctl_t
     INTEGER(c_int) :: advect_cpl_zeta_and_press_modules
     INTEGER(c_int) :: met_press_level_def
     INTEGER(c_int) :: advect_vert_coord
     INTEGER(c_int) :: met_vert_coord
     INTEGER(c_int) :: met_clams
     INTEGER(c_int) :: nq
     CHARACTER(c_char), DIMENSION(length,nqq) :: qnt_name
     CHARACTER(c_char), DIMENSION(length,nqq) :: qnt_longname
     CHARACTER(c_char), DIMENSION(length,nqq) :: qnt_unit
     CHARACTER(c_char), DIMENSION(length,nqq) :: qnt_format
     INTEGER(c_int) :: qnt_idx
     INTEGER(c_int) :: qnt_ens
     INTEGER(c_int) :: qnt_stat
     INTEGER(c_int) :: qnt_m
     INTEGER(c_int) :: qnt_vmr
     INTEGER(c_int) :: qnt_rp
     INTEGER(c_int) :: qnt_rhop
     INTEGER(c_int) :: qnt_ps
     INTEGER(c_int) :: qnt_ts
     INTEGER(c_int) :: qnt_zs
     INTEGER(c_int) :: qnt_us
     INTEGER(c_int) :: qnt_vs
     INTEGER(c_int) :: qnt_lsm
     INTEGER(c_int) :: qnt_sst
     INTEGER(c_int) :: qnt_pbl
     INTEGER(c_int) :: qnt_pt
     INTEGER(c_int) :: qnt_tt
     INTEGER(c_int) :: qnt_zt
     INTEGER(c_int) :: qnt_h2ot
     INTEGER(c_int) :: qnt_zg
     INTEGER(c_int) :: qnt_p
     INTEGER(c_int) :: qnt_t
     INTEGER(c_int) :: qnt_rho
     INTEGER(c_int) :: qnt_u
     INTEGER(c_int) :: qnt_v
     INTEGER(c_int) :: qnt_w
     INTEGER(c_int) :: qnt_h2o
     INTEGER(c_int) :: qnt_o3
     INTEGER(c_int) :: qnt_lwc
     INTEGER(c_int) :: qnt_rwc
     INTEGER(c_int) :: qnt_iwc
     INTEGER(c_int) :: qnt_swc
     INTEGER(c_int) :: qnt_cc
     INTEGER(c_int) :: qnt_pct
     INTEGER(c_int) :: qnt_pcb
     INTEGER(c_int) :: qnt_cl
     INTEGER(c_int) :: qnt_plcl
     INTEGER(c_int) :: qnt_plfc
     INTEGER(c_int) :: qnt_pel
     INTEGER(c_int) :: qnt_cape
     INTEGER(c_int) :: qnt_cin
     INTEGER(c_int) :: qnt_o3c
     INTEGER(c_int) :: qnt_hno3
     INTEGER(c_int) :: qnt_oh
     INTEGER(c_int) :: qnt_h2o2
     INTEGER(c_int) :: qnt_ho2
     INTEGER(c_int) :: qnt_o1d
     INTEGER(c_int) :: qnt_mloss_oh
     INTEGER(c_int) :: qnt_mloss_h2o2
     INTEGER(c_int) :: qnt_mloss_kpp
     INTEGER(c_int) :: qnt_mloss_wet
     INTEGER(c_int) :: qnt_mloss_dry
     INTEGER(c_int) :: qnt_mloss_decay
     INTEGER(c_int) :: qnt_loss_rate
     INTEGER(c_int) :: qnt_psat
     INTEGER(c_int) :: qnt_psice
     INTEGER(c_int) :: qnt_pw
     INTEGER(c_int) :: qnt_sh
     INTEGER(c_int) :: qnt_rh
     INTEGER(c_int) :: qnt_rhice
     INTEGER(c_int) :: qnt_theta
     INTEGER(c_int) :: qnt_zeta
     INTEGER(c_int) :: qnt_zeta_d
     INTEGER(c_int) :: qnt_tvirt
     INTEGER(c_int) :: qnt_lapse
     INTEGER(c_int) :: qnt_vh
     INTEGER(c_int) :: qnt_vz
     INTEGER(c_int) :: qnt_pv
     INTEGER(c_int) :: qnt_tdew
     INTEGER(c_int) :: qnt_tice
     INTEGER(c_int) :: qnt_tsts
     INTEGER(c_int) :: qnt_tnat
     INTEGER(c_int) :: qnt_Cx
     INTEGER(c_int) :: qnt_Ch2o
     INTEGER(c_int) :: qnt_Co3
     INTEGER(c_int) :: qnt_Cco
     INTEGER(c_int) :: qnt_Coh
     INTEGER(c_int) :: qnt_Ch
     INTEGER(c_int) :: qnt_Cho2
     INTEGER(c_int) :: qnt_Ch2o2
     INTEGER(c_int) :: qnt_Co1d
     INTEGER(c_int) :: qnt_Co3p
     INTEGER(c_int) :: qnt_Cccl4
     INTEGER(c_int) :: qnt_Cccl3f
     INTEGER(c_int) :: qnt_Cccl2f2
     INTEGER(c_int) :: qnt_Cn2o
     INTEGER(c_int) :: qnt_Csf6
     INTEGER(c_int) :: qnt_aoa
     INTEGER(c_int) :: direction
     REAL(c_double) :: t_start
     REAL(c_double) :: t_stop
     REAL(c_double) :: dt_mod
     CHARACTER(c_char), DIMENSION(length) :: metbase
     REAL(c_double) :: dt_met
     INTEGER(c_int) :: met_convention
     INTEGER(c_int) :: met_type
     INTEGER(c_int) :: met_nc_scale
     INTEGER(c_int) :: met_nc_level
     INTEGER(c_int) :: met_nc_quant
     INTEGER(c_int) :: met_zfp_prec
     REAL(c_double) :: met_zfp_tol_t
     REAL(c_double) :: met_zfp_tol_z
     INTEGER(c_int) :: met_cms_batch
     INTEGER(c_int) :: met_cms_heur
     REAL(c_double) :: met_cms_eps_z
     REAL(c_double) :: met_cms_eps_t
     REAL(c_double) :: met_cms_eps_u
     REAL(c_double) :: met_cms_eps_v
     REAL(c_double) :: met_cms_eps_w
     REAL(c_double) :: met_cms_eps_pv
     REAL(c_double) :: met_cms_eps_h2o
     REAL(c_double) :: met_cms_eps_o3
     REAL(c_double) :: met_cms_eps_lwc
     REAL(c_double) :: met_cms_eps_rwc
     REAL(c_double) :: met_cms_eps_iwc
     REAL(c_double) :: met_cms_eps_swc
     REAL(c_double) :: met_cms_eps_cc
     INTEGER(c_int) :: met_dx
     INTEGER(c_int) :: met_dy
     INTEGER(c_int) :: met_dp
     INTEGER(c_int) :: met_sx
     INTEGER(c_int) :: met_sy
     INTEGER(c_int) :: met_sp
     REAL(c_double) :: met_detrend
     INTEGER(c_int) :: met_np
     REAL(c_double), DIMENSION(ep) :: met_p
     INTEGER(c_int) :: met_geopot_sx
     INTEGER(c_int) :: met_geopot_sy
     INTEGER(c_int) :: met_relhum
     INTEGER(c_int) :: met_cape
     INTEGER(c_int) :: met_pbl
     INTEGER(c_int) :: met_pbl_min
     INTEGER(c_int) :: met_pbl_max
     INTEGER(c_int) :: met_tropo
     REAL(c_double) :: met_tropo_pv
     REAL(c_double) :: met_tropo_theta
     INTEGER(c_int) :: met_tropo_spline
     REAL(c_double) :: met_dt_out
     INTEGER(c_int) :: met_cache
     INTEGER(c_int) :: met_mpi_share
     REAL(c_double) :: sort_dt
     INTEGER(c_int) :: isosurf
     CHARACTER(c_char), DIMENSION(length) :: balloon
     INTEGER(c_int) :: advect
     INTEGER(c_int) :: reflect
     INTEGER(c_int) :: rng_type
     REAL(c_double) :: turb_dx_trop
     REAL(c_double) :: turb_dx_strat
     REAL(c_double) :: turb_dz_trop
     REAL(c_double) :: turb_dz_strat
     REAL(c_double) :: turb_mesox
     REAL(c_double) :: turb_mesoz
     REAL(c_double) :: conv_cape
     REAL(c_double) :: conv_cin
     REAL(c_double) :: conv_dt
     INTEGER(c_int) :: conv_mix_bot
     INTEGER(c_int) :: conv_mix_top
     REAL(c_double) :: bound_mass
     REAL(c_double) :: bound_mass_trend
     REAL(c_double) :: bound_vmr
     REAL(c_double) :: bound_vmr_trend
     REAL(c_double) :: bound_lat0
     REAL(c_double) :: bound_lat1
     REAL(c_double) :: bound_p0
     REAL(c_double) :: bound_p1
     REAL(c_double) :: bound_dps
     REAL(c_double) :: bound_dzs
     REAL(c_double) :: bound_zetas
     INTEGER(c_int) :: bound_pbl
     CHARACTER(c_char), DIMENSION(length) :: species
     REAL(c_double) :: molmass
     REAL(c_double) :: tdec_trop
     REAL(c_double) :: tdec_strat
     CHARACTER(c_char), DIMENSION(length) :: clim_photo
     CHARACTER(c_char), DIMENSION(length) :: clim_hno3_filename
     CHARACTER(c_char), DIMENSION(length) :: clim_oh_filename
     CHARACTER(c_char), DIMENSION(length) :: clim_h2o2_filename
     CHARACTER(c_char), DIMENSION(length) :: clim_ho2_filename
     CHARACTER(c_char), DIMENSION(length) :: clim_o1d_filename
     CHARACTER(c_char), DIMENSION(length) :: clim_o3_filename
     CHARACTER(c_char), DIMENSION(length) :: clim_ccl4_timeseries
     CHARACTER(c_char), DIMENSION(length) :: clim_ccl3f_timeseries
     CHARACTER(c_char), DIMENSION(length) :: clim_ccl2f2_timeseries
     CHARACTER(c_char), DIMENSION(length) :: clim_n2o_timeseries
     CHARACTER(c_char), DIMENSION(length) :: clim_sf6_timeseries
     REAL(c_double) :: mixing_dt
     REAL(c_double) :: mixing_trop
     REAL(c_double) :: mixing_strat
     INTEGER(c_int) :: mixing_nz
     REAL(c_double) :: mixing_z0
     REAL(c_double) :: mixing_z1
     INTEGER(c_int) :: mixing_nx
     REAL(c_double) :: mixing_lon0
     REAL(c_double) :: mixing_lon1
     INTEGER(c_int) :: mixing_ny
     REAL(c_double) :: mixing_lat0
     REAL(c_double) :: mixing_lat1
     INTEGER(c_int) :: chemgrid_nz
     REAL(c_double) :: chemgrid_z0
     REAL(c_double) :: chemgrid_z1
     INTEGER(c_int) :: chemgrid_nx
     REAL(c_double) :: chemgrid_lon0
     REAL(c_double) :: chemgrid_lon1
     INTEGER(c_int) :: chemgrid_ny
     REAL(c_double) :: chemgrid_lat0
     REAL(c_double) :: chemgrid_lat1
     INTEGER(c_int) :: oh_chem_reaction
     REAL(c_double), DIMENSION(4) :: oh_chem
     REAL(c_double) :: oh_chem_beta
     INTEGER(c_int) :: h2o2_chem_reaction
     INTEGER(c_int) :: kpp_chem
     REAL(c_double) :: dt_kpp
     INTEGER(c_int) :: tracer_chem
     REAL(c_double), DIMENSION(2) :: wet_depo_pre
     REAL(c_double) :: wet_depo_bc_a
     REAL(c_double) :: wet_depo_bc_b
     REAL(c_double) :: wet_depo_ic_a
     REAL(c_double) :: wet_depo_ic_b
     REAL(c_double), DIMENSION(2) :: wet_depo_ic_h
     REAL(c_double), DIMENSION(2) :: wet_depo_bc_h
     REAL(c_double) :: wet_depo_so2_ph
     REAL(c_double) :: wet_depo_ic_ret_ratio
     REAL(c_double) :: wet_depo_bc_ret_ratio
     REAL(c_double) :: dry_depo_dp
     REAL(c_double) :: dry_depo_vdep
     REAL(c_double) :: psc_h2o
     REAL(c_double) :: psc_hno3
     CHARACTER(c_char), DIMENSION(length) :: atm_basename
     CHARACTER(c_char), DIMENSION(length) :: atm_gpfile
     REAL(c_double) :: atm_dt_out
     INTEGER(c_int) :: atm_filter
     INTEGER(c_int) :: atm_stride
     INTEGER(c_int) :: atm_type
     INTEGER(c_int) :: atm_type_out
     INTEGER(c_int) :: atm_nc_level
     INTEGER(c_int), DIMENSION(nqq) :: atm_nc_quant
     INTEGER(c_int) :: obs_type
     CHARACTER(c_char), DIMENSION(length) :: csi_basename
     CHARACTER(c_char), DIMENSION(length) :: csi_kernel
     REAL(c_double) :: csi_dt_out
     CHARACTER(c_char), DIMENSION(length) :: csi_obsfile
     REAL(c_double) :: csi_obsmin
     REAL(c_double) :: csi_modmin
     INTEGER(c_int) :: csi_nz
     REAL(c_double) :: csi_z0
     REAL(c_double) :: csi_z1
     INTEGER(c_int) :: csi_nx
     REAL(c_double) :: csi_lon0
     REAL(c_double) :: csi_lon1
     INTEGER(c_int) :: csi_ny
     REAL(c_double) :: csi_lat0
     REAL(c_double) :: csi_lat1
     CHARACTER(c_char), DIMENSION(length) :: ens_basename
     REAL(c_double) :: ens_dt_out
     CHARACTER(c_char), DIMENSION(length) :: grid_basename
     CHARACTER(c_char), DIMENSION(length) :: grid_kernel
     CHARACTER(c_char), DIMENSION(length) :: grid_gpfile
     REAL(c_double) :: grid_dt_out
     INTEGER(c_int) :: grid_sparse
     INTEGER(c_int) :: grid_nc_level
     INTEGER(c_int), DIMENSION(nqq) :: grid_nc_quant
     INTEGER(c_int) :: grid_stddev
     INTEGER(c_int) :: grid_nz
     REAL(c_double) :: grid_z0
     REAL(c_double) :: grid_z1
     INTEGER(c_int) :: grid_nx
     REAL(c_double) :: grid_lon0
     REAL(c_double) :: grid_lon1
     INTEGER(c_int) :: grid_ny
     REAL(c_double) :: grid_lat0
     REAL(c_double) :: grid_lat1
     INTEGER(c_int) :: grid_type
     CHARACTER(c_char), DIMENSION(length) :: prof_basename
     CHARACTER(c_char), DIMENSION(length) :: prof_obsfile
     INTEGER(c_int) :: prof_nz
     REAL(c_double) :: prof_z0
     REAL(c_double) :: prof_z1
     INTEGER(c_int) :: prof_nx
     REAL(c_double) :: prof_lon0
     REAL(c_double) :: prof_lon1
     INTEGER(c_int) :: prof_ny
     REAL(c_double) :: prof_lat0
     REAL(c_double) :: prof_lat1
     CHARACTER(c_char), DIMENSION(length) :: sample_basename
     CHARACTER(c_char), DIMENSION(length) :: sample_kernel
     CHARACTER(c_char), DIMENSION(length) :: sample_obsfile
     REAL(c_double) :: sample_dx
     REAL(c_double) :: sample_dz
     CHARACTER(c_char), DIMENSION(length) :: stat_basename
     REAL(c_double) :: stat_lon
     REAL(c_double) :: stat_lat
     REAL(c_double) :: stat_r
     REAL(c_double) :: stat_t0
     REAL(c_double) :: stat_t1
     CHARACTER(c_char), DIMENSION(length) :: vtk_basename
     REAL(c_double) :: vtk_dt_out
     INTEGER(c_int) :: vtk_stride
     REAL(c_double) :: vtk_scale
     REAL(c_double) :: vtk_offset
     INTEGER(c_int) :: vtk_sphere 
  END TYPE ctl_t
  
  TYPE, bind(c) :: met_t
     REAL(c_double) :: time
     INTEGER(c_int) :: nx
     INTEGER(c_int) :: ny
     INTEGER(c_int) :: np
     INTEGER(c_int) :: npl
     REAL(c_double), DIMENSION(ex) :: lon
     REAL(c_double), DIMENSION(ey) :: lat
     REAL(c_double), DIMENSION(ep) :: p
     REAL(c_double), DIMENSION(ep) :: hybrid
     REAL(c_float), DIMENSION(ey,ex) :: ps
     REAL(c_float), DIMENSION(ey,ex) :: ts
     REAL(c_float), DIMENSION(ey,ex) :: zs
     REAL(c_float), DIMENSION(ey,ex) :: us
     REAL(c_float), DIMENSION(ey,ex) :: vs
     REAL(c_float), DIMENSION(ey,ex) :: lsm
     REAL(c_float), DIMENSION(ey,ex) :: sst
     REAL(c_float), DIMENSION(ey,ex) :: pbl
     REAL(c_float), DIMENSION(ey,ex) :: pt
     REAL(c_float), DIMENSION(ey,ex) :: tt
     REAL(c_float), DIMENSION(ey,ex) :: zt
     REAL(c_float), DIMENSION(ey,ex) :: h2ot
     REAL(c_float), DIMENSION(ey,ex) :: pct
     REAL(c_float), DIMENSION(ey,ex) :: pcb
     REAL(c_float), DIMENSION(ey,ex) :: cl
     REAL(c_float), DIMENSION(ey,ex) :: plcl
     REAL(c_float), DIMENSION(ey,ex) :: plfc
     REAL(c_float), DIMENSION(ey,ex) :: pel
     REAL(c_float), DIMENSION(ey,ex) :: cape
     REAL(c_float), DIMENSION(ey,ex) :: cin
     REAL(c_float), DIMENSION(ey,ex) :: o3c
     REAL(c_float), DIMENSION(ep,ey,ex) :: z
     REAL(c_float), DIMENSION(ep,ey,ex) :: t
     REAL(c_float), DIMENSION(ep,ey,ex) :: u
     REAL(c_float), DIMENSION(ep,ey,ex) :: v
     REAL(c_float), DIMENSION(ep,ey,ex) :: w
     REAL(c_float), DIMENSION(ep,ey,ex) :: pv
     REAL(c_float), DIMENSION(ep,ey,ex) :: h2o
     REAL(c_float), DIMENSION(ep,ey,ex) :: o3
     REAL(c_float), DIMENSION(ep,ey,ex) :: lwc
     REAL(c_float), DIMENSION(ep,ey,ex) :: rwc
     REAL(c_float), DIMENSION(ep,ey,ex) :: iwc
     REAL(c_float), DIMENSION(ep,ey,ex) :: swc
     REAL(c_float), DIMENSION(ep,ey,ex) :: cc
     REAL(c_float), DIMENSION(ep,ey,ex) :: pl
     REAL(c_float), DIMENSION(ep,ey,ex) :: ul
     REAL(c_float), DIMENSION(ep,ey,ex) :: vl
     REAL(c_float), DIMENSION(ep,ey,ex) :: wl
     REAL(c_float), DIMENSION(ep,ey,ex) :: zetal
     REAL(c_float), DIMENSION(ep,ey,ex) :: zeta_dotl
  END TYPE met_t

END MODULE mptrac_struct

! MPTRAC functions
MODULE mptrac_func
  INTERFACE

     SUBROUTINE mptrac_get_met(ctl, clim, t, met0, met1) &
          bind(c,name='get_met')
       USE iso_c_binding
       USE mptrac_struct, ONLY : ctl_t, clim_t, met_t
       IMPLICIT NONE
       TYPE(ctl_t), INTENT(in), TARGET :: ctl
       TYPE(clim_t), INTENT(in), TARGET :: clim
       REAL(c_double), INTENT(in), VALUE :: t
       TYPE(met_t), INTENT(inout), POINTER :: met0, met1
     END SUBROUTINE mptrac_get_met

     SUBROUTINE mptrac_module_advect(ctl, met0, met1, atm, dt) &
          bind(c,name='module_advect')
       USE iso_c_binding
       USE mptrac_struct, ONLY : ctl_t, met_t, atm_t, npp
       IMPLICIT NONE
       TYPE(ctl_t), INTENT(in), TARGET :: ctl
       TYPE(met_t), INTENT(in), TARGET :: met0, met1
       TYPE(atm_t), INTENT(inout), TARGET :: atm
       REAL(c_double), DIMENSION(npp) :: dt
     END SUBROUTINE mptrac_module_advect

     SUBROUTINE mptrac_module_timesteps(ctl, met0, atm, dt, t) &
          bind(c,name='module_timesteps')
       USE iso_c_binding
       USE mptrac_struct, ONLY : ctl_t, met_t, atm_t, npp
       TYPE(ctl_t), INTENT(in), TARGET :: ctl
       TYPE(met_t), INTENT(in), TARGET :: met0
       TYPE(atm_t), INTENT(in), TARGET :: atm
       REAL(c_double), DIMENSION(npp) :: dt
       REAL(c_double), INTENT(in), VALUE :: t
     END SUBROUTINE mptrac_module_timesteps

     SUBROUTINE mptrac_module_timesteps_init(ctl, atm) &
          bind(c,name='module_timesteps_init')
       USE iso_c_binding
       USE mptrac_struct, ONLY : ctl_t, atm_t
       TYPE(ctl_t), INTENT(in), TARGET :: ctl
       TYPE(atm_t), INTENT(out), TARGET :: atm
     END SUBROUTINE mptrac_module_timesteps_init

     SUBROUTINE mptrac_read_atm(filename,ctl,atm) &
          bind(c,name='read_atm')
       USE iso_c_binding
       USE mptrac_struct, ONLY : ctl_t, atm_t
       IMPLICIT NONE
       CHARACTER(c_char), INTENT(in) :: filename
       TYPE(ctl_t), INTENT(in), TARGET :: ctl
       TYPE(atm_t), INTENT(out), TARGET :: atm
     END SUBROUTINE mptrac_read_atm

     SUBROUTINE mptrac_read_clim(ctl,clim) &
          bind(c,name='read_clim')
       USE iso_c_binding
       USE mptrac_struct, ONLY : ctl_t, clim_t, clim_ts_t, clim_zm_t, clim_photo_t
       IMPLICIT NONE
       TYPE(ctl_t), INTENT(in), TARGET :: ctl
       TYPE(clim_t), INTENT(out), TARGET :: clim
     END SUBROUTINE mptrac_read_clim

     SUBROUTINE mptrac_read_ctl(filename,argc,argv,ctl) &
          bind(c,name='read_ctl')
       USE iso_c_binding
       USE mptrac_struct, ONLY : ctl_t
       IMPLICIT NONE
       CHARACTER(c_char), INTENT(in) :: filename
       INTEGER(c_int), VALUE :: argc
       TYPE(c_ptr), DIMENSION(5) :: argv
       TYPE(ctl_t), INTENT(out), TARGET :: ctl
     END SUBROUTINE mptrac_read_ctl

     SUBROUTINE mptrac_read_met(filename,ctl,clim,met,atm) &
          bind(c,name='read_met')
       USE iso_c_binding
       USE mptrac_struct, ONLY : ctl_t, clim_t, met_t, atm_t
       IMPLICIT NONE
       CHARACTER(c_char), DIMENSION(*), INTENT(in) :: filename
       TYPE(ctl_t), INTENT(in), TARGET :: ctl
       TYPE(clim_t), INTENT(in), TARGET :: clim
       TYPE(atm_t), INTENT(in), TARGET :: atm
       TYPE(met_t), INTENT(out), TARGET :: met
     END SUBROUTINE mptrac_read_met

     SUBROUTINE mptrac_write_output(dirname, ctl, met0, met1, atm, t) &
          bind(c,name='write_output')
       USE iso_c_binding
       USE mptrac_struct, ONLY : ctl_t, met_t, atm_t
       CHARACTER(c_char), INTENT(in) :: dirname
       TYPE(ctl_t), INTENT(in), TARGET :: ctl
       TYPE(met_t), INTENT(in), TARGET :: met0, met1
       TYPE(atm_t), INTENT(in), TARGET :: atm
       REAL(c_double), INTENT(in), VALUE :: t
     END SUBROUTINE mptrac_write_output

  END INTERFACE
END MODULE mptrac_func

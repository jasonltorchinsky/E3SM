!
! Theta-l: Namelist for Held-Suarez
! https://e3sm.atlassian.net/wiki/spaces/DOC/pages/1044644202/EAM+s+HOMME+dycore
!_______________________________________________________________________
&ctl_nl
  nthreads               = 1
  partmethod             = 4               ! Mesh Parition Method: 4 = Space Filling Curve
  topology               = "cube"          ! Mesh Type: Cubed Sphere
  test_case              = "held_suarez"   ! Test Identifier
  sub_case               = 0               ! 0 = Original (Held & Suarez 1994), 1 = Stratosphere Modification (Polvani & Kushner 2002)
  ne                     = 8               ! Number of Elements Along Edge of Each Cube Face
  qsize                  = 0               ! Number of Tracer Fields
  ndays                  = 1200
  statefreq              = 1440            ! Number of Steps Between Screen Dumps
  restartfreq            = -1              ! Don't write restart files if < 0
  restartdir             = "./restart/"
  restartfile            = "restart/R0001"
  runtype                = 0               ! 0 = New Run
  tstep                  = 900             ! Largest Timestep in Seconds
  integration            = 'explicit'      ! Time Integration - 'explicit', 'implicit'
  tstep_type             = 9               ! IMEX Scheme (Default: 9)
  rsplit                 = -1              ! Remapping Frequency (Default: 6)
  qsplit                 = -1
  dt_tracer_factor       = 1
  dt_remap_factor        = 2
  nu                     = 3.4e-8
  nu_div                 = -1
  nu_s                   = -1
  nu_p                   = -1
  nu_top                 = 2.5e5
  limiter_option         = 9
  vert_remap_q_alg       = 10
  hypervis_order         = 2               ! 2 = Hyperviscosity
  hypervis_subcycle_q    = 1               ! Default: Equal to dt_tracer_factor
  hypervis_subcycle      = 1               ! 1 = No Hyperviscosity Subcycling
  hypervis_subcycle_tom  = 1               !
  hypervis_scaling       = 3.0
  hv_ref_profiles        = 2
  hv_theta_correction    = 1
  moisture               = 'wet'
  pgrad_correction       = 1
  se_ftype               = 2               ! Default: 2
  theta_hydrostatic_mode = .false.         ! (Non-)Hydrostatic Mode
  theta_advect_form      = 1
/

&vert_nl
  vfile_mid     = "../vcoord/sabm-128.ascii"
  vfile_int     = "../vcoord/sabi-128.ascii"
/

&analysis_nl
  infilenames       = ""              ! Topography file path
  output_dir        = "./movies/"     ! Destination directory for NetCDF file
  output_timeunits  = 2,              ! 0 = Timesteps, 1 = Days, 2 = Hours, 3 = Seconds
  output_frequency  = 24,             ! Every N time units
  output_start_time = 0,              ! Output start time
  output_end_time   = 99999,          ! Output end time
  output_varnames1  = 'u','v','T','zeta','ps','p','pnh','w'   ! Variables to write to file
  interp_type       = 1               ! 0 = Native grid, 1 = Bilinear
  interp_gridtype   = 1
  interp_lon0       = -180.0          ! Shift longitude range to [-180, +180)
  output_type       = 'netcdf'        ! NetCDF or PNetCDF
  num_io_procs      = 16
/

&prof_inparm
  profile_outpe_num   = 100
  profile_single_file = .true.
/

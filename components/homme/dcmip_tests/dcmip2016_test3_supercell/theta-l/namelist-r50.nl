!
! Theta: Namelist for DCMIP 2016 Test 3 - Supercell Storm (Small Planet X=120)
!_______________________________________________________________________
&ctl_nl
  nthreads          = 1
  partmethod        = 4                         ! Mesh Parition Method: 4 = Space Filling Curve
  topology          = "cube"                    ! Mesh Type: Cubed Sphere
  test_case         = "dcmip2016_test3"         ! Test Identifier
  ne                = 60                        ! Number of Elements per Cube Face
  qsize             = 3                         ! Number of Tracer Fields
  nmax              = 7200                      ! 7200s(120min)/tstep
  statefreq         = 360                       ! Number of Steps Between Screen Dumps
  restartfreq       = -1                        ! Don't write restart files if < 0
  runtype           = 0                         ! 0 = New Run
  tstep             = 1.0                       ! Largest Timestep in Seconds
  integration       = 'explicit'                ! Explicit Time Integration
  tstep_type        = 9                         ! IMEX Scheme (Default: 7, but is worse than 9?)
  rsplit            = 2                         ! Remapping Frequency (Default: 6)
  qsplit            = 2                         ! dt_tracer <= 3s  
  nu                = 6.3e7                     ! default= 1e15/(120)^3 *(ne30/ne60)**3.2
  nu_s              = 0
  nu_p              = 0
  nu_q              = 0
  nu_top            = 0                         ! 2.5e5/(120)^(1)
  vert_remap_q_alg  = -1
  limiter_option    = 4
  dcmip16_mu        = 500.0d0                   ! Additional Uniform Viscosity
  dcmip16_mu_s      = 1500.0d0
  hypervis_order    = 2                         ! 2 = Hyperviscosity
  hypervis_subcycle = 1                         ! 1 = No Hyperviscosity Subcycling
  rearth            = 53133                     ! 6.376E6  / 120
  omega             = 0
  se_ftype          = 0
  moisture          = 'wet'
  theta_hydrostatic_mode = .false.
  dcmip16_prec_type = 0                         ! 0 = Kessler
  dcmip16_pbl_type  = 0                         ! 0 = Basic
  dcmip16_phys_type = 1                         ! 0 = Isochoric, 1 = Isobaric 
/

&vert_nl
  vanalytic         = 1                         ! Set vcoords in Initialization Routine
/
&analysis_nl
  output_dir        = "./movies/"               ! Destination directory for NetCDF file
  output_timeunits  = 3                         ! 0 = Timesteps, 1 = Days, 2 = Hours, 3 = Seconds
  output_frequency  = 600                       ! Every N time units
  output_varnames1  ='T','p','ps','pnh','geo','u','v','w','Th','Q','Q2','Q3','precl'   ! variables to write to file
!  interp_nlon       = 360
!  interp_nlat       = 181
  interp_gridtype   = 1
  interp_type       = 1                         ! 0 = Native grid, 1 = Bilinear
  interp_lon0       = -180.0                    ! Shift longitude range to [-180, +180)
  output_type       ='netcdf'                   ! NetCDF or PNetCDF
  num_io_procs      = 16         
/
&prof_inparm
  profile_outpe_num   = 100
  profile_single_file	= .true.
/

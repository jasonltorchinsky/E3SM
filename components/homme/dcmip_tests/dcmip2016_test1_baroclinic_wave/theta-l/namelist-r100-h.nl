!
! Theta-l: Namelist for DCMIP 2016 Test 1 -  Moist Baroclinic Wave
!_______________________________________________________________________
&ctl_nl
  nthreads          = -1                        ! Use OMP_NUM_THREADS
  partmethod        = 4                         ! Mesh Partition Method: 4 = Space Filling Curve
  topology          = "cube"                    ! Mesh Type: Cubed Sphere
  test_case         = "dcmip2016_test1"         ! Test Identifier
  ne                = 30                        ! Number of Elements per Cube Face
  qsize             = 6                         ! Number of Tracer Fields
  ndays             = 30
  statefreq         = 72                        ! Number of Steps Between Screen Dumps
  restartfreq       = -1                        ! Don't write restart files if < 0
  runtype           = 0                         ! 0 = New Run
  tstep             = 300                       ! Largest Timestep in Seconds
  integration       = 'explicit'                ! Explicit Time Integration
  tstep_type        = 9                         ! IMEX Scheme (Default: 7, but is worse than 9?)
  rsplit            = 2                         ! Remapping Frequency (Default: 6)
  qsplit            = 1
  nu                = 1e15                      ! Default: 1e15*(ne30/ne30)**3.2 = 1e15
  nu_s              = 1e15
  nu_p              = 1e15  
  nu_top            = 0                         ! Default: 2.5e5
  limiter_option    = 9
  hypervis_order    = 2                         ! 2 = Hyperviscosity
  hypervis_subcycle = 1                         ! 1 = No Hyperviscosity Subcycling
  moisture          = 'wet'
  theta_hydrostatic_mode = .true.
  dcmip16_prec_type = 1                         ! 0 = Kessler,   1 = Reed-Jablonowski
  dcmip16_pbl_type  = -1                        ! 0 = Basic,     1 = Bryan [Bugged], -1 = None
  dcmip16_phys_type = 1                         ! 0 = Isochoric, 1 = Isobaric
/
&vert_nl
  vfile_mid         = "../vcoord/camm-30.ascii"
  vfile_int         = "../vcoord/cami-30.ascii"
/
&analysis_nl
!  output_prefix     = "r100-"
  output_dir        = "./movies/"               ! Destination directory for NetCDF file
  output_timeunits  = 2,                        ! 0 = Timesteps, 1 = Days, 2 = Hours, 3 = Seconds
  output_frequency  = 24                        ! Every N hours
!  output_varnames1  ='T','ps','pnh','geo','u','v','w','omega','Th','Q','Q2','Q3','Q4','Q5','rho','precl','zeta'   ! variables to write to file
  output_varnames1  ='T','ps','pnh','geo','u','Th','Q','Q2','Q3','Q4','Q5','precl','zeta'   ! variables to write to file
  interp_type       = 1                         ! 0 = Native grid, 1 = Bilinear
  output_type       ='netcdf'                   ! NetCDF or PNetCDF
  num_io_procs      = 16
!  interp_nlon       = 360
!  interp_nlat       = 181
  interp_gridtype   = 1
/
&prof_inparm
  profile_outpe_num   = 100
  profile_single_file	= .true.
/

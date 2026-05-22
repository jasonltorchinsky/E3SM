#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module held_suarez_mod

  use coordinate_systems_mod, only: spherical_polar_t
  use dimensions_mod,         only: nlev,np,qsize,nlevp
  use element_mod,            only: element_t
  use element_state,          only: timelevels
  use element_ops,            only: set_thermostate, get_temperature
  use hybrid_mod,             only: hybrid_t
  use hybvcoord_mod,          only: hvcoord_t
  use kinds,                  only: real_kind,iulog
  use parallel_mod,           only: abortmp
  use physical_constants,     only: p0,kappa,g,dd_pi,Rgas,rearth0
  use physics_mod,            only: prim_condense
  use time_mod,               only: secpday
#ifndef HOMME_WITHOUT_PIOLIBRARY
  use common_io_mod,          only: infilenames
#endif

implicit none
private

  ! From (Held & Suarez 1994), doi: 10.1175/1520-0477(1994)075<1825:APFTIO>2.0.CO;2
  real (kind=real_kind), public, parameter :: hs94_sigma_b  = 0.70D0
  real (kind=real_kind), public, parameter :: hs94_k_a      = 1.0D0/(40.0D0*secpday)
  real (kind=real_kind), public, parameter :: hs94_k_f      = 1.0D0/(1.0D0*secpday)
  real (kind=real_kind), public, parameter :: hs94_k_s      = 1.0D0/(4.0D0*secpday)
  real (kind=real_kind), public, parameter :: hs94_dT_y     = 60.0D0
  real (kind=real_kind), public, parameter :: hs94_dtheta_z = 10.0D0
  
  ! U.S. Standard Atmosphere 1976 parameters
  ! Derived from (NOAA 1976), doi: unknown
  real (kind=real_kind), public, parameter :: us76_M_0 = 28.9644D0 ! Mean molecular weight of dry air at sea-level [kg kmol^{-1}]
  real (kind=real_kind), public, parameter :: us76_R_star = 8.31432D3 ! Universal gas constant [N m K^{-1} kmol^{-1}]
  real (kind=real_kind), public, parameter :: us76_H_b(0:7) &
    = (/ 0.0D3, 11.0D3, 20.0D3, 32.0D3, 47.0D3, 51.0D3, 71.0D3, 84.8520D3 /) ! Geopotential height at reference levels [m']
  real (kind=real_kind), public, parameter :: us76_LM_b(0:6) &
    = (/ -6.5D-3, 0.0D-3, 1.0D-3, 2.8D-3, 0.0D-3, -2.8D-3, -2.0D-3 /) ! Molecular-scale temperature gradient between reference levels [K m'^{-1}]
  real (kind=real_kind), public, parameter :: us76_TM_b(0:7) &
    = (/ 288.15D0, 216.65D0, 216.65D0, 228.65D0, 270.65D0, 270.65D0, 214.65D0, 186.946D0 /) ! Molecular-scale temperature at reference levels [K]
  real (kind=real_kind), public, parameter :: us76_p_b(0:7) &
    = (/ 1.01325000D5, 2.26320640D4, 5.47488867D3, 8.68018685D2, 1.10906306D2, &
       6.69388731D3, 3.95642043D4, 3.73383590D5 /) ! Pressure at reference levels [Pa]

  ! From (Polvani & Kushner 2002), doi: 10.1029/2001GL014284
  ! Typos corrected in, e.g., (Kushner & Polvani 2006), doi: 10.7916/D8G451FW
  ! Modification from (Held & Suarez 1994)
  real (kind=real_kind), public, parameter :: pk02_sigma_b   = 0.70D0 ! [N/A]
  real (kind=real_kind), public, parameter :: pk02_k_a       = 1.0D0/(40.0D0*secpday) ! [day^{-1}] => [sec^{-1}]
  real (kind=real_kind), public, parameter :: pk02_k_s       = 1.0D0/(4.0D0*secpday)  ! [day^{-1}] => [sec^{-1}]
  real (kind=real_kind), public, parameter :: pk02_k_max     = 1.0D0/(2.0D0*secpday)  ! [day^{-1}] => [sec^{-1}]
  real (kind=real_kind), public, parameter :: pk02_p_sp      = 0.5D2  ! [Pa]
  real (kind=real_kind), public, parameter :: pk02_phi_0     = -50.0D0*(dd_pi/180.0D0)  ! [deg latitude] => [radians]
  real (kind=real_kind), public, parameter :: pk02_delta_phi = 10.0D0*(dd_pi/180.0D0)   ! [deg latitude] => [radians]
  real (kind=real_kind), public, parameter :: pk02_p_T       = 100.0D2  ! [Pa]
  real (kind=real_kind), public, parameter :: pk02_T_T       = 216.65D0 ! [K]
  real (kind=real_kind), public, parameter :: pk02_T_0       = 315.0D0  ! [K]
  real (kind=real_kind), public, parameter :: pk02_p_0       = 1000.0D2 ! [Pa]
  real (kind=real_kind), public, parameter :: pk02_kappa     = 2.0D0 / 7.0D0 ! [N/A]
  real (kind=real_kind), public, parameter :: pk02_delta_y   = 60.0D0   ! [K]
  real (kind=real_kind), public, parameter :: pk02_delta_z   = 10.0D0   ! [K]
  real (kind=real_kind), public, parameter :: pk02_epsilon   = 10.0D0   ! [K]
  real (kind=real_kind), public, parameter :: pk02_gamma     = 2.0D-3   ! [K m^{-1}]
   
  public :: hs_init_state
  public :: hs_forcing

contains

  subroutine hs_forcing(elemin,hvcoord,nm1,nm1_Q,dt,sub_case)

    type (element_t)                  :: elemin
    type (hvcoord_t)                  :: hvcoord
    integer                           :: nm1,nm1_Q  ! timelevel to use
    real (kind=real_kind)             :: dt
    integer,              intent(in)  :: sub_case

    select case (sub_case)
    case(0)
      call hs0_forcing(elemin,hvcoord,nm1,nm1_Q,dt)
    case(1)
      call hs1_forcing(elemin,hvcoord,nm1,nm1_Q,dt)
    case default
      call abortmp('invalid forcing sub_case: only subcase = 0 supported')
    end select

  end subroutine hs_forcing

  subroutine hs0_forcing(elemin,hvcoord,nm1,nm1_Q,dt)

    type (element_t)      :: elemin
    type (hvcoord_t)      :: hvcoord
    integer               :: nm1,nm1_Q  ! timelevel to use
    real (kind=real_kind) :: dt

    ! local
    real (kind=real_kind) :: pmid,r0,r1,dtf_q,dp,rdp,FQ
    real (kind=real_kind) :: psfrc(np,np)
    real (kind=real_kind) :: temperature(np,np,nlev)
    real (kind=real_kind) :: v(np,np,3,nlev)
    real (kind=real_kind) :: fv(np,np,3,nlev)
    integer               :: i,j,k,q

    dtf_q = dt
    call get_temperature(elemin,temperature,hvcoord,nm1)
        
    do j=1,np
       do i=1,np
          psfrc(i,j) = (elemin%state%ps_v(i,j,nm1))
       end do
    end do

    elemin%derived%FT(:,:,:) = elemin%derived%FT(:,:,:) + &
         hs0_T_forcing(hvcoord,psfrc(1,1),temperature,elemin%spherep,np,nlev)

    v(:,:,1:2,:) = elemin%state%v(:,:,1:2,:,nm1)
#if ( defined MODEL_THETA_L ) 
    v(:,:,3,:) = elemin%state%w_i(:,:,1:nlev,nm1)  ! dont apply at surface
#else
    v(:,:,3,:) = 0
#endif

    fv = hs0_v_forcing(hvcoord,psfrc(1,1),v,np,nlev)

#if ( defined MODEL_THETA_L ) 
    elemin%derived%FM(:,:,1:3,:) = elemin%derived%FM(:,:,1:3,:) + fv(:,:,1:3,:)
#else
    elemin%derived%FM(:,:,1:2,:) = elemin%derived%FM(:,:,1:2,:) + fv(:,:,1:2,:)
#endif

    if (qsize>=1) then
       ! HS with tracer  (Galewsky type forcing, with flux of  2.3e-5 kg/m^2/s
       ! MASS in kg/m^2   = < Q dp_in_Pa / g >   
       ! flux in kg/m^2/s = < FQ dp_in_Pa / g >   
       ! We want < FQ dp_in_Pa / g > = 2.3e-5  so:  FQ = 2.3e-5*g/dp_in_Pa 

       ! lowest layer thickness, in Pa
       dp = ( hvcoord%hyai(nlev+1) - hvcoord%hyai(nlev) ) + &
               ( hvcoord%hybi(nlev+1) - hvcoord%hybi(nlev) )*1000*100
       rdp = 1./ dp
       q = 1
       do j = 1,np
          do i = 1,np
             FQ = rdp * g * 2.3E-5 * COS(elemin%spherep(i,j)%lat)**2
             elemin%derived%FQ(i,j,nlev,q) = elemin%derived%FQ(i,j,nlev,q) + FQ
          end do
       end do

       do j=1,np
          do i=1,np
             do k=1,nlev
                pmid = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*(elemin%state%ps_v(i,j,nm1))
                r0 = elemin%state%Q(i,j,k,q)
                r1 = r0
                call Prim_Condense(r1,temperature(i,j,k),pmid)
                elemin%derived%FQ(i,j,k,q) = elemin%derived%FQ(i,j,k,q) + &
                     (r1-r0)/(dtf_q)
             end do
          end do
       end do
    end if

  end subroutine hs0_forcing

  function hs0_v_forcing(hvcoord,ps,v,npts,nlevels) result(hs_v_frc)

    integer,               intent(in) :: npts
    integer,               intent(in) :: nlevels
    type (hvcoord_t),      intent(in) :: hvcoord
    real (kind=real_kind), intent(in) :: ps(npts,npts)

    real (kind=real_kind), intent(in) :: v(npts,npts,3,nlevels)
    real (kind=real_kind)             :: hs_v_frc(npts,npts,3,nlevels)

    ! Local variables
    integer i,j,k
    real (kind=real_kind) :: k_v
    real (kind=real_kind) :: p,eta

    do k=1,nlevels
       do j=1,npts
          do i=1,npts
             p = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*ps(i,j)
             eta = hvcoord%hyam(k) + hvcoord%hybm(k)
             k_v = hs94_k_f*MAX(0.0_real_kind,(eta - hs94_sigma_b )/(1.0_real_kind - hs94_sigma_b))
             hs_v_frc(i,j,1,k) = -k_v*v(i,j,1,k)
             hs_v_frc(i,j,2,k) = -k_v*v(i,j,2,k)

             eta = hvcoord%hyai(k) + hvcoord%hybi(k)
             k_v = hs94_k_f*MAX(0.0_real_kind,(eta - hs94_sigma_b )/(1.0_real_kind - hs94_sigma_b))
             hs_v_frc(i,j,3,k) = -k_v*v(i,j,3,k)
          end do
       end do
    end do

  end function hs0_v_forcing

  function hs0_T_forcing(hvcoord,ps,T,sphere,npts,nlevels) result(hs_T_frc)

    integer, intent(in) :: npts
    integer, intent(in) :: nlevels

    type (hvcoord_t),         intent(in) :: hvcoord
    real (kind=real_kind),    intent(in) :: ps(npts,npts)
    real (kind=real_kind),    intent(in) :: T(npts,npts,nlevels)
    type (spherical_polar_t), intent(in) :: sphere(npts,npts)

    real (kind=real_kind)                :: hs_T_frc(npts,npts,nlevels)

    ! Local variables
    real (kind=real_kind) :: p,logprat,pratk,Teq
    real (kind=real_kind) :: logps0,etam
    real (kind=real_kind) :: lat,snlat

    real (kind=real_kind) :: k_t(npts,npts)
    real (kind=real_kind) :: snlatsq(npts,npts)
    real (kind=real_kind) :: cslatsq(npts,npts)

    real (kind=real_kind) :: rec_one_minus_sigma_b

    integer i,j,k

    logps0 = LOG(hvcoord%ps0)

    do j=1,npts
       do i=1,npts
         snlat        = SIN(sphere(i,j)%lat)
         snlatsq(i,j) = snlat*snlat
         cslatsq(i,j) = 1.0D0 - snlatsq(i,j)
       end do
    end do

    rec_one_minus_sigma_b = 1.0D0/(1.0D0 - hs94_sigma_b)

    do k=1,nlevels
       do j=1,npts
          do i=1,npts
             p         = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*ps(i,j)
             logprat   = LOG(p)-logps0
             pratk     = EXP(kappa*(logprat))
             etam      = hvcoord%hyam(k) + hvcoord%hybm(k)

             k_t(i,j)= hs94_k_a + (hs94_k_s-hs94_k_a)*cslatsq(i,j)*cslatsq(i,j)* &
                       MAX(0.0D0,(etam - hs94_sigma_b)/(1.0D0 - hs94_sigma_b))
             Teq     = MAX(200.0D0,(315.0D0 - hs94_dT_y*snlatsq(i,j) - hs94_dtheta_z*logprat*cslatsq(i,j))*pratk)

#if 0
             ! ======================================
             ! This is a smooth forcing 
             ! for debugging purposes only...
             ! ======================================

             k_t(i,j) = hs94_k_a 
             pratk    = EXP(0.081*(logprat))
             Teq      = (315.0D0 - hs94_dT_y*snlatsq(i,j))*pratk
#endif
             hs_T_frc(i,j,k)= -k_t(i,j)*(T(i,j,k)-Teq)
          end do
       end do
    end do
      
  end function hs0_T_forcing

  subroutine hs1_forcing(elemin,hvcoord,nm1,nm1_Q,dt)

    type (element_t)      :: elemin
    type (hvcoord_t)      :: hvcoord
    integer               :: nm1,nm1_Q  ! timelevel to use
    real (kind=real_kind) :: dt

    ! local
    real (kind=real_kind) :: pmid,r0,r1,dtf_q,dp,rdp,FQ
    real (kind=real_kind) :: psfrc(np,np)
    real (kind=real_kind) :: temperature(np,np,nlev)
    real (kind=real_kind) :: v(np,np,3,nlev)
    real (kind=real_kind) :: fv(np,np,3,nlev)
    integer               :: i,j,k,q

    dtf_q = dt
    call get_temperature(elemin,temperature,hvcoord,nm1)
        
    do j = 1,np
       do i = 1,np
          psfrc(i,j) = (elemin%state%ps_v(i,j,nm1))
       end do
    end do

    elemin%derived%FT(:,:,:) = elemin%derived%FT(:,:,:) & 
         + hs1_T_forcing(hvcoord,psfrc(1,1),temperature,elemin%spherep,np,nlev)

    v(:,:,1:2,:) = elemin%state%v(:,:,1:2,:,nm1)
#if ( defined MODEL_THETA_L ) 
    v(:,:,3,:) = elemin%state%w_i(:,:,1:nlev,nm1)  ! dont apply at surface
#else
    v(:,:,3,:) = 0
#endif

    fv = hs1_v_forcing(hvcoord,psfrc(1,1),v,np,nlev)

#if ( defined MODEL_THETA_L ) 
    elemin%derived%FM(:,:,1:3,:) = elemin%derived%FM(:,:,1:3,:) + fv(:,:,1:3,:)
#else
    elemin%derived%FM(:,:,1:2,:) = elemin%derived%FM(:,:,1:2,:) + fv(:,:,1:2,:)
#endif

    if (qsize >= 1) then
       ! HS with tracer  (Galewsky type forcing, with flux of  2.3e-5 kg/m^2/s
       ! MASS in kg/m^2   = < Q dp_in_Pa / g >   
       ! flux in kg/m^2/s = < FQ dp_in_Pa / g >   
       ! We want < FQ dp_in_Pa / g > = 2.3e-5  so:  FQ = 2.3e-5*g/dp_in_Pa 

       ! lowest layer thickness, in Pa
       dp = ( hvcoord%hyai(nlev+1) - hvcoord%hyai(nlev) ) + &
               ( hvcoord%hybi(nlev+1) - hvcoord%hybi(nlev) )*1000*100
       rdp = 1./ dp
       q = 1
       do j = 1,np
          do i = 1,np
             FQ = rdp * g * 2.3E-5 * COS(elemin%spherep(i,j)%lat)**2
             elemin%derived%FQ(i,j,nlev,q) = elemin%derived%FQ(i,j,nlev,q) + FQ
          end do
       end do

       do j = 1,np
          do i = 1,np
             do k = 1,nlev
                pmid = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*(elemin%state%ps_v(i,j,nm1))
                r0 = elemin%state%Q(i,j,k,q)
                r1 = r0
                call Prim_Condense(r1,temperature(i,j,k),pmid)
                elemin%derived%FQ(i,j,k,q) = elemin%derived%FQ(i,j,k,q) + &
                     (r1 - r0) / (dtf_q)
             end do
          end do
       end do
    end if

  end subroutine hs1_forcing

  function hs1_v_forcing(hvcoord,ps,v,npts,nlevels) result(hs_v_frc)

    integer,               intent(in) :: npts
    integer,               intent(in) :: nlevels
    type (hvcoord_t),      intent(in) :: hvcoord
    real (kind=real_kind), intent(in) :: ps(npts,npts)

    real (kind=real_kind), intent(in) :: v(npts,npts,3,nlevels)
    real (kind=real_kind)             :: hs_v_frc(npts,npts,3,nlevels)

    ! Local variables
    integer i,j,k
    real (kind=real_kind) :: k_sp
    real (kind=real_kind) :: p

    do k = 1,nlevels
       do j = 1,npts
          do i = 1,npts
             ! Damp horizontal velocities
             p = hvcoord%hyam(k) * hvcoord%ps0 + hvcoord%hybm(k) * ps(i,j) ! Pressure at midpoint [Pa]
             if (p >= pk02_p_sp) then
               k_sp = 0.0D0
             else ! p < pk02_p_sp
               k_sp = pk02_k_max * ((pk02_p_sp - p) / pk02_p_sp)**2
             end if
             hs_v_frc(i,j,1,k) = -k_sp * v(i,j,1,k)
             hs_v_frc(i,j,2,k) = -k_sp * v(i,j,2,k)

             ! Damp vertical velocity
             p = hvcoord%hyai(k) * hvcoord%ps0 + hvcoord%hybi(k) * ps(i,j) ! Pressure at interface [Pa]
             if (p >= pk02_p_sp) then
               k_sp = 0.0D0
             else ! p < pk02_p_sp
               k_sp = pk02_k_max * ((pk02_p_sp - p) / pk02_p_sp)**2
             end if
             hs_v_frc(i,j,3,k) = -k_sp * v(i,j,3,k)
          end do
       end do
    end do

  end function hs1_v_forcing

  function hs1_T_forcing(hvcoord,ps,T,sphere,npts,nlevels) result(hs_T_frc)

    integer, intent(in) :: npts
    integer, intent(in) :: nlevels

    type (hvcoord_t),         intent(in) :: hvcoord
    real (kind=real_kind),    intent(in) :: ps(npts,npts)
    real (kind=real_kind),    intent(in) :: T(npts,npts,nlevels)
    type (spherical_polar_t), intent(in) :: sphere(npts,npts)

    real (kind=real_kind) :: hs_T_frc(npts,npts,nlevels)

    ! Local variables
    real (kind=real_kind) :: p,logprat,pratk,delta_T,Teq
    real (kind=real_kind) :: logps0,etam
    real (kind=real_kind) :: lat

    real (kind=real_kind) :: k_t(npts,npts)
    real (kind=real_kind) :: snlat(npts,npts),snlatsq(npts,npts),cslatsq(npts,npts)
    real (kind=real_kind) :: W(npts,npts),one_minus_W(npts,npts)

    real (kind=real_kind) :: rec_one_minus_sigma_b
    real (kind=real_kind) :: T_US,T_PV

    integer i,j,k

    ! Comments reflect notation in (Polvani & Kushner 2002), doi: 10.1029/2001GL014284
    logps0 = LOG(pk02_p_0) ! log(p_0)

    do j = 1,npts
       do i = 1,npts
         snlat(i,j) = SIN(sphere(i,j)%lat)      ! sin(phi)
         snlatsq(i,j) = snlat(i,j) * snlat(i,j) ! sin^2(phi)
         cslatsq(i,j) = 1.0D0 - snlatsq(i,j)    ! cos^2(phi)
         W(i,j) = 0.5D0 * (1.0D0 - TANH((sphere(i,j)%lat - pk02_phi_0)/ pk02_delta_phi))
         one_minus_W(i,j) = 1.0D0 - W(i,j)
       end do
    end do

    rec_one_minus_sigma_b = 1.0D0 / (1.0D0 - pk02_sigma_b) ! 1 / (1 - sigma_b)

    do k = 1,nlevels
       etam = hvcoord%hyam(k) + hvcoord%hybm(k) ! sigma
       do j = 1,npts
          do i = 1,npts
             k_t(i,j) = pk02_k_a + (pk02_k_s - hs94_k_a) &
                * MAX(0.0D0, (etam - pk02_sigma_b) / (1.0D0 - pk02_sigma_b)) &
                * cslatsq(i,j) * cslatsq(i,j)
             
             p = hvcoord%hyam(k) * hvcoord%ps0 + hvcoord%hybm(k) * ps(i,j) ! p; pressure at midpoint [Pa]
             if (p >= pk02_p_T) then ! Troposphere equilibrium temperature profile
                logprat = LOG(p) - logps0 ! log(p / p0)
                pratk = EXP(pk02_kappa * logprat) ! (p / p0)^kappa

                delta_T =  pk02_delta_y * snlatsq(i,j) &
                  + pk02_epsilon * snlat(i,j) &
                  + pk02_delta_z * logprat * cslatsq(i,j)
                Teq = MAX(pk02_T_T, (pk02_T_0 - delta_T) * pratk)
             else ! p < pk02_p_T, stratosphere equilibrium temperature profile
               T_US = us_std_T(us_std_z(p))
               T_PV = pk02_T_T * (p / pk02_p_T)**(Rgas * pk02_gamma / g)
               Teq = one_minus_W(i,j) * T_US + W(i,j) * T_PV
             end if

             hs_T_frc(i,j,k) = -k_t(i,j) * (T(i,j,k) - Teq)
          end do
       end do
    end do
      
  end function hs1_T_forcing

  function us_std_z(p) result(z)
    real (kind=real_kind), intent(in)  :: p
    real (kind=real_kind) :: z

    ! local variables
    real (kind=real_kind) :: h

    if ((p <= us76_p_b(0)) .and. (p > us76_p_b(1))) then
         h = us76_H_b(0) + (us76_TM_b(0) / us76_LM_b(0)) &
            * ((us76_p_b(0) / p)**((us76_R_star * us76_LM_b(0)) / (g * us76_M_0)) - 1.0D0)
    else if ((p <= us76_p_b(1)) .and. (p > us76_p_b(2))) then
         h = us76_H_b(1) - ((us76_R_star * us76_TM_b(1)) / (g * us76_M_0)) * LOG(p / us76_p_b(1))
    else if ((p <= us76_p_b(2)) .and. (p > us76_p_b(3))) then
         h = us76_H_b(2) + (us76_TM_b(2) / us76_LM_b(2)) &
            * ((us76_p_b(2) / p)**((us76_R_star * us76_LM_b(2)) / (g * us76_M_0)) - 1.0D0)
    else if ((p <= us76_p_b(3)) .and. (p > us76_p_b(4))) then
         h = us76_H_b(3) + (us76_TM_b(3) / us76_LM_b(3)) &
            * ((us76_p_b(3) / p)**((us76_R_star * us76_LM_b(3)) / (g * us76_M_0)) - 1.0D0)
    else if ((p <= us76_p_b(4)) .and. (p > us76_p_b(5))) then
         h = us76_H_b(4) - ((us76_R_star * us76_TM_b(4)) / (g * us76_M_0)) * LOG(p / us76_p_b(4))
    else if ((p <= us76_p_b(5)) .and. (p > us76_p_b(6))) then
         h = us76_H_b(5) + (us76_TM_b(5) / us76_LM_b(5)) &
            * ((us76_p_b(5) / p)**((us76_R_star * us76_LM_b(5)) / (g * us76_M_0)) - 1.0D0)
    else if ((p <= us76_p_b(6)) .and. (p >= us76_p_b(7))) then
         h = us76_H_b(6) + (us76_TM_b(6) / us76_LM_b(6)) &
            * ((us76_p_b(6) / p)**((us76_R_star * us76_LM_b(6)) / (g * us76_M_0)) - 1.0D0)
    end if

    z = (rearth0 * h) / (rearth0 - h) ! Geometric height [m]

  end function us_std_z

  function us_std_T(z) result (T)
    real (kind=real_kind), intent(in)  :: z
    real (kind=real_kind) :: T

    ! local variables
    real (kind=real_kind) :: h

    h = (rearth0 * z) / (rearth0 + z) ! Geopotential height [m']
    if ((h >= us76_H_b(0)) .and. (h < us76_H_b(1))) then
         T = us76_TM_b(0) + us76_LM_b(0) * (h - us76_H_b(0))
    else if ((h >= us76_H_b(1)) .and. (h < us76_H_b(2))) then
         T = us76_TM_b(1) + us76_LM_b(1) * (h - us76_H_b(1))
    else if ((h >= us76_H_b(2)) .and. (h < us76_H_b(3))) then
         T = us76_TM_b(2) + us76_LM_b(2) * (h - us76_H_b(2))
    else if ((h >= us76_H_b(3)) .and. (h < us76_H_b(4))) then
         T = us76_TM_b(3) + us76_LM_b(3) * (h - us76_H_b(3))
    else if ((h >= us76_H_b(4)) .and. (h < us76_H_b(5))) then
         T = us76_TM_b(4) + us76_LM_b(4) * (h - us76_H_b(4))
    else if ((h >= us76_H_b(5)) .and. (h < us76_H_b(6))) then
         T = us76_TM_b(5) + us76_LM_b(5) * (h - us76_H_b(5))
    else if ((h >= us76_H_b(6)) .and. (h <= us76_H_b(7))) then
         T = us76_TM_b(6) + us76_LM_b(6) * (h - us76_H_b(6))
    end if

  end function us_std_T

  subroutine hs_init_state(elem, hybrid, hvcoord,nets,nete,Tinit,sub_case)

    type(element_t),        intent(inout) :: elem(:)
    type(hybrid_t),         intent(in)    :: hybrid ! hybrid parallel structure
    type (hvcoord_t),       intent(in)    :: hvcoord
    integer,                intent(in)    :: nets
    integer,                intent(in)    :: nete
    real (kind=real_kind),  intent(in)    :: Tinit
    integer,                intent(in)    :: sub_case

    if (hybrid%masterthread) write(iulog,*) 'initializing Held-Suarez primitive equations test sub_case ', sub_case

    select case (sub_case)
    case(0)
      call hs0_init_state(elem,hybrid,hvcoord,nets,nete,Tinit)
    case(1)
      call hs1_init_state(elem,hybrid,hvcoord,nets,nete,Tinit)
    case default
      call abortmp('invalid initialization sub_case: only subcase = 0, 1 supported')
    end select

  end subroutine hs_init_state

  subroutine hs0_init_state(elem, hybrid, hvcoord,nets,nete,Tinit)

    type(element_t),        intent(inout) :: elem(:)
    type(hybrid_t),         intent(in)    :: hybrid ! hybrid parallel structure
    type (hvcoord_t),       intent(in)    :: hvcoord
    integer,                intent(in)    :: nets
    integer,                intent(in)    :: nete
    real (kind=real_kind),  intent(in)    :: Tinit

    ! Local variables
    integer ie,i,j,k,q,tl
    integer :: nm1 
    integer :: n0 
    integer :: np1
    real (kind=real_kind) :: lat_mtn,lon_mtn,r_mtn,h_mtn,rsq,lat,lon
    real (kind=real_kind) :: temperature(np,np,nlev),p(np,np),exner(np,np),ps(np,np)

    nm1 = 1
    n0 = 2
    np1 = 3

    do ie = nets,nete
       elem(ie)%state%ps_v(:,:,n0) = hvcoord%ps0
       elem(ie)%state%ps_v(:,:,nm1) = hvcoord%ps0
       elem(ie)%state%ps_v(:,:,np1) = hvcoord%ps0

       elem(ie)%state%v(:,:,:,:,n0) = 0.0D0
       elem(ie)%state%v(:,:,:,:,nm1) = elem(ie)%state%v(:,:,:,:,n0)
       elem(ie)%state%v(:,:,:,:,np1) = elem(ie)%state%v(:,:,:,:,n0)

#ifdef MODEL_THETA_L
       elem(ie)%state%w_i = 0.0
#endif

       temperature(:,:,:) = Tinit

       ! if topo file was given in the namelist, PHIS was initilized in prim_main
       ! otherwise assume 0
#ifndef HOMME_WITHOUT_PIOLIBRARY
       if (infilenames(1)=='') then
          elem(ie)%state%phis(:,:) = 0.0D0
       end if
#endif

#undef HS_TOPO1
#ifdef HS_TOPO1
       lat_mtn = dd_pi/6
       lon_mtn = 3*dd_pi/2
       r_mtn = dd_pi/9
       h_mtn = 4000
       do i = 1,np
          do j = 1,np
             lat = elem(ie)%spherev(i,j)%lat
             lon = elem(ie)%spherev(i,j)%lon
             rsq = MIN((lat-lat_mtn)**2 + (lon - lon_mtn)**2,R_mtn**2)
             elem(ie)%state%phis(i,j) = g * h_mtn * (1.0D0 - SQRT(rsq) / R_mtn)
          end do
       end do
#endif

       ! initialize surface pressure to be consistent with topo
       elem(ie)%state%ps_v(:,:,n0) = elem(ie)%state%ps_v(:,:,n0) * &
            exp(-elem(ie)%state%phis(:,:) / (Rgas*Tinit))
       elem(ie)%state%ps_v(:,:,nm1) = elem(ie)%state%ps_v(:,:,n0)
       elem(ie)%state%ps_v(:,:,np1) = elem(ie)%state%ps_v(:,:,n0)

#if 0
       do k=1,nlev
          p(:,:) = hvcoord%hyam(k) * hvcoord%ps0 + hvcoord%hybm(k) * elem(ie)%state%ps_v(:,:,n0)
          exner(:,:) = (p(:,:)/hvcoord%ps0)**kappa
          temperature(:,:,k) = (Tinit - 150) + 150 * exner(:,:)
       end do
#endif

       if (qsize >= 1) then
          q=1
          elem(ie)%state%Q(:,:,:,q) = 0  ! moist HS tracer IC=0
          do q = 2,qsize
             elem(ie)%state%Q(:,:,:,q) = temperature(:,:,:) / 400
          end do
       end if
       ps = elem(ie)%state%ps_v(:,:,n0)
       call set_thermostate(elem(ie),ps,temperature,hvcoord)

    end do

  end subroutine hs0_init_state

  subroutine hs1_init_state(elem,hybrid,hvcoord,nets,nete,Tinit)

    type(element_t),        intent(inout) :: elem(:)
    type(hybrid_t),         intent(in)    :: hybrid  ! hybrid parallel structure
    type (hvcoord_t),       intent(in)    :: hvcoord
    integer,                intent(in)    :: nets
    integer,                intent(in)    :: nete
    real (kind=real_kind),  intent(in)    :: Tinit

    call hs0_init_state(elem,hybrid,hvcoord,nets,nete,Tinit)

  end subroutine hs1_init_state

end module held_suarez_mod


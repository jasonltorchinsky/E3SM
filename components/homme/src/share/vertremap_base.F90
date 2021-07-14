#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!SUBROUTINES:
!
!Notes on Lagrange+REMAP advection
!dynamics will compute mean fluxes, so that (i.e. for qsplit=3)
!
!    dp(t+3)-dp(t) = -3dt div(Udp_sum/3) - 3dt d(eta_dot_dpdn_sum/3)  + 3dt D(dpdiss_sum/3)
!
!Where the floating lagrangian component:
!    dp_star(t+3) = dp(t)  -3dt div(Udp_sum/3)  + 3dt D(dpdiss_sum/3)
!OR:
!    dp_star(t+3) = dp(t+1) + 3dt d( eta_dot_dpdn_ave(t) )


module vertremap_base

  !**************************************************************************************
  !
  !  Purpose:
  !        Construct sub-grid-scale polynomials using piecewise spline method with
  !        monotone filters.
  !
  !  References: PCM - Zerroukat et al., Q.J.R. Meteorol. Soc., 2005. (ZWS2005QJR)
  !              PSM - Zerroukat et al., Int. J. Numer. Meth. Fluids, 2005. (ZWS2005IJMF)
  !
  !**************************************************************************************

  use kinds, only                  : real_kind,int_kind
  use dimensions_mod, only         : np,nlev,qsize,nlevp,npsq
  use hybvcoord_mod, only          : hvcoord_t
  use perf_mod, only               : t_startf, t_stopf  ! _EXTERNAL
  use parallel_mod, only           : abortmp, parallel_t

! options:
! remap1                  ! remap any field, splines, monotone
! remap1_nofilter         ! remap any field, splines, no filter
! remap_q_ppm             ! remap state%Q, PPM, monotone

  contains


!=======================================================================================================!

!remap_calc_grids computes the vertical pressures and pressure differences for one vertical column for the reference grid
!and for the deformed Lagrangian grid. This was pulled out of each routine since it was a repeated task.
subroutine remap_calc_grids( hvcoord , ps , dt , eta_dot_dpdn , p_lag , p_ref , dp_lag , dp_ref )
  implicit none
  type(hvcoord_t)      , intent(in   ) :: hvcoord               !Derived type to hold vertical sigma grid parameters
  real(kind=real_kind) , intent(in   ) :: ps                    !Surface pressure for this column
  real(kind=real_kind) , intent(in   ) :: dt                    !Time step
  real(kind=real_kind) , intent(in   ) :: eta_dot_dpdn(nlev+1)  !Looks like a vertical pressure flux
                                                                !to compute deformed grid spacing
  real(kind=real_kind) , intent(  out) :: p_lag(nlev+1)         !Pressures at interfaces of the Lagrangian deformed grid
  real(kind=real_kind) , intent(  out) :: p_ref(nlev+1)         !Pressures at interfaces of the reference grid
  real(kind=real_kind) , intent(  out) :: dp_lag(nlev)          !Pressure differences on Lagrangian deformed grid
  real(kind=real_kind) , intent(  out) :: dp_ref(nlev)          !Pressure differences on reference grid
  integer :: k                                                  !Iterator
  p_ref(1) = 0  !Both grids have a model top pressure of zero
  p_lag(1) = 0  !Both grids have a model top pressure of zero
  do k = 1 , nlev
    dp_ref(k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) ) * hvcoord%ps0 + &
         ( hvcoord%hybi(k+1) - hvcoord%hybi(k) ) * ps  !Reference pressure difference
    ! Lagrangian pressure difference (flux in - flux out over the time step)
    dp_lag(k) = dp_ref(k) + dt * ( eta_dot_dpdn(k+1) - eta_dot_dpdn(k) )
    p_ref(k+1) = p_ref(k) + dp_ref(k) !Pressure at interfaces accumulated using difference over each cell
    p_lag(k+1) = p_lag(k) + dp_lag(k) !Pressure at interfaces accumulated using difference over each cell
  enddo
end subroutine remap_calc_grids

!=======================================================================================================!



subroutine remap1(Qdp,nx,qsize,dp1,dp2,remap_alg)
  ! remap 1 field
  ! input:  Qdp   field to be remapped (NOTE: MASS, not MIXING RATIO)
  !         dp1   layer thickness (source)
  !         dp2   layer thickness (target)
  !
  ! output: remaped Qdp, conserving mass, monotone on Q=Qdp/dp
  !
  implicit none
  integer, intent(in) :: nx,qsize,remap_alg
  real (kind=real_kind), intent(inout) :: Qdp(nx,nx,nlev,qsize)
  real (kind=real_kind), intent(in) :: dp1(nx,nx,nlev),dp2(nx,nx,nlev)
  ! ========================
  ! Local Variables
  ! ========================

  real (kind=real_kind), dimension(nlev+1)    :: rhs,lower_diag,diag,upper_diag,q_diag,zgam,z1c,z2c,zv
  real (kind=real_kind), dimension(nlev)      :: h,Qcol,dy,za0,za1,za2,zarg,zhdp,dp_star,dp_np1
  real (kind=real_kind)  :: f_xm,level1,level2,level3,level4,level5, &
                            peaks_min,peaks_max,tmp_cal,xm,xm_d,zv1,zv2, &
                            zero = 0,one = 1,tiny = 1e-12,qmax = 1d50
  integer(kind=int_kind) :: zkr(nlev+1),filter_code(nlev),peaks,im1,im2,im3,ip1,ip2, &
                            lt1,lt2,lt3,t0,t1,t2,t3,t4,tm,tp,ie,i,ilev,j,jk,k,q
  logical :: abrtf=.false.

  q = remap_alg
  if (         (q .ne. -1) .and. (q .ne.  0) .and. (q .ne.  1) &
       & .and. (q .ne. 10) .and. (q .ne. 11) .and. (q .ne. 12) ) then
     call abortmp('Bad remap_alg value. Use -1, 0, 1, 10, 11, or 12.')
  end if

  if (remap_alg == -1) then
     call remap1_nofilter(qdp,nx,qsize,dp1,dp2)
     return
  endif
  if (remap_alg >= 1) then
     call remap_Q_ppm(qdp,nx,qsize,dp1,dp2,remap_alg)
     return
  endif

  call t_startf('remap_Q_noppm')
#if (defined COLUMN_OPENMP)
!$omp parallel do private(q,i,j,z1c,z2c,zv,k,dp_np1,dp_star,Qcol,zkr,ilev) &
!$omp    private(jk,zgam,zhdp,h,zarg,rhs,lower_diag,diag,upper_diag,q_diag,tmp_cal,filter_code) &
!$omp    private(dy,im1,im2,im3,ip1,t1,t2,t3,za0,za1,za2,xm_d,xm,f_xm,t4,tm,tp,peaks,peaks_min) &
!$omp    private(peaks_max,ip2,level1,level2,level3,level4,level5,lt1,lt2,lt3,zv1,zv2)
#endif
  do q=1,qsize
  do i=1,nx
    do j=1,nx

      z1c(1)=0 ! source grid
      z2c(1)=0 ! target grid
      do k=1,nlev
         z1c(k+1)=z1c(k)+dp1(i,j,k)
         z2c(k+1)=z2c(k)+dp2(i,j,k)
      enddo

      zv(1)=0
      do k=1,nlev
        Qcol(k)=Qdp(i,j,k,q)!  *(z1c(k+1)-z1c(k)) input is mass
        zv(k+1) = zv(k)+Qcol(k)
      enddo

      if (ABS(z2c(nlev+1)-z1c(nlev+1)).GE.0.000001) then
        write(6,*) 'SURFACE PRESSURE IMPLIED BY ADVECTION SCHEME'
        write(6,*) 'NOT CORRESPONDING TO SURFACE PRESSURE IN    '
        write(6,*) 'DATA FOR MODEL LEVELS'
        write(6,*) 'PLEVMODEL=',z2c(nlev+1)
        write(6,*) 'PLEV     =',z1c(nlev+1)
        write(6,*) 'DIFF     =',z2c(nlev+1)-z1c(nlev+1)
        abrtf=.true.
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! quadratic splines with UK met office monotonicity constraints !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      zkr  = 0 !99
      ilev = 2
      zkr(1) = 1
      zkr(nlev+1) = nlev
      kloop: do k = 2,nlev
        do jk = ilev,nlev+1
          if (z1c(jk).ge.z2c(k)) then
            ilev      = jk
            zkr(k)   = jk-1
            cycle kloop
          endif
        enddo
      enddo kloop

      !if(any(zkr==0)) then
      !  write(6,*) "Invalid zkr values in remap1."
      !  abrtf=.true.
      !endif

      zgam  = (z2c(1:nlev+1)-z1c(zkr)) / (z1c(zkr+1)-z1c(zkr))
      zgam(1)      = 0.0
      zgam(nlev+1) = 1.0
      zhdp = z1c(2:nlev+1)-z1c(1:nlev)


      h = 1/zhdp
      zarg = Qcol * h
      rhs = 0
      lower_diag = 0
      diag = 0
      upper_diag = 0

      rhs(1)=3*zarg(1)
      rhs(2:nlev) = 3*(zarg(2:nlev)*h(2:nlev) + zarg(1:nlev-1)*h(1:nlev-1))
      rhs(nlev+1)=3*zarg(nlev)

      lower_diag(1)=1
      lower_diag(2:nlev) = h(1:nlev-1)
      lower_diag(nlev+1)=1

      diag(1)=2
      diag(2:nlev) = 2*(h(2:nlev) + h(1:nlev-1))
      diag(nlev+1)=2

      upper_diag(1)=1
      upper_diag(2:nlev) = h(2:nlev)
      upper_diag(nlev+1)=0

      q_diag(1)=-upper_diag(1)/diag(1)
      rhs(1)= rhs(1)/diag(1)

      do k=2,nlev+1
        tmp_cal    =  1/(diag(k)+lower_diag(k)*q_diag(k-1))
        q_diag(k) = -upper_diag(k)*tmp_cal
        rhs(k) =  (rhs(k)-lower_diag(k)*rhs(k-1))*tmp_cal
      enddo
      do k=nlev,1,-1
        rhs(k)=rhs(k)+q_diag(k)*rhs(k+1)
      enddo

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!  monotonicity modifications  !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      filter_code = 0
      dy(1:nlev-1) = zarg(2:nlev)-zarg(1:nlev-1)
      dy(nlev) = dy(nlev-1)

      dy = merge(zero, dy, abs(dy) < tiny )

      do k=1,nlev
        im1=MAX(1,k-1)
        im2=MAX(1,k-2)
        im3=MAX(1,k-3)
        ip1=MIN(nlev,k+1)
        t1 = merge(1,0,(zarg(k)-rhs(k))*(rhs(k)-zarg(im1)) >= 0)
        t2 = merge(1,0,dy(im2)*(rhs(k)-zarg(im1)) > 0 .AND. dy(im2)*dy(im3) > 0 &
             .AND. dy(k)*dy(ip1) > 0 .AND. dy(im2)*dy(k) < 0 )
        t3 = merge(1,0,ABS(rhs(k)-zarg(im1)) > ABS(rhs(k)-zarg(k)))

        filter_code(k) = merge(0,1,t1+t2 > 0)
        rhs(k) = (1-filter_code(k))*rhs(k)+filter_code(k)*(t3*zarg(k)+(1-t3)*zarg(im1))
        filter_code(im1) = MAX(filter_code(im1),filter_code(k))
      enddo

      rhs = merge(qmax,rhs,rhs > qmax)
      rhs = merge(zero,rhs,rhs < zero)

      za0 = rhs(1:nlev)
      za1 = -4*rhs(1:nlev) - 2*rhs(2:nlev+1) + 6*zarg
      za2 =  3*rhs(1:nlev) + 3*rhs(2:nlev+1) - 6*zarg

      dy(1:nlev) = rhs(2:nlev+1)-rhs(1:nlev)
      dy = merge(zero, dy, abs(dy) < tiny )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Compute the 3 quadratic spline coeffients {za0, za1, za2}				   !!
      !! knowing the quadratic spline parameters {rho_left,rho_right,zarg}		   !!
      !! Zerroukat et.al., Q.J.R. Meteorol. Soc., Vol. 128, pp. 2801-2820 (2002).   !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      h = rhs(2:nlev+1)

      do k=1,nlev
        xm_d = merge(one,2*za2(k),abs(za2(k)) < tiny)
        xm = merge(zero,-za1(k)/xm_d, abs(za2(k)) < tiny)
        f_xm = za0(k) + za1(k)*xm + za2(k)*xm**2

        t1 = merge(1,0,ABS(za2(k)) > tiny)
        t2 = merge(1,0,xm <= zero .OR. xm >= 1)
        t3 = merge(1,0,za2(k) > zero)
        t4 = merge(1,0,za2(k) < zero)
        tm = merge(1,0,t1*((1-t2)+t3) .EQ. 2)
        tp = merge(1,0,t1*((1-t2)+(1-t3)+t4) .EQ. 3)

        peaks=0
        peaks = merge(-1,peaks,tm .EQ. 1)
        peaks = merge(+1,peaks,tp .EQ. 1)
        peaks_min = merge(f_xm,MIN(za0(k),za0(k)+za1(k)+za2(k)),tm .EQ. 1)
        peaks_max = merge(f_xm,MAX(za0(k),za0(k)+za1(k)+za2(k)),tp .EQ. 1)

        im1=MAX(1,k-1)
        im2=MAX(1,k-2)
        ip1=MIN(nlev,k+1)
        ip2=MIN(nlev,k+2)

        t1 = merge(abs(peaks),0,(dy(im2)*dy(im1) <= tiny) .OR. &
             (dy(ip1)*dy(ip2) <= tiny) .OR. (dy(im1)*dy(ip1) >= tiny) .OR. &
             (dy(im1)*float(peaks) <= tiny))

        filter_code(k) = merge(1,t1+(1-t1)*filter_code(k),(rhs(k) >= qmax) .OR. &
             (rhs(k) <= zero) .OR. (peaks_max > qmax) .OR. (peaks_min < tiny))

        if (filter_code(k) > 0) then
          level1 = rhs(k)
          level2 = (2*rhs(k)+h(k))/3
          level3 = 0.5*(rhs(k)+h(k))
          level4 = (1/3d0)*rhs(k)+2*(1/3d0)*h(k)
          level5 = h(k)

          t1 = merge(1,0,h(k) >= rhs(k))
          t2 = merge(1,0,zarg(k) <= level1 .OR.  zarg(k) >= level5)
          t3 = merge(1,0,zarg(k) >  level1 .AND. zarg(k) <  level2)
          t4 = merge(1,0,zarg(k) >  level4 .AND. zarg(k) <  level5)

          lt1 = t1*t2
          lt2 = t1*(1-t2+t3)
          lt3 = t1*(1-t2+1-t3+t4)

          za0(k) = merge(zarg(k),za0(k),lt1 .EQ. 1)
          za1(k) = merge(zero,za1(k),lt1 .EQ. 1)
          za2(k) = merge(zero,za2(k),lt1 .EQ. 1)

          za0(k) = merge(rhs(k),za0(k),lt2 .EQ. 2)
          za1(k) = merge(zero,za1(k),lt2 .EQ. 2)
          za2(k) = merge(3*(zarg(k)-rhs(k)),za2(k),lt2 .EQ. 2)

          za0(k) = merge(-2*h(k)+3*zarg(k),za0(k),lt3 .EQ. 3)
          za1(k) = merge(+6*h(k)-6*zarg(k),za1(k),lt3 .EQ. 3)
          za2(k) = merge(-3*h(k)+3*zarg(k),za2(k),lt3 .EQ. 3)

          t2 = merge(1,0,zarg(k) >= level1 .OR.  zarg(k) <= level5)
          t3 = merge(1,0,zarg(k) <  level1 .AND. zarg(k) >  level2)
          t4 = merge(1,0,zarg(k) <  level4 .AND. zarg(k) >  level5)

          lt1 = (1-t1)*t2
          lt2 = (1-t1)*(1-t2+t3)
          lt3 = (1-t1)*(1-t2+1-t3+t4)

          za0(k) = merge(zarg(k),za0(k),lt1 .EQ. 1)
          za1(k) = merge(zero,za1(k),lt1 .EQ. 1)
          za2(k) = merge(zero,za2(k),lt1 .EQ. 1)

          za0(k) = merge(rhs(k),za0(k),lt2 .EQ. 2)
          za1(k) = merge(zero,za1(k),lt2 .EQ. 2)
          za2(k) = merge(3*(zarg(k)-rhs(k)),za2(k),lt2 .EQ. 2)

          za0(k) = merge(-2*h(k)+3*zarg(k),za0(k),lt3 .EQ. 3)
          za1(k) = merge(+6*h(k)-6*zarg(k),za1(k),lt3 .EQ. 3)
          za2(k) = merge(-3*h(k)+3*zarg(k),za2(k),lt3 .EQ. 3)
        endif
      enddo

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! start iteration from top to bottom of atmosphere !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      zv1 = 0
      do k=1,nlev
        if (zgam(k+1)>1d0) then
          WRITE(*,*) 'r not in [0:1]', zgam(k+1)
          abrtf=.true.
        endif
        zv2 = zv(zkr(k+1))+(za0(zkr(k+1))*zgam(k+1)+(za1(zkr(k+1))/2)*(zgam(k+1)**2)+ &
             (za2(zkr(k+1))/3)*(zgam(k+1)**3))*zhdp(zkr(k+1))
        Qdp(i,j,k,q) = (zv2 - zv1) ! / (z2c(k+1)-z2c(k) ) dont convert back to mixing ratio
        zv1 = zv2
      enddo
    enddo
  enddo
  enddo ! q loop
  if (abrtf) call abortmp('Bad levels in remap1.  usually CFL violatioin')
  call t_stopf('remap_Q_noppm')

end subroutine remap1

subroutine remap1_nofilter(Qdp,nx,qsize,dp1,dp2)
  ! remap 1 field
  ! input:  Qdp   field to be remapped (NOTE: MASS, not MIXING RATIO)
  !         dp1   layer thickness (source)
  !         dp2   layer thickness (target)
  !
  ! output: remaped Qdp, conserving mass
  !
  implicit none
  integer, intent(in) :: nx,qsize
  real (kind=real_kind), intent(inout) :: Qdp(nx,nx,nlev,qsize)
  real (kind=real_kind), intent(in) :: dp1(nx,nx,nlev),dp2(nx,nx,nlev)
  ! ========================
  ! Local Variables
  ! ========================

  real (kind=real_kind), dimension(nlev+1)    :: rhs,lower_diag,diag,upper_diag,q_diag,zgam,z1c,z2c,zv
  real (kind=real_kind), dimension(nlev)      :: h,Qcol,dy,za0,za1,za2,zarg,zhdp,dp_star,dp_np1
  real (kind=real_kind)  :: f_xm,level1,level2,level3,level4,level5, &
                            peaks_min,peaks_max,tmp_cal,xm,xm_d,zv1,zv2, &
                            zero = 0,one = 1,tiny = 1e-12,qmax = 1d50
  integer(kind=int_kind) :: zkr(nlev+1),filter_code(nlev),peaks,im1,im2,im3,ip1,ip2, &
                            lt1,lt2,lt3,t0,t1,t2,t3,t4,tm,tp,ie,i,ilev,j,jk,k,q
  logical :: abrtf=.false.
!   call t_startf('remap1_nofilter')

#if (defined COLUMN_OPENMP)
!$omp parallel do private(q,i,j,z1c,z2c,zv,k,dp_np1,dp_star,Qcol,zkr,ilev) &
!$omp    private(jk,zgam,zhdp,h,zarg,rhs,lower_diag,diag,upper_diag,q_diag,tmp_cal,filter_code) &
!$omp    private(dy,im1,im2,im3,ip1,t1,t2,t3,za0,za1,za2,xm_d,xm,f_xm,t4,tm,tp,peaks,peaks_min) &
!$omp    private(peaks_max,ip2,level1,level2,level3,level4,level5,lt1,lt2,lt3,zv1,zv2)
#endif
  do q=1,qsize
  do i=1,nx
    do j=1,nx

      z1c(1)=0 ! source grid
      z2c(1)=0 ! target grid
      do k=1,nlev
         z1c(k+1)=z1c(k)+dp1(i,j,k)
         z2c(k+1)=z2c(k)+dp2(i,j,k)
      enddo

      zv(1)=0
      do k=1,nlev
        Qcol(k)=Qdp(i,j,k,q)!  *(z1c(k+1)-z1c(k)) input is mass
        zv(k+1) = zv(k)+Qcol(k)
      enddo

      if (ABS(z2c(nlev+1)-z1c(nlev+1)).GE.0.000001) then
        write(6,*) 'SURFACE PRESSURE IMPLIED BY ADVECTION SCHEME'
        write(6,*) 'NOT CORRESPONDING TO SURFACE PRESSURE IN    '
        write(6,*) 'DATA FOR MODEL LEVELS'
        write(6,*) 'PLEVMODEL=',z2c(nlev+1)
        write(6,*) 'PLEV     =',z1c(nlev+1)
        write(6,*) 'DIFF     =',z2c(nlev+1)-z1c(nlev+1)
        abrtf=.true.
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! quadratic splies with UK met office monotonicity constraints  !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      zkr  = 99
      ilev = 2
      zkr(1) = 1
      zkr(nlev+1) = nlev
      kloop: do k = 2,nlev
        do jk = ilev,nlev+1
          if (z1c(jk).ge.z2c(k)) then
            ilev      = jk
            zkr(k)   = jk-1
            cycle kloop
          endif
        enddo
      enddo kloop

      zgam  = (z2c(1:nlev+1)-z1c(zkr)) / (z1c(zkr+1)-z1c(zkr))
      zgam(1)      = 0.0
      zgam(nlev+1) = 1.0
      zhdp = z1c(2:nlev+1)-z1c(1:nlev)


      h = 1/zhdp
      zarg = Qcol * h
      rhs = 0
      lower_diag = 0
      diag = 0
      upper_diag = 0

      rhs(1)=3*zarg(1)
      rhs(2:nlev) = 3*(zarg(2:nlev)*h(2:nlev) + zarg(1:nlev-1)*h(1:nlev-1))
      rhs(nlev+1)=3*zarg(nlev)

      lower_diag(1)=1
      lower_diag(2:nlev) = h(1:nlev-1)
      lower_diag(nlev+1)=1

      diag(1)=2
      diag(2:nlev) = 2*(h(2:nlev) + h(1:nlev-1))
      diag(nlev+1)=2

      upper_diag(1)=1
      upper_diag(2:nlev) = h(2:nlev)
      upper_diag(nlev+1)=0

      q_diag(1)=-upper_diag(1)/diag(1)
      rhs(1)= rhs(1)/diag(1)

      do k=2,nlev+1
        tmp_cal    =  1/(diag(k)+lower_diag(k)*q_diag(k-1))
        q_diag(k) = -upper_diag(k)*tmp_cal
        rhs(k) =  (rhs(k)-lower_diag(k)*rhs(k-1))*tmp_cal
      enddo
      do k=nlev,1,-1
        rhs(k)=rhs(k)+q_diag(k)*rhs(k+1)
      enddo

      za0 = rhs(1:nlev)
      za1 = -4*rhs(1:nlev) - 2*rhs(2:nlev+1) + 6*zarg
      za2 =  3*rhs(1:nlev) + 3*rhs(2:nlev+1) - 6*zarg


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! start iteration from top to bottom of atmosphere !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      zv1 = 0
      do k=1,nlev
        if (zgam(k+1)>1d0) then
          WRITE(*,*) 'r not in [0:1]', zgam(k+1)
          abrtf=.true.
        endif
        zv2 = zv(zkr(k+1))+(za0(zkr(k+1))*zgam(k+1)+(za1(zkr(k+1))/2)*(zgam(k+1)**2)+ &
             (za2(zkr(k+1))/3)*(zgam(k+1)**3))*zhdp(zkr(k+1))
        Qdp(i,j,k,q) = (zv2 - zv1) ! / (z2c(k+1)-z2c(k) ) dont convert back to mixing ratio
        zv1 = zv2
      enddo
    enddo
  enddo
  enddo ! q loop
  if (abrtf) call abortmp('Bad levels in remap1_nofilter.  usually CFL violatioin')
!   call t_stopf('remap1_nofilter')
end subroutine remap1_nofilter

!=======================================================================================================!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This uses the exact same model and reference grids and data as remap_Q, but it interpolates
!using PPM instead of splines.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! remap_alg == 1  means mirrored values in ghost cells (i.e., no flux)
!! remap_alg == 2  means 1st order reconstruction in ghost cells
!! remap_alg == 10 linear extrapolation with global bounds preservation
!! remap_alg == 11 linear extrapolation 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine remap_Q_ppm(Qdp,nx,qsize,dp1,dp2,remap_alg)
  ! remap 1 field
  ! input:  Qdp   field to be remapped (NOTE: MASS, not MIXING RATIO)
  !         dp1   layer thickness (source)
  !         dp2   layer thickness (target)
  !
  ! output: remaped Qdp, conserving mass
  !
  implicit none
  integer,intent(in) :: nx,qsize, remap_alg
  real (kind=real_kind), intent(inout) :: Qdp(nx,nx,nlev,qsize)
  real (kind=real_kind), intent(in) :: dp1(nx,nx,nlev),dp2(nx,nx,nlev)
  ! Local Variables
  integer, parameter :: gs = 2                              !Number of cells to place in the ghost region
  real(kind=real_kind), dimension(       nlev+2 ) :: pio    !Pressure at interfaces for old grid
  real(kind=real_kind), dimension(       nlev+1 ) :: pin    !Pressure at interfaces for new grid
  real(kind=real_kind), dimension(       nlev+1 ) :: masso  !Accumulate mass up to each interface
  real(kind=real_kind), dimension(  1-gs:nlev+gs) :: ao     !Tracer value on old grid
  real(kind=real_kind), dimension(  1-gs:nlev+gs) :: dpo    !change in pressure over a cell for old grid
  real(kind=real_kind), dimension(  1-gs:nlev+gs) :: dpn    !change in pressure over a cell for old grid
  real(kind=real_kind), dimension(3,     nlev   ) :: coefs  !PPM coefficients within each cell
  real(kind=real_kind), dimension(       nlev   ) :: z1, z2
  real(kind=real_kind) :: ppmdx(10,0:nlev+1)  !grid spacings
  real(kind=real_kind) :: mymass, massn1, massn2, ext(2)
  integer :: i, j, k, q, kk, kid(nlev)

  call t_startf('remap_Q_ppm')
  do j = 1 , nx
    do i = 1 , nx

      pin(1)=0
      pio(1)=0
      do k=1,nlev
         dpn(k)=dp2(i,j,k)
         dpo(k)=dp1(i,j,k)
         pin(k+1)=pin(k)+dpn(k)
         pio(k+1)=pio(k)+dpo(k)
      enddo



      pio(nlev+2) = pio(nlev+1) + 1.  !This is here to allow an entire block of k threads to run in the remapping phase.
                                      !It makes sure there's an old interface value below the domain that is larger.
      pin(nlev+1) = pio(nlev+1)       !The total mass in a column does not change.
                                      !Therefore, the pressure of that mass cannot either.
      !Fill in the ghost regions with mirrored values. if remap_alg is defined, this is of no consequence.
      do k = 1 , gs
        dpo(1   -k) = dpo(       k)
        dpo(nlev+k) = dpo(nlev+1-k)
      enddo

      !Compute remapping intervals once for all tracers. Find the old grid cell index in which the
      !k-th new cell interface resides. Then integrate from the bottom of that old cell to the new
      !interface location. In practice, the grid never deforms past one cell, so the search can be
      !simplified by this. Also, the interval of integration is usually of magnitude close to zero
      !or close to dpo because of minimial deformation.
      !Numerous tests confirmed that the bottom and top of the grids match to machine precision, so
      !I set them equal to each other.
      do k = 1 , nlev
        kk = k  !Keep from an order n^2 search operation by assuming the old cell index is close.
        !Find the index of the old grid cell in which this new cell's bottom interface resides.
        if (pio(kk) <= pin(k+1)) then
           do while ( pio(kk) <= pin(k+1) )
              kk = kk + 1
           enddo
           kk = kk - 1                   !kk is now the cell index we're integrating over.
        else
           call binary_search(pio, pin(k+1), kk)
        end if
        if (kk == nlev+1) kk = nlev   !This is to keep the indices in bounds.
                                      !Top bounds match anyway, so doesn't matter what coefficients are used
        kid(k) = kk                   !Save for reuse
        z1(k) = -0.5D0                !This remapping assumes we're starting from the left interface of an old grid cell
                                      !In fact, we're usually integrating very little or almost all of the cell in question
        z2(k) = ( pin(k+1) - ( pio(kk) + pio(kk+1) ) * 0.5 ) / dpo(kk)  !PPM interpolants are normalized to an independent
                                                                        !coordinate domain [-0.5,0.5].
      enddo

      !This turned out a big optimization, remembering that only parts of the PPM algorithm depends on the data, namely the
      !limiting. So anything that depends only on the grid is pre-computed outside the tracer loop.
      ppmdx(:,:) = compute_ppm_grids( dpo )

      !From here, we loop over tracers for only those portions which depend on tracer data, which includes PPM limiting and
      !mass accumulation

      do q = 1 , qsize
        !Accumulate the old mass up to old grid cell interface locations to simplify integration
        !during remapping. Also, divide out the grid spacing so we're working with actual tracer
        !values and can conserve mass. The option for ifndef ZEROHORZ I believe is there to ensure
        !tracer consistency for an initially uniform field. I copied it from the old remap routine.
        masso(1) = 0.
        do k = 1 , nlev
          ao(k) = Qdp(i,j,k,q)
          masso(k+1) = masso(k) + ao(k) !Accumulate the old mass. This will simplify the remapping
          ao(k) = ao(k) / dpo(k)        !Divide out the old grid spacing because we want the tracer mixing ratio, not mass.
        enddo
        !Fill in ghost values.
        if (remap_alg == 10) then
           ext(1) = minval(ao(1:nlev))
           ext(2) = maxval(ao(1:nlev))
           call linextrap(dpo(2), dpo(1), dpo(0), dpo(-1), ao(2), ao(1), ao(0), ao(-1), 1,ext(1), ext(2))
           call linextrap(dpo(nlev-1), dpo(nlev), dpo(nlev+1), dpo(nlev+2),&
                ao(nlev-1), ao(nlev), ao(nlev+1), ao(nlev+2), 1, ext(1), ext(2))
        else if (remap_alg == 11 .or. remap_alg == 12) then
           call linextrap(dpo(2), dpo(1), dpo(0), dpo(-1), ao(2), ao(1), ao(0), ao(-1), 0,ext(1), ext(2))
           call linextrap(dpo(nlev-1), dpo(nlev), dpo(nlev+1), dpo(nlev+2),&
                ao(nlev-1), ao(nlev), ao(nlev+1), ao(nlev+2), 0, ext(1), ext(2))
        else
           !values for gs>1 are not used
           !For remap_alg == 2 values might be used in computations that are then overwritten,
           ! so they should be set but wont effect answer
           do k = 1 , gs
              ao(1   -k) = ao(1)
              ao(nlev+k) = ao(nlev)
           enddo
        endif
        !Compute monotonic and conservative PPM reconstruction over every cell
        coefs(:,:) = compute_ppm( ao , ppmdx, remap_alg )
        !Compute tracer values on the new grid by integrating from the old cell bottom to the new
        !cell interface to form a new grid mass accumulation. Taking the difference between
        !accumulation at successive interfaces gives the mass inside each cell. Since Qdp is
        !supposed to hold the full mass this needs no normalization.
        massn1 = 0.
        do k = 1 , nlev
          kk = kid(k)
          massn2 = masso(kk) + integrate_parabola( coefs(:,kk) , z1(k) , z2(k) ) * dpo(kk)
          Qdp(i,j,k,q) = massn2 - massn1
          massn1 = massn2
        enddo
      enddo
    enddo
  enddo
  call t_stopf('remap_Q_ppm')
end subroutine remap_Q_ppm


!THis compute grid-based coefficients from Collela & Woodward 1984.
function compute_ppm_grids( dx )   result(rslt)
  implicit none
  real(kind=real_kind), intent(in) :: dx(-1:nlev+2)  !grid spacings
  real(kind=real_kind)             :: rslt(10,0:nlev+1)  !grid spacings
  integer :: j
  integer :: indB, indE

  !Calculate grid-based coefficients for stage 1 of compute_ppm
  do j = 0 , nlev+1
    rslt( 1,j) = dx(j) / ( dx(j-1) + dx(j) + dx(j+1) )
    rslt( 2,j) = ( 2.*dx(j-1) + dx(j) ) / ( dx(j+1) + dx(j) )
    rslt( 3,j) = ( dx(j) + 2.*dx(j+1) ) / ( dx(j-1) + dx(j) )
  enddo

  !Caculate grid-based coefficients for stage 2 of compute_ppm
  do j = 0 , nlev
    rslt( 4,j) = dx(j) / ( dx(j) + dx(j+1) )
    rslt( 5,j) = 1. / sum( dx(j-1:j+2) )
    rslt( 6,j) = ( 2. * dx(j+1) * dx(j) ) / ( dx(j) + dx(j+1 ) )
    rslt( 7,j) = ( dx(j-1) + dx(j  ) ) / ( 2. * dx(j  ) + dx(j+1) )
    rslt( 8,j) = ( dx(j+2) + dx(j+1) ) / ( 2. * dx(j+1) + dx(j  ) )
    rslt( 9,j) = dx(j  ) * ( dx(j-1) + dx(j  ) ) / ( 2.*dx(j  ) +    dx(j+1) )
    rslt(10,j) = dx(j+1) * ( dx(j+1) + dx(j+2) ) / (    dx(j  ) + 2.*dx(j+1) )
  enddo
end function compute_ppm_grids

!=======================================================================================================!



!This computes a limited parabolic interpolant using a net 5-cell stencil, but the stages of computation are broken up into 3 stages
function compute_ppm( a , dx, remap_alg )    result(coefs)
  implicit none
  real(kind=real_kind), intent(in) :: a    (    -1:nlev+2)  !Cell-mean values
  real(kind=real_kind), intent(in) :: dx   (10,  0:nlev+1)  !grid spacings
  integer :: remap_alg
  real(kind=real_kind) ::             coefs(0:2,   nlev  )  !PPM coefficients (for parabola)
  real(kind=real_kind) :: ai (0:nlev  )                     !fourth-order accurate, then limited interface values
  real(kind=real_kind) :: dma(0:nlev+1)                     !An expression from Collela's '84 publication
  real(kind=real_kind) :: da                                !Ditto
  ! Hold expressions based on the grid (which are cumbersome).
  real(kind=real_kind) :: dx1, dx2, dx3, dx4, dx5, dx6, dx7, dx8, dx9, dx10
  real(kind=real_kind) :: al, ar                            !Left and right interface values for cell-local limiting
  integer :: j
  integer :: indB, indE

  ! for vert_remap_q_alg=2, use 1st order reconstruction in
  ! the specified number of layers near model top and surface
  integer :: vert_remap_tom = 1   
  integer :: vert_remap_bl = 1


  ! Stage 1: Compute dma for each cell, allowing a 1-cell ghost stencil below and above the domain
  do j = 0 , nlev+1
    da = dx(1,j) * ( dx(2,j) * ( a(j+1) - a(j) ) + dx(3,j) * ( a(j) - a(j-1) ) )
    dma(j) = minval( (/ abs(da) , 2. * abs( a(j) - a(j-1) ) , 2. * abs( a(j+1) - a(j) ) /) ) * sign(1.D0,da)
    if ( ( a(j+1) - a(j) ) * ( a(j) - a(j-1) ) <= 0. ) dma(j) = 0.
  enddo

  ! Stage 2: Compute ai for each cell interface in the physical domain (dimension nlev+1)
  do j = 0 , nlev
    ai(j) = a(j) + dx(4,j) * ( a(j+1) - a(j) ) + dx(5,j) * ( dx(6,j) * ( dx(7,j) - dx(8,j) ) &
         * ( a(j+1) - a(j) ) - dx(9,j) * dma(j+1) + dx(10,j) * dma(j) )
  enddo

  ! Stage 3: Compute limited PPM interpolant over each cell in the physical domain
  ! (dimension nlev) using ai on either side and ao within the cell.
  do j = 1 , nlev
    al = ai(j-1)
    ar = ai(j  )
    if ( (ar - a(j)) * (a(j) - al) <= 0. ) then
      al = a(j)
      ar = a(j)
    endif
    if ( (ar - al) * (a(j) - (al + ar)/2.) >  (ar - al)**2/6. ) al = 3.*a(j) - 2. * ar
    if ( (ar - al) * (a(j) - (al + ar)/2.) < -(ar - al)**2/6. ) ar = 3.*a(j) - 2. * al
    !Computed these coefficients from the edge values and cell mean in Maple. Assumes normalized coordinates: xi=(x-x0)/dx
    coefs(0,j) = 1.5 * a(j) - ( al + ar ) / 4.
    coefs(1,j) = ar - al
    ! coefs(2,j) = -6. * a(j) + 3. * ( al + ar )
    coefs(2,j) = 3. * (-2. * a(j) + ( al + ar ))
  enddo

  ! Switch to piecewise constant near the boundaries
  if (remap_alg == 2 .or. remap_alg == 12) then
     do j=1,vert_remap_tom
        coefs(0,j) = a(j)
        coefs(1:2,j) = 0.D0
     enddo
     do j=nlev+1-vert_remap_bl,nlev
        coefs(0,j) = a(j)
        coefs(1:2,j) = 0.D0
     enddo
  endif

end function compute_ppm

!=======================================================================================================!


!Simple function computes the definite integral of a parabola in normalized coordinates, xi=(x-x0)/dx,
!given two bounds. Make sure this gets inlined during compilation.
function integrate_parabola( a , x1 , x2 )    result(mass)
  implicit none
  real(kind=real_kind), intent(in) :: a(0:2)  !Coefficients of the parabola
  real(kind=real_kind), intent(in) :: x1      !lower domain bound for integration
  real(kind=real_kind), intent(in) :: x2      !upper domain bound for integration
  real(kind=real_kind)             :: mass
  mass = a(0) * (x2 - x1) + a(1) * (x2 * x2 - x1 * x1) / 0.2D1 + a(2) * (x2 * x2 * x2 - x1 * x1 * x1) / 0.3D1
end function integrate_parabola


!=============================================================================================!

  ! Find k such that pio(k) <= pivot < pio(k+1). Provide a reasonable input
  ! value for k.
  subroutine binary_search(pio, pivot, k)
    real(kind=real_kind), intent(in) :: pio(nlev+2), pivot
    integer, intent(inout) :: k
    integer :: lo, hi, mid

    if (pio(k) > pivot) then
       lo = 1
       hi = k
    else
       lo = k
       hi = nlev+2
    end if
    do while (hi > lo + 1)
       k = (lo + hi)/2
       if (pio(k) > pivot) then
          hi = k
       else
          lo = k
       end if
    end do
    k = lo
  end subroutine binary_search


  subroutine linextrap(dx1,dx2,dx3,dx4,y1,y2,y3,y4,limit,lo,hi)
    real(kind=real_kind), intent(in) :: dx1,dx2,dx3,dx4,y1,y2,lo,hi
    real(kind=real_kind), intent(out) :: y3,y4

    real(kind=real_kind) :: den,num,a
    real(kind=real_kind) :: z3,z4
    integer :: limit

    ! In exact arithmetic, the following is equivalent to
    !   x1 = half*dx1
    !   x2 = x1 + half*(dx1 + dx2)
    !   x3 = x2 + half*(dx2 + dx3)
    !   x4 = x3 + half*(dx3 + dx4)
    !   a  = (x3-x1)/(x2-x1)
    !   y3 = (1-a)*y1 + a*y2
    !   a  = (x4-x1)/(x2-x1)
    !   y4 = (1-a)*y1 + a*y2
    ! In FP, -fp-model-fast -O3 produces different results depending
    ! on whether -qopenmp is present.

    den = (dx1 + dx2)/2
    num = den + (dx2 + dx3)/2
    a = num/den
    y3 = (1-a)*y1 + a*y2

    num = num + (dx3 + dx4)/2
    a  = num/den
    y4 = (1-a)*y1 + a*y2

    if (limit==1) then
       y3 = max(lo, min(hi, y3))
       y4 = max(lo, min(hi, y4))
    endif
  end subroutine linextrap

!=============================================================================
! "Corrects" remapped solution to fall within the bounds specified via 
! mass-borrowing. Essentially, the bounds specified should be the global
! bounds of the original solution, and we shift mass from either side of 
! the domain toward the interior (bot->top first, top-bot second).
!
! Qdp - Remapped field to be adjusted to obey bounds
!       (NOTE: MASS not MIXING RATIO).
! dp2 - Layer thickness (target).
! bnds - (Min,Max) for each field q (NOTE: MIXING RATIO not MASS).
!
! Output - Qdp now obeying global bounds.
!
! TO-DO: Increase performance, parallelize with OpenMP.
!=============================================================================

  subroutine mass_borrow(Qdp, nx, nlev, qsize, dp2, bnds)
    
    implicit none

    integer (kind=int_kind), intent(in)  :: nx, nlev, qsize
    real (kind=real_kind), intent(inout) :: Qdp(nx, nx, nlev, qsize)
    real (kind=real_kind), intent(in)    :: dp2(nx, nx, nlev)
    real (kind=real_kind), intent(in)    :: bnds(2, nx, nx, nlev, qsize) ! Min/Max for each q
    
    !================
    ! Local variables
    !================
    real (kind=real_kind)   :: mr ! Mixing ratio
    real (kind=real_kind)   :: lbnd, ubnd ! Lower, uppoer bounds for column
    real (kind=real_kind)   :: mr_diff ! Difference between mixing ratio and
    ! bound of interest.
    real (kind=real_kind)   :: ms_diff ! Mass associated with mr
    integer (kind=int_kind) :: ii, jj, kk, ll ! Iterators.
    integer (kind=int_kind) :: count_bot, count_top 
    ! Keep track of how many layers we change from the bottom and from the top
    real (kind=real_kind), parameter :: tol = 3.5e-14 ! Tolerance for floating point
    ! comparisons.

    do ll = 1, qsize
       do ii = 1, nx
          do jj = 1, nx
             count_bot = 0
             count_top = 0
             do kk = 1, (nlev-1)/2
                lbnd = bnds(1, ii, jj, kk, ll)
                ubnd = bnds(2, ii, jj, kk, ll)
                mr = Qdp(ii, jj, kk, ll)/dp2(ii, jj, kk)
                ! bot->top
                if (mr - ubnd > tol) then ! Too big, give to
                   ! next layer up.
                   mr_diff = mr - ubnd
                   ms_diff = mr_diff * dp2(ii, jj, kk)
                   Qdp(ii, jj, kk, ll) = Qdp(ii, jj, kk, ll) - ms_diff
                   Qdp(ii, jj, kk+1, ll) = Qdp(ii, jj, kk+1, ll) + ms_diff
                   count_bot = count_bot + 1
                else if (lbnd - mr > tol) then ! Too small, take
                   ! from next layer up.
                   mr_diff = lbnd - mr
                   ms_diff = mr_diff * dp2(ii, jj, kk)
                   Qdp(ii, jj, kk, ll) = Qdp(ii, jj, kk, ll) + ms_diff
                   Qdp(ii, jj, kk+1, ll) = Qdp(ii, jj, kk+1, ll) - ms_diff
                   count_bot = count_bot + 1
                end if

                ! top->bot
                lbnd = bnds(1, ii, jj, (nlev-kk+1), ll)
                ubnd = bnds(2, ii, jj, (nlev-kk+1), ll)
                mr = Qdp(ii, jj, (nlev-kk+1), ll)/dp2(ii, jj, (nlev-kk+1))
                if (mr - ubnd > tol) then ! Too big,
                   ! give to next layer down.
                   mr_diff = mr - ubnd
                   ms_diff = mr_diff * dp2(ii, jj, (nlev-kk+1))
                   Qdp(ii, jj, (nlev-kk+1), ll) &
                      & = Qdp(ii, jj, (nlev-kk+1), ll) - ms_diff
                   Qdp(ii, jj, (nlev-kk), ll) &
                      & = Qdp(ii, jj, (nlev-kk), ll) + ms_diff
                   count_top = count_top + 1
                else if (lbnd - mr > tol)  then
                   ! Too small, take from next layer down.
                   mr_diff = lbnd - mr
                   ms_diff = mr_diff * dp2(ii, jj, (nlev-kk+1))
                   Qdp(ii, jj, (nlev-kk+1), ll) &
                      & = Qdp(ii, jj, (nlev-kk+1), ll) + ms_diff
                   Qdp(ii, jj, (nlev-kk), ll) &
                      & = Qdp(ii, jj, (nlev-kk), ll) - ms_diff
                   count_top = count_top + 1
                end if
             end do
         
          ! print '(A, I4, I4, A, I4)', ' ~~ count_top, bot: ', count_top, &
          !   & count_bot, ' of ', nlev
          end do
       end do
    end do

  end subroutine mass_borrow

!=============================================================================
! Returns a convex combination of q_alg = 10, q_alg = 11 to obey global
! bounds.
!=============================================================================

  subroutine conv_comb(Qdp10, Qdp11, nx, nlev, qsize, dp2, bnds)
    
    implicit none

    integer (kind=int_kind), intent(in)  :: nx, nlev, qsize
    real (kind=real_kind), intent(in)    :: Qdp10(nx, nx, nlev, qsize)
    real (kind=real_kind), intent(inout) :: Qdp11(nx, nx, nlev, qsize)
    real (kind=real_kind), intent(in)    :: dp2(nx, nx, nlev)
    real (kind=real_kind), intent(in)    :: bnds(2, nx, nx, qsize) ! Min/Max for each q
    
    !================
    ! Local variables
    !================
    real (kind=real_kind)   :: mr10(nx, nx, nlev, qsize), &
         & mr11(nx, nx, nlev, qsize) ! Mixing ratios
    real (kind=real_kind)   :: mr_diff ! Difference between mixing ratio and
    ! bound of interest.
    integer (kind=int_kind) :: ii, jj, kk, ll ! Iterators.
    real (kind=real_kind)   :: theta ! Coefficient for combination 
    ! Keep track of how many layers we change from the bottom and from the top
    real (kind=real_kind), parameter :: tol = 3.5e-14 ! Tolerance for floating point
    ! comparisons.
    real (kind=real_kind)   :: maxmr, minmr ! Maximum and minimum mixing ratios of
    ! original data.
    real (kind=real_kind)   :: max11, min11 ! Maximum and minimum of remapped qalg11
    real (kind=real_kind)   :: ext10 ! Extrememum for mr10
    integer (kind=int_kind) :: theta_cnt(4) ! Count of theta = 0, 1, ne 0 & 1, and total

    ! Get mixing ratios
    do ii = 1, qsize
       mr10(:,:,:,ii) = Qdp10(:,:,:,ii)/dp2(:,:,:)
       mr11(:,:,:,ii) = Qdp11(:,:,:,ii)/dp2(:,:,:)
    end do

    theta_cnt = 0

    do ii = 1, qsize
       do jj = 1, nx
          do kk = 1, nx
             max11 = maxval(mr11(kk, jj, :, ii))
             min11 = minval(mr11(kk, jj, :, ii))

             maxmr = bnds(2, kk, jj, ii)
             minmr = bnds(1, kk, jj, ii)
             
             if ((max11 - maxmr < tol) .and. (minmr - min11 < tol)) then
                ! Column within bounds, skip
                theta = 1
             else
                if (max11 - maxmr > minmr - min11) then
                   ! Upper bound is exceeded more than lower bound
                   ext10 = maxval(mr10(kk, jj, :, ii))
                   if (abs(max11 - ext10) .lt. tol) then
                      ! Extrema are equal, go with unlimited solution
                      theta = 1
                   else
                      theta = (maxmr - ext10) / (max11 - ext10)
                   end if
                else
                   ext10 = minval(mr10(kk, jj, :, ii))
                   if (abs(ext10 - min11) .lt. tol) then
                      ! Extrema are equal, go with unlimited solution
                      theta = 1
                   else
                      theta = (minmr - ext10) / (min11 - ext10)
                   end if
                end if            
             end if

             theta_cnt(4) = theta_cnt(4) + 1
             if (abs(theta) .lt. tol) then
               theta_cnt(1) = theta_cnt(1) + 1
             else if (abs(theta - 1) .lt. tol) then
               theta_cnt(2) = theta_cnt(2) + 1
             else
               theta_cnt(3) = theta_cnt(3) + 1
             end if

             if (theta .gt. 1) then
               theta = 1
             else if (theta .lt. 0) then
               theta = 0
             end if

             Qdp11(kk, jj, :, ii) = theta * Qdp11(kk, jj, :, ii) &
                  & + (1 - theta) * Qdp10(kk, jj, :, ii)

          end do
       end do
    end do
 
    ! print '(A,4I5)', ' ~~ theta_cnt: 0, 1, ~= 0 and 1, total: ', theta_cnt
    
  end subroutine conv_comb

!=============================================================================
! Returns the domain of dependence of each cell in the target grid
!=============================================================================

  subroutine get_dom_of_dep(ncell, dp1, dp2, ext, cell_dod)

    implicit none

    integer (kind=int_kind), intent(in)  :: ncell
    real (kind=real_kind),   intent(in)  :: dp1(ncell), dp2(ncell)
    integer (kind=int_kind), intent(in)  :: ext ! How far to extend DoD past
    ! cell intersections.
    integer (kind=int_kind), intent(out) :: cell_dod(2, ncell)

    real (kind=real_kind), parameter :: tol = 3.5e-15
    real (kind=real_kind) :: grid1(0:ncell), grid2(0:ncell)
    integer (kind=int_kind) :: jj, kk, kidx

    ! Reconstruct grid1, grid2
    grid1(0) = 0.0_real_kind
    grid2(0) = 0.0_real_kind
    do kk = 1, ncell
       grid1(kk) = dp1(kk) + grid1(kk-1)
       grid2(kk) = dp2(kk) + grid2(kk-1)
    end do

    ! Iterate through the source and target grid to get the location
    ! of the left and right target cell interfaces on the old grid.
    cell_dod(1, 1) = 1 ! Left bound of first new cell is in first old cell
    kidx = 1 ! Cell of old grid were are checking the right boundary of
    do jj = 1, ncell-1 ! Loop through all cells of new grid
       do kk = kidx, ncell ! Loop through all cell boundaries of new grid
          if (grid2(jj) .le. grid1(kk)) then ! In the kkth cell
             cell_dod(2, jj) = kk
             cell_dod(1, jj+1) = kk ! Left bound of next cell is same as 
                                    ! right boundary of current cell
             kidx = kk
             exit
          end if
       end do
    end do
    cell_dod(2, ncell) = ncell ! Right bound of last new cell in last old cell.

    ! Domain of dependence might be further out than intersections.
    do kk = 1, ncell
       ! If we can move the dod by the ext amount while staying in the 
       ! domain, we do so.
       if (cell_dod(1, kk) - ext .ge. 1) then
          cell_dod(1, kk) = cell_dod(1, kk) - ext
       else
          cell_dod(1, kk) = 1
       end if

       if (cell_dod(2, kk) + ext .le. ncell) then
          cell_dod(2, kk) = cell_dod(2, kk) + ext
       else
          cell_dod(2, kk) = ncell
       end if
    end do

      
  end subroutine get_dom_of_dep

!=============================================================================
! Checks that the remapped solution obeys the global bounds. 
!=============================================================================

  subroutine is_bounded(hybrid, Qdp, nx, nlev, qsize, dp2, bnds, res)
    
    use hybrid_mod, only: hybrid_t

    implicit none

    type (hybrid_t),         intent(in)  :: hybrid  ! distributed parallel structure (shared)
    integer (kind=int_kind), intent(in)  :: nx, nlev, qsize
    real (kind=real_kind),   intent(in)  :: Qdp(nx, nx, nlev, qsize)
    real (kind=real_kind),   intent(in)  :: dp2(nx, nx, nlev)
    real (kind=real_kind),   intent(in)  :: bnds(2, nx, nx, nlev, qsize) ! Min/Max for each q
    logical,                 intent(out) :: res
       
    !================
    ! Local variables
    !================
    real (kind=real_kind)   :: mr
    real (kind=real_kind)   :: lbnd, ubnd
    logical                 :: flag
    integer (kind=int_kind) :: ii, jj, kk, ll ! Iterators.
    real (kind=real_kind), parameter :: tol = 3.5e-14 ! Tolerance for floating point
    ! comparisons.

    res = .true.

    do ll = 1, qsize
       do ii = 1, nx
          do jj = 1, nx
             flag = .true.
             do kk = 1, nlev
                lbnd = bnds(1, ii, jj, kk, ll)
                ubnd = bnds(2, ii, jj, kk, ll)

                mr = Qdp(ii, jj, kk, ll) / dp2(ii, jj, kk)
                
                if ( ((mr - ubnd) > tol) .or. ((lbnd - mr) > tol) ) then
                   res = .false.
                   flag = .false.
                end if
             end do
             if (hybrid%masterthread .and. .not. flag) then
                print *, "~~ Violation found in column: ", ii, jj
             end if
          end do
       end do
    end do

  end subroutine is_bounded

end module vertremap_base

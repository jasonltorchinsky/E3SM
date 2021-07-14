#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module vertremap_mod
  use vertremap_base, only: remap1, remap1_nofilter, mass_borrow, conv_comb, &
                            get_dom_of_dep, is_bounded

  use kinds, only                  : real_kind,int_kind
  use dimensions_mod, only         : np,nlev,qsize,nlevp,npsq
  use hybvcoord_mod, only          : hvcoord_t
  use element_mod, only            : element_t
  use perf_mod, only               : t_startf, t_stopf  ! _EXTERNAL
  use parallel_mod, only           : abortmp, parallel_t
  use control_mod, only            : vert_remap_q_alg
  use eos, only : phi_from_eos
  implicit none
  private
  public :: vertical_remap

contains


  subroutine vertical_remap(hybrid,elem,hvcoord,dt,np1,np1_qdp,nets,nete)

  ! This routine is called at the end of the vertically Lagrangian
  ! dynamics step to compute the vertical flux needed to get back
  ! to reference eta levels
  !
  ! input:
  !     derived%dp()  delta p on levels at beginning of timestep
  !     state%dp3d(np1)  delta p on levels at end of timestep
  ! output:
  !     state%ps_v(np1)          surface pressure at time np1
  !     derived%eta_dot_dpdn()   vertical flux from final Lagrangian
  !                              levels to reference eta levels
  !
  use kinds,          only: real_kind
  use hybvcoord_mod,  only: hvcoord_t
  use control_mod,    only: rsplit
  use hybrid_mod,     only: hybrid_t
  use physical_constants, only : Cp,g

  type (hybrid_t),  intent(in)    :: hybrid  ! distributed parallel structure (shared)
  type (element_t), intent(inout) :: elem(:)
  type (hvcoord_t)                :: hvcoord
  real (kind=real_kind)           :: dt

  integer :: ie,i,j,k,np1,nets,nete,np1_qdp
  integer :: q

  real (kind=real_kind), dimension(np,np,nlev)  :: dp,dp_star
  real (kind=real_kind), dimension(np,np,nlevp) :: phi_ref
  real (kind=real_kind), dimension(np,np,nlev,5)  :: ttmp

  real (kind=real_kind), dimension(2,np,np,nlev) :: theta_ext ! Min, max of temperature before remapping
  real (kind=real_kind), dimension(np,np,nlev)  :: theta_tmp ! Additional array for the linear combination
  integer (kind = int_kind), dimension(2, nlev) :: cell_dod ! Domain of dependence of each target cell
  logical :: bflag = .true.
  call t_startf('vertical_remap')

  ! reference levels:
  !   dp(k) = (hyai(k+1)-hyai(k))*ps0 + (hybi(k+1)-hybi(k))*ps_v(i,j)
  !   hybi(1)=0          pure pressure at top of atmosphere
  !   hyai(1)=ptop
  !   hyai(nlev+1) = 0   pure sigma at bottom
  !   hybi(nlev+1) = 1
  !
  ! sum over k=1,nlev
  !  sum(dp(k)) = (hyai(nlev+1)-hyai(1))*ps0 + (hybi(nlev+1)-hybi(1))*ps_v
  !             = -ps0 + ps_v
  !  ps_v =  ps0+sum(dp(k))
  !
  ! reference levels:
  !    dp(k) = (hyai(k+1)-hyai(k))*ps0 + (hybi(k+1)-hybi(k))*ps_v
  ! floating levels:
  !    dp_star(k) = dp(k) + dt_q*(eta_dot_dpdn(i,j,k+1) - eta_dot_dpdn(i,j,k) )
  ! hence:
  !    (dp_star(k)-dp(k))/dt_q = (eta_dot_dpdn(i,j,k+1) - eta_dot_dpdn(i,j,k) )
  !
   do ie=nets,nete
     ! update final ps_v
     elem(ie)%state%ps_v(:,:,np1) = hvcoord%hyai(1)*hvcoord%ps0 + &
          sum(elem(ie)%state%dp3d(:,:,:,np1),3)
     do k=1,nlev
        dp(:,:,k) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
             ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps_v(:,:,np1)
        if (rsplit==0) then
           dp_star(:,:,k) = dp(:,:,k) + dt*(elem(ie)%derived%eta_dot_dpdn(:,:,k+1) -&
                elem(ie)%derived%eta_dot_dpdn(:,:,k))
        else
           dp_star(:,:,k) = elem(ie)%state%dp3d(:,:,k,np1)
        endif
     enddo
     if (minval(dp_star)<0) then
        do k=1,nlev
        do i=1,np
        do j=1,np
           if (dp_star(i,j,k ) < 0) then
              print *,'index ie,i,j,k = ',ie,i,j,k
              print *,'dp,dp_star = ',dp(i,j,k),dp_star(i,j,k)
              print *,'eta_dot_dpdn = ',elem(ie)%derived%eta_dot_dpdn(i,j,k+1),elem(ie)%derived%eta_dot_dpdn(i,j,k)
              print *,"column location lat,lon (radians):",elem(ie)%spherep(i,j)%lat,elem(ie)%spherep(i,j)%lon
           endif
        enddo
        enddo
        enddo
        call abortmp('negative layer thickness.  timestep or remap time too large')
     endif

     if (rsplit>0) then
        ! remove hydrostatic phi befor remap
        call phi_from_eos(hvcoord,elem(ie)%state%phis,elem(ie)%state%vtheta_dp(:,:,:,np1),dp_star,phi_ref)
        elem(ie)%state%phinh_i(:,:,:,np1)=&
             elem(ie)%state%phinh_i(:,:,:,np1) -phi_ref(:,:,:)
 
        !  REMAP u,v,T from levels in dp3d() to REF levels
        ttmp(:,:,:,1)=elem(ie)%state%v(:,:,1,:,np1)*dp_star
        ttmp(:,:,:,2)=elem(ie)%state%v(:,:,2,:,np1)*dp_star
        ttmp(:,:,:,3)=elem(ie)%state%vtheta_dp(:,:,:,np1)
        do k=1,nlev
           ttmp(:,:,k,4)=elem(ie)%state%phinh_i(:,:,k+1,np1)-&
                elem(ie)%state%phinh_i(:,:,k,np1) 
           ttmp(:,:,k,5)=elem(ie)%state%w_i(:,:,k+1,np1)-&
                elem(ie)%state%w_i(:,:,k,np1)
        enddo
        !ttmp(:,:,:,4)=ttmp(:,:,:,4) !*dp_star
        !ttmp(:,:,:,5)=ttmp(:,:,:,5) !*dp_star
    
        call t_startf('vertical_remap1_1')
        if (hybrid%masterthread) then
           print *, "Original T col 1,1: ", ttmp(1,1,:,3)/dp(1,1,:)
        endif
        if (vert_remap_q_alg == 20) then
          ! Use q_alg = 11 for temperature only
          call remap1(ttmp(:,:,:,1:2),np,2,dp_star,dp,10)
          call remap1(ttmp(:,:,:,4:5),np,2,dp_star,dp,10)
          call remap1(ttmp(:,:,:,3),np,1,dp_star,dp,11)
        else if (vert_remap_q_alg == 21) then
          ! Use q_alg = 11 + local-bounds mass-borrowing for temperature only
          call remap1(ttmp(:,:,:,1:2),np,2,dp_star,dp,10)
          call remap1(ttmp(:,:,:,4:5),np,2,dp_star,dp,10)
          call get_dom_of_dep(nlev, dp_star, dp, 2, cell_dod)
          do i = 1, np
            do j = 1, np
              do k = 1, nlev
                theta_ext(1,i,j,k) = minval( ttmp(i,j,cell_dod(1,k):cell_dod(2,k),3) &
                                             & / dp_star(i,j,cell_dod(1,k):cell_dod(2,k)) )
                theta_ext(2,i,j,k) = maxval( ttmp(i,j,cell_dod(1,k):cell_dod(2,k),3) &
                                             & / dp_star(i,j,cell_dod(1,k):cell_dod(2,k)) )
              end do
            end do
          end do
          call remap1(ttmp(:,:,:,3),np,1,dp_star,dp,11)
          call mass_borrow(ttmp(:,:,:,3),np,nlev,1,dp,theta_ext)

          call is_bounded(hybrid, ttmp(:,:,:,3), np, nlev, 1, dp, theta_ext, bflag)

          if (hybrid%masterthread) then
             print *, 'Remapped solution is bounded: ', bflag
          end if

        else if (vert_remap_q_alg == 22) then
          ! Use q_alg = 11 + global-bounds mass-borrowing for temperature only
          call remap1(ttmp(:,:,:,1:2),np,2,dp_star,dp,10)
          call remap1(ttmp(:,:,:,4:5),np,2,dp_star,dp,10)
          do i = 1, np
            do j = 1, np
              theta_ext(1,i,j,:) = minval( ttmp(i,j,:,3) / dp_star(i,j,:) )
              theta_ext(2,i,j,:) = maxval( ttmp(i,j,:,3) / dp_star(i,j,:) )
            end do
          end do

          if (hybrid%masterthread) then
             print *, " Extrema for column 1: ", theta_ext(:,1,1,1)
          end if

          call remap1(ttmp(:,:,:,3),np,1,dp_star,dp,11)
          call mass_borrow(ttmp(:,:,:,3),np,nlev,1,dp,theta_ext)

          call is_bounded(hybrid, ttmp(:,:,:,3), np, nlev, 1, dp, theta_ext, bflag)

          if (hybrid%masterthread) then
             print *, 'Remapped solution is bounded: ', bflag
          end if

        else if (vert_remap_q_alg == 23) then
          ! Use q_alg = 11 + linear combination for temperature only
          vert_remap_q_alg = 10
          call remap1(ttmp(:,:,:,1:2),np,2,dp_star,dp,vert_remap_q_alg)
          call remap1(ttmp(:,:,:,4:5),np,2,dp_star,dp,vert_remap_q_alg)

          do i = 1, np
            do j = 1, np
              theta_ext(1,i,j,:) = minval(ttmp(i,j,:,3)/dp_star(i,j,:))
              theta_ext(2,i,j,:) = maxval(ttmp(i,j,:,3)/dp_star(i,j,:))
            end do
          end do
          theta_tmp = ttmp(:,:,:,3)
          call remap1(theta_tmp,np,1,dp_star,dp,10)
          call remap1(ttmp(:,:,:,3),np,1,dp_star,dp,11)
          call conv_comb(theta_tmp,ttmp(:,:,:,3),np,nlev,1,dp,theta_ext)
        else
          do i = 1, np
            do j = 1, np
              theta_ext(1,i,j,:) = minval( ttmp(i,j,:,3) / dp_star(i,j,:) )
              theta_ext(2,i,j,:) = maxval( ttmp(i,j,:,3) / dp_star(i,j,:) )
            end do
          end do

          call remap1(ttmp,np,5,dp_star,dp,vert_remap_q_alg)

          call is_bounded(hybrid, ttmp(:,:,:,3), np, nlev, 1, dp, theta_ext, bflag)

          if (hybrid%masterthread) then
             print *, 'Remapped solution is bounded: ', bflag
          end if
        end if
        if (hybrid%masterthread) then
           print *, "Remapped T col 1,1: ", ttmp(1,1,:,3)/dp(1,1,:)
        endif
        call t_stopf('vertical_remap1_1')

        elem(ie)%state%v(:,:,1,:,np1)=ttmp(:,:,:,1)/dp
        elem(ie)%state%v(:,:,2,:,np1)=ttmp(:,:,:,2)/dp
        elem(ie)%state%vtheta_dp(:,:,:,np1)=ttmp(:,:,:,3)
      
        
        do k=nlev,1,-1
           elem(ie)%state%phinh_i(:,:,k,np1)=&
                elem(ie)%state%phinh_i(:,:,k+1,np1)-ttmp(:,:,k,4)  !/dp(:,:,k)
           elem(ie)%state%w_i(:,:,k,np1)=&
                elem(ie)%state%w_i(:,:,k+1,np1)-ttmp(:,:,k,5)  !/dp(:,:,k)
        enddo

        ! depends on theta, so do this after updating theta:
        call phi_from_eos(hvcoord,elem(ie)%state%phis,elem(ie)%state%vtheta_dp(:,:,:,np1),dp,phi_ref)
        elem(ie)%state%phinh_i(:,:,:,np1)=&
             elem(ie)%state%phinh_i(:,:,:,np1)+phi_ref(:,:,:)


        ! since u changed, update w b.c.:
        elem(ie)%state%w_i(:,:,nlevp,np1) = (elem(ie)%state%v(:,:,1,nlev,np1)*elem(ie)%derived%gradphis(:,:,1) + &
             elem(ie)%state%v(:,:,2,nlev,np1)*elem(ie)%derived%gradphis(:,:,2))/g
       
     endif

     ! remap the gll tracers from lagrangian levels (dp_star)  to REF levels dp
     if (qsize>0 .and. np1_qdp > 0) then

       call t_startf('vertical_remap1_3')
       call remap1(elem(ie)%state%Qdp(:,:,:,:,np1_qdp),np,qsize,dp_star,dp,vert_remap_q_alg)
       call t_stopf('vertical_remap1_3')

       !dir$ simd
       do q=1,qsize
          elem(ie)%state%Q(:,:,:,q)=elem(ie)%state%Qdp(:,:,:,q,np1_qdp)/dp(:,:,:)
       enddo
     endif

     ! reinitialize dp3d after remap
     elem(ie)%state%dp3d(:,:,:,np1)=dp(:,:,:)

  enddo
  call t_stopf('vertical_remap')
  end subroutine vertical_remap


end module 





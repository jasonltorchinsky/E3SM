#ifndef SHOC_DIAG_THIRD_SHOC_MOMENTS_IMPL_HPP
#define SHOC_DIAG_THIRD_SHOC_MOMENTS_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc diag_third_shoc_moments. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::diag_third_shoc_moments(
  const MemberType&            team,
  const Int&                   nlev,
  const Int&                   nlevi,
  const uview_1d<const Spack>& w_sec,
  const uview_1d<const Spack>& thl_sec,
  const uview_1d<const Spack>& wthl_sec,
  const uview_1d<const Spack>& isotropy,
  const uview_1d<const Spack>& brunt,
  const uview_1d<const Spack>& thetal,
  const uview_1d<const Spack>& tke,
  const uview_1d<const Spack>& dz_zt,
  const uview_1d<const Spack>& dz_zi,
  const uview_1d<const Spack>& zt_grid,
  const uview_1d<const Spack>& zi_grid,
  const uview_1d<Spack>&       w_sec_zi,
  const uview_1d<Spack>&       isotropy_zi,
  const uview_1d<Spack>&       brunt_zi,
  const uview_1d<Spack>&       thetal_zi,
  const uview_1d<Spack>&       w3)
{
  // Constants
  const auto largeneg = SC::largeneg;
  const auto mintke = SC::mintke;

  // Interpolate variables onto the interface levels
  linear_interp(team,zt_grid,zi_grid,isotropy,isotropy_zi,nlev,nlevi,0);
  linear_interp(team,zt_grid,zi_grid,brunt,brunt_zi,nlev,nlevi,largeneg);
  linear_interp(team,zt_grid,zi_grid,w_sec,w_sec_zi,nlev,nlevi,sp(2.0/3.0)*mintke);
  linear_interp(team,zt_grid,zi_grid,thetal,thetal_zi,nlev,nlevi,0);
  team.team_barrier();

  //Diagnose the third moment of the vertical-velocity
  compute_diag_third_shoc_moment(team,nlev,nlevi,w_sec,thl_sec,wthl_sec,
                                 tke, dz_zt, dz_zi,isotropy_zi, brunt_zi,
                                 w_sec_zi,thetal_zi,w3);

  // Perform clipping to prevent unrealistically large values from occuring
  clipping_diag_third_shoc_moments(team,nlevi,w_sec_zi,w3);
}

} // namespace shoc
} // namespace scream

#endif

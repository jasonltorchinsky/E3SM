#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "share/util/scream_utils.hpp"
#include "share/scream_kokkos.hpp"
#include "share/scream_pack.hpp"
#include "physics/p3/p3_functions.hpp"
#include "physics/p3/p3_functions_f90.hpp"
#include "physics/p3/p3_f90.hpp"
#include "share/util/scream_kokkos_utils.hpp"

#include "p3_unit_tests_common.hpp"

#include <thread>
#include <array>
#include <algorithm>
#include <random>

namespace scream {
namespace p3 {
namespace unit_test {

// Derive rough bounds on table domain:
//
// First, mu_r is in [0,9].
//
// Second, start with
//     dum1 = (1 + mu_r)/lamr
// Let's get the lower limit on lamr as a function of mu first. So
//     1 = rdumii = (dum1 1e6 + 5) inv_dum3, inv_dum3 = 0.1
//       = ((mu_r+1)/lamr 1e6 + 5) 0.1
//     (1/0.1 - 5) 1e-6 lamr = 1 + mu_r
//     lamr = 0.2e6 (1 + mu_r).
// Now the upper limit:
//     300 = (dum1 1e6 + 5) inv_dum3, inv_dum3 = 1/30
//         = ((mu_r+1)/lamr 1e6 - 195) 1/30 + 20
//     => lamr = 1e6 (1 + mu_r) / (30 (300 - 20) + 195)
//             ~ 116.3 (1 + mu_r)
//
// This test checks for smoothness (C^0 function) of the interpolant both inside
// the domains over which the tables are defined and outside them. The primary
// tool is measuring the maximum slope magnitude as a function of mesh
// refinement, where the mesh is a 1D mesh transecting the table domain.

template <typename D>
struct UnitWrap::UnitTest<D>::TestTable3 {

  KOKKOS_FUNCTION static Scalar calc_lamr (const Scalar& mu_r, const Scalar& alpha) {
    // Parameters for lower and upper bounds, derived above, multiplied by
    // factors so we go outside of the bounds a bit. Using these, we map alpha
    // in [0,1] -> meaningful lamr.
    const Scalar lamr_lo = 0.1*116.3, lamr_hi = 1.5*0.2e6;
    return (lamr_lo + alpha*(lamr_hi - lamr_lo))*(1 + mu_r);
  }

  KOKKOS_FUNCTION static Scalar calc_mu_r (const Scalar& alpha) {
    return alpha*10;
  }

  // Perform the table lookup and interpolation operations for (mu_r, lamr).
  KOKKOS_FUNCTION static Spack interp (const view_2d_table& table, const Scalar& mu_r,
                                       const Scalar& lamr) {
    // Init the pack to all the same value, and compute in every pack slot.
    Smask qr_gt_small(true);
    Spack mu_r_p(mu_r), lamr_p(lamr);
    Table3 t3;
    Functions::lookup(qr_gt_small, mu_r_p, lamr_p, t3);
    Spack val(qr_gt_small, Functions::apply_table(qr_gt_small, table, t3));
    return val;
  }

  static void run () {
    // This test doesn't use mu_r_table, as that is not a table3 type. It
    // doesn't matter whether we use vm_table or vn_table, as the table values
    // don't matter in what we are testing; we are testing interpolation
    // procedures that are indepennt of particular table values.
    view_1d_table mu_r_table;
    view_2d_table vn_table, vm_table, revap_table;
    view_dnu_table dnu;
    Functions::init_kokkos_tables(vn_table, vm_table, revap_table, mu_r_table, dnu);

    // Estimate two maximum slope magnitudes for two meshes, the second 10x
    // refined w.r.t. the first.
    Real slopes[2];
    const Int nslopes = sizeof(slopes)/sizeof(*slopes);
    Int N;

    // Study the mu_r direction.
    // For a sequence of refined meshes:
    N = 1000;
    for (Int refine = 0; refine < nslopes; ++refine) {
      // Number of cells in the mesh.
      N *= 10;
      // Cell size relative to a parameter domain of 1.
      const Scalar delta = 1.0/N;

      // Compute the slope magnitude at a specific (mu_r, lamr) in the mu_r
      // direction.
      const Scalar lamr = calc_lamr(4.5, 0.5);
      const auto get_max_slope = KOKKOS_LAMBDA (const Int& i, Scalar& slope) {
        // Interpolate at a specific (mu_r, lamr).
        const auto eval = [&] (const Int& i) {
          const auto alpha = double(i)/N;
          const auto mu_r = calc_mu_r(alpha);
          const auto val = interp(vm_table, mu_r, lamr);
          return std::log(val[0]);
        };
        slope = util::max(slope, std::abs((eval(i+1) - eval(i))/delta));
      };

      Scalar max_slope;
      Kokkos::parallel_reduce(RangePolicy(0, N), get_max_slope,
                              Kokkos::Max<Scalar>(max_slope));
      Kokkos::fence();
      slopes[refine] = max_slope;
    }

    // Now that we have collected slopes as a function of 10x mesh refinement,
    // determine whether the slope estimates are converging. If they are, we can
    // conclude there are no discontinuities. If they are not, then we are sensing
    // a discontinuity, which is a bug.
    //   In detail, for a 10x mesh refinement, a good slope growth rate is right
    // around 1, and a bad one is is roughly 10. We set the threshold at 1.1.
    const auto check_growth = [&] (const std::string& label, const Scalar& growth) {
      bool bad_growth = growth > 1.1;
      if (bad_growth) {
        std::cout << "Table3 FAIL: Slopes in the " << label << " direction are "
        << slopes[0] << " and " << slopes[1]
        << ", which grows by factor " << growth
        << ". Near 1 is good; near 10 is bad.\n";
      }
      REQUIRE(!bad_growth);
    };
    check_growth("mu_r", slopes[1]/slopes[0]);

    // Study the lamr direction.
    N = 4000;
    for (Int refine = 0; refine < nslopes; ++refine) {
      N *= 2;
      const Scalar delta = 1.0/N;

      // Compute the slope magnitude at a specific (mu_r, lamr) in the lamr
      // direction.
      const Scalar mu_r = 3.5;
      const auto get_max_slope = KOKKOS_LAMBDA (const Int& i, Scalar& slope) {
        // Interpolate at a specific (mu_r, lamr).
        const auto eval = [&] (const Int& i) {
          const auto alpha = double(i)/N;
          const auto lamr = calc_lamr(mu_r, alpha);
          const auto val = interp(vm_table, mu_r, lamr);
          return std::log(val[0]);
        };
        slope = util::max(slope, std::abs((eval(i+1) - eval(i))/delta));
      };

      Scalar max_slope;
      Kokkos::parallel_reduce(RangePolicy(0, N), get_max_slope,
                              Kokkos::Max<Scalar>(max_slope));
      Kokkos::fence();
      slopes[refine] = max_slope;
    }
    check_growth("lamr", slopes[1]/slopes[0]);
  }
};

}
}
}

namespace {

TEST_CASE("p3_tables", "[p3_functions]")
{
  scream::p3::p3_init(true); // need fortran table data

  scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestTable3::run();
}

} // namespace

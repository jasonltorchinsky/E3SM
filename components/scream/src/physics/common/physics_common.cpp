#include "physics_common.hpp"
#include "ekat/scream_assert.hpp"
#include "ekat/util/scream_utils.hpp"
#include "ekat/util/scream_kokkos_utils.hpp"
#include "ekat/scream_pack_kokkos.hpp"

#include <random>

using scream::Real;
using scream::Int;

namespace scream {

namespace physics {

//
// A C++ interface to physics-related fortran calls and vice versa
//

// Cuda implementations of std math routines are not necessarily BFB
// with the host.
template <typename ScalarT, typename DeviceT>
struct CudaWrap
{
  using Scalar = ScalarT;

  static Scalar cxx_pow(Scalar base, Scalar exp)
  {
    Scalar result;
    Kokkos::parallel_reduce(1, KOKKOS_LAMBDA(const Int&, Scalar& value) {
        value = std::pow(base, exp);
    }, result);

    return result;
  }

#define cuda_wrap_single_arg(wrap_name, func_call)      \
static Scalar wrap_name(Scalar input) {                 \
  Scalar result;                                        \
  Kokkos::parallel_reduce(1, KOKKOS_LAMBDA(const Int&, Scalar& value) { \
    value = func_call(input);                                         \
  }, result);                                                         \
  return result;                                                      \
}

  cuda_wrap_single_arg(cxx_gamma, std::tgamma)
  cuda_wrap_single_arg(cxx_sqrt, std::sqrt)
  cuda_wrap_single_arg(cxx_cbrt, std::cbrt)
  cuda_wrap_single_arg(cxx_log, std::log)
  cuda_wrap_single_arg(cxx_log10, std::log10)
  cuda_wrap_single_arg(cxx_exp, std::exp)
  cuda_wrap_single_arg(cxx_tanh, std::tanh)

#undef cuda_wrap_single_arg
};

extern "C" {

Real cxx_pow(Real base, Real exp)
{
#ifdef KOKKOS_ENABLE_CUDA
  return CudaWrap<Real, DefaultDevice>::cxx_pow(base, exp);
#else
  return std::pow(base, exp);
#endif
}

Real cxx_gamma(Real input)
{
#ifdef KOKKOS_ENABLE_CUDA
  return CudaWrap<Real, DefaultDevice>::cxx_gamma(input);
#else
  return std::tgamma(input);
#endif
}

Real cxx_cbrt(Real input)
{
#ifdef KOKKOS_ENABLE_CUDA
  return CudaWrap<Real, DefaultDevice>::cxx_cbrt(input);
#else
  return std::cbrt(input);
#endif
}

Real cxx_sqrt(Real input)
{
#ifdef KOKKOS_ENABLE_CUDA
  return CudaWrap<Real, DefaultDevice>::cxx_sqrt(input);
#else
  return std::sqrt(input);
#endif
}

Real cxx_log(Real input)
{
#ifdef KOKKOS_ENABLE_CUDA
  return CudaWrap<Real, DefaultDevice>::cxx_log(input);
#else
  return std::log(input);
#endif
}

Real cxx_log10(Real input)
{
#ifdef KOKKOS_ENABLE_CUDA
  return CudaWrap<Real, DefaultDevice>::cxx_log10(input);
#else
  return std::log10(input);
#endif
}

Real cxx_exp(Real input)
{
#ifdef KOKKOS_ENABLE_CUDA
  return CudaWrap<Real, DefaultDevice>::cxx_exp(input);
#else
  return std::exp(input);
#endif
}

Real cxx_tanh(Real input)
{
#ifdef KOKKOS_ENABLE_CUDA
  return CudaWrap<Real, DefaultDevice>::cxx_tanh(input);
#else
  return std::tanh(input);
#endif
}

} // end extern "C"

} // namespace physics
} // namespace scream

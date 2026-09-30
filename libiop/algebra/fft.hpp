/**@file
 *****************************************************************************
 Implementation of Gao-Mateer for the additive FFT/IFFT,
 implementation of Nlog(d) Cooley-Tukey for the multiplicative FFT,
 and wrappers to libfqfft for the multiplicative IFFT.
 *****************************************************************************
 * @author     This file is part of libiop (see AUTHORS)
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#ifndef LIBIOP_ALGEBRA_FFT_HPP_
#define LIBIOP_ALGEBRA_FFT_HPP_

#include <cstdint>
#include <vector>

#include "libiop/algebra/field_subset/field_subset.hpp"
#include "libiop/algebra/field_subset/subspace.hpp"
#include "libiop/algebra/field_subset/subgroup.hpp"
#include <libff/common/utils.hpp>

namespace libiop {

/** When the evaluation domain uses the Cantor special basis, additive_FFT_wrapper
 *  can call either the LCH implementation or cantor:: (depends/additive-fft/C++/Cantor).
 *  Default is LCH. Benchmarks may switch per thread without recompiling. */
enum class additive_fft_cantor_implementation : std::uint8_t {
    lch_afft = 0,
    cantor_afft = 1,
};

inline additive_fft_cantor_implementation &additive_fft_cantor_impl_storage()
{
    thread_local additive_fft_cantor_implementation v = additive_fft_cantor_implementation::lch_afft;
    return v;
}

inline void set_additive_fft_cantor_implementation(additive_fft_cantor_implementation impl)
{
    additive_fft_cantor_impl_storage() = impl;
}

inline additive_fft_cantor_implementation get_additive_fft_cantor_implementation()
{
    return additive_fft_cantor_impl_storage();
}

/* Performs naive computation of the polynomial evaluation
   problem. Mostly useful for testing. */
template<typename FieldT>
std::vector<FieldT> naive_FFT(const std::vector<FieldT> &poly_coeffs,
                              const field_subset<FieldT> &domain);

template<typename FieldT>
std::vector<FieldT> additive_FFT(const std::vector<FieldT> &poly_coeffs,
                                 const affine_subspace<FieldT> &domain);

template<typename FieldT>
std::vector<FieldT> additive_IFFT(const std::vector<FieldT> &evals,
                                  const affine_subspace<FieldT> &domain);

/* Calls additive_FFT but adds trace data */
template<typename FieldT>
std::vector<FieldT> additive_FFT_wrapper(const std::vector<FieldT> &v,
                                         const affine_subspace<FieldT> &H);

/* Calls additive_IFFT but adds trace data */
template<typename FieldT>
std::vector<FieldT> additive_IFFT_wrapper(const std::vector<FieldT> &v,
                                          const affine_subspace<FieldT> &H);

template<typename FieldT>
std::vector<FieldT> multiplicative_FFT(const std::vector<FieldT> &poly_coeffs,
                                       const multiplicative_coset<FieldT> &domain);

template<typename FieldT>
std::vector<FieldT> multiplicative_IFFT(const std::vector<FieldT> &evals,
                                        const multiplicative_coset<FieldT> &domain);

template<typename FieldT>
std::vector<FieldT> multiplicative_FFT_wrapper(const std::vector<FieldT> &v,
                                               const multiplicative_coset<FieldT> &H);

template<typename FieldT>
std::vector<FieldT> multiplicative_IFFT_wrapper(const std::vector<FieldT> &v,
                                                const multiplicative_coset<FieldT> &H);

template<typename FieldT>
std::vector<FieldT> FFT_over_field_subset(const std::vector<typename libff::enable_if<libff::is_multiplicative<FieldT>::value, FieldT>::type> &coeffs,
                                            field_subset<FieldT> domain);

template<typename FieldT>
std::vector<FieldT> FFT_over_field_subset(const std::vector<typename libff::enable_if<libff::is_additive<FieldT>::value, FieldT>::type> &coeffs,
                                            field_subset<FieldT> domain);

template<typename FieldT>
std::vector<FieldT> IFFT_over_field_subset(const std::vector<typename libff::enable_if<libff::is_multiplicative<FieldT>::value, FieldT>::type> &evals,
                                             field_subset<FieldT> domain);

template<typename FieldT>
std::vector<FieldT> IFFT_over_field_subset(const std::vector<typename libff::enable_if<libff::is_additive<FieldT>::value, FieldT>::type> &evals,
                                             field_subset<FieldT> domain);

template<typename FieldT>
std::vector<FieldT> IFFT_of_known_degree_over_field_subset(
    const std::vector<typename libff::enable_if<libff::is_multiplicative<FieldT>::value, FieldT>::type> evals,
    size_t degree_bound,
    field_subset<FieldT> domain);

template<typename FieldT>
std::vector<FieldT> IFFT_of_known_degree_over_field_subset(
    const std::vector<typename libff::enable_if<libff::is_additive<FieldT>::value, FieldT>::type> evals,
    size_t degree_bound,
    field_subset<FieldT> domain);

} // namespace libiop

#include "libiop/algebra/fft.tcc"

#endif // LIBIOP_ALGEBRA_FFT_HPP_

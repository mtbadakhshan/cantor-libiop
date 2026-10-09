#ifndef LIBIOP_ALGEBRA_POLYNOMIALS_LCH_POLYNOMIAL_HPP_
#define LIBIOP_ALGEBRA_POLYNOMIALS_LCH_POLYNOMIAL_HPP_
#include <cstddef>
#include <vector>

#include "libiop/algebra/polynomials/lch_standard.hpp"

namespace libiop
{
                    
    template <typename FieldT>
    FieldT lch_polynomial_evaluation_at_point(const std::vector<FieldT> &lch_coeffs,
                                        const FieldT &x);

    /** LCH / Cantor identity: X_{t + 2^eta j'}(x) = X_t(x) * X_{j'}(W_eta(x)).
     *  The FRI fold of f at x_i is therefore the coefficient split
     *  c'_j' = sum_{t < 2^eta} X_t(x_i) * c_{t + 2^eta j'}.
     *  bit0_scale multiplies the weights of odd t, i.e. the bit-0 factor
     *  of the input basis (used when the input is a previous tail vector). */
    template <typename FieldT>
    std::vector<FieldT> lch_fri_fold(const std::vector<FieldT> &lch_coeffs,
                                     const std::size_t eta,
                                     const FieldT &x_i,
                                     const FieldT &bit0_scale = FieldT::one());

    /** Same grouping, with X_t from a per-round basis table (later FRI images). */
    template <typename FieldT>
    std::vector<FieldT> lch_fri_fold(const std::vector<FieldT> &lch_coeffs,
                                     const std::size_t eta,
                                     const FieldT &x_i,
                                     const lch_basis_tables<FieldT> &tables,
                                     const FieldT &bit0_scale = FieldT::one());

    /** LCH subspace-basis coefficients to monomials (LCH inv_basis_conversion). */
    template <typename FieldT>
    std::vector<FieldT> lch_to_monomial(const std::vector<FieldT> &lch_coeffs);

} // namespace libiop

#include "libiop/algebra/polynomials/lch_polynomial.tcc"

#endif // LIBIOP_ALGEBRA_POLYNOMIALS_LCH_POLYNOMIAL_HPP_

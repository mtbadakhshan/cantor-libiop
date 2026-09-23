#ifndef LIBIOP_ALGEBRA_POLYNOMIALS_LCH_VANISHING_POLYNOMIAL_HPP_
#define LIBIOP_ALGEBRA_POLYNOMIALS_LCH_VANISHING_POLYNOMIAL_HPP_
#include <cstddef>
#include <vector>

#include "libiop/algebra/polynomials/polynomial.hpp"

namespace libiop
{
    // Returns the quotient and remainder of P / Z in the LCH basis. The first element in 
    // the pair is the quotient, the latter is the remainder.

    template <typename FieldT>
    std::pair<polynomial<FieldT>,
              polynomial<FieldT>>
    lch_basis_polynomial_over_vanishing_polynomial(const polynomial<FieldT> &P, // polynomial in LCH basis
                                                     const vanishing_polynomial<FieldT> &Z); // vanishing polynomial

    template <typename FieldT>
    FieldT lch_vanishing_polynomial_evaluation_at_point(const FieldT &evalpoint,
                                                         const vanishing_polynomial<FieldT> &Z);

} // namespace libiop

#include "libiop/algebra/polynomials/lch_vanishing_polynomial.tcc"

#endif // LIBIOP_ALGEBRA_POLYNOMIALS_LCH_VANISHING_POLYNOMIAL_HPP_

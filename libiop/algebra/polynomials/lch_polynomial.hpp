#ifndef LIBIOP_ALGEBRA_POLYNOMIALS_LCH_POLYNOMIAL_HPP_
#define LIBIOP_ALGEBRA_POLYNOMIALS_LCH_POLYNOMIAL_HPP_
#include <cstddef>
#include <vector>


namespace libiop
{
                    
    template <typename FieldT>
    FieldT lch_polynomial_evaluation_at_point(const std::vector<FieldT> &lch_coeffs,
                                        const FieldT &x);

} // namespace libiop

#include "libiop/algebra/polynomials/lch_polynomial.tcc"

#endif // LIBIOP_ALGEBRA_POLYNOMIALS_LCH_POLYNOMIAL_HPP_

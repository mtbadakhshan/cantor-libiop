#include <cstddef>
#include <vector>
#include "libiop/algebra/polynomials/polynomial.hpp"
#include "libiop/algebra/polynomials/poly_basis.hpp"

namespace libiop
{

    // template <typename FieldT>
    // std::pair<polynomial<FieldT>,
    //           polynomial<FieldT>>
    // lch_basis_polynomial_over_vanishing_polynomial(const polynomial<FieldT> &P,
    //                                                const vanishing_polynomial<FieldT> &Z)
    // {
    //     if (!Z.is_cantor_basis())
    //     {
    //         throw std::invalid_argument("lch_basis_polynomial_over_vanishing_polynomial is only implemented for the Cantor basis");
    //     }
    //     std::vector<FieldT> P_coeffs = P.coefficients();
    //     const size_t chunk_size = Z.degree();
    //     const size_t n_chunks = (size_t)ceil((double)(P_coeffs.size()) / (double)chunk_size); // we assume these are divisible dividend_size = n_chunks * chunk_size
    //     const size_t m = (size_t)ceil(libff::log2(n_chunks));

    //     size_t actual_quotient_size = P_coeffs.size() - chunk_size;
    //     std::vector<FieldT> quotient_coeffs(P_coeffs.begin() + chunk_size, P_coeffs.end());
    //     std::vector<FieldT> remainder_coeffs(P_coeffs.begin(), P_coeffs.begin() + chunk_size);

    //     for (size_t r = 0; r < m; ++r)
    //     {
    //         size_t s = 1 << r;
    //         for (size_t offset = 0; offset < n_chunks; offset += 2 * s)
    //         {
    //             if (offset + s > n_chunks - 1)
    //                 break;
    //             size_t src_chunk_index = (offset + s - 1) * chunk_size;
    //             for (size_t j = 0; j < chunk_size; ++j)
    //             {
    //                 for (size_t i = 0; i < s - 1; ++i)
    //                 {
    //                     quotient_coeffs[(offset + i) * chunk_size + j] += quotient_coeffs[src_chunk_index + j];
    //                 }
    //             }
    //         }
    //     }
    //     return std::make_pair(polynomial<FieldT>(std::move(quotient_coeffs)),
    //                           polynomial<FieldT>(std::move(remainder_coeffs)));
    // }

    template <typename FieldT>
std::pair<polynomial<FieldT>, polynomial<FieldT>>
lch_basis_polynomial_over_vanishing_polynomial(const polynomial<FieldT> &P,
                                               const vanishing_polynomial<FieldT> &Z)
{
    if (!Z.is_cantor_basis()) {
        throw std::invalid_argument(
            "lch_basis_polynomial_over_vanishing_polynomial is only implemented for the Cantor basis");
    }

    std::vector<FieldT> P_coeffs = P.coefficients();
    const size_t chunk_size = Z.degree();
    if (chunk_size == 0 || (chunk_size & (chunk_size - 1)) != 0) {
        throw std::invalid_argument("LCH divide: |H| must be a power of two");
    }
    if (P_coeffs.size() <= chunk_size) {
        P_coeffs.resize(chunk_size, FieldT::zero());
        return std::make_pair(polynomial<FieldT>(),
                              polynomial<FieldT>(std::move(P_coeffs)));
    }

    const size_t n_chunks =
        (P_coeffs.size() + chunk_size - 1) / chunk_size; // ceil
    const size_t m = (n_chunks <= 1) ? 0 : libff::log2(n_chunks - 1) + 1; // ceil(log2(n_chunks))

    std::vector<FieldT> remainder_coeffs(P_coeffs.begin(),
                                         P_coeffs.begin() + chunk_size);

    // Preon: quotient_size >= n_chunks * chunk_size, then copy dividend[|H|..]
    std::vector<FieldT> quotient_coeffs(n_chunks * chunk_size, FieldT::zero());
    std::copy(P_coeffs.begin() + chunk_size, P_coeffs.end(),
              quotient_coeffs.begin());

    for (size_t r = 0; r < m; ++r) {
        const size_t s = size_t(1) << r;
        for (size_t offset = 0; offset < n_chunks; offset += 2 * s) {
            if (offset + s > n_chunks - 1) {
                break;
            }
            const size_t src_chunk_index = (offset + s - 1) * chunk_size;
            for (size_t j = 0; j < chunk_size; ++j) {
                const FieldT src = quotient_coeffs[src_chunk_index + j];
                for (size_t i = 0; i < s - 1; ++i) {
                    quotient_coeffs[(offset + i) * chunk_size + j] += src;
                }
            }
        }
    }

    return std::make_pair(polynomial<FieldT>(std::move(quotient_coeffs)),
                          polynomial<FieldT>(std::move(remainder_coeffs)));
}

    /** Evaluate Z_H at a point for an m-dimensional Cantor/LCH linear subspace.
     *  For the Cantor basis, W_i(beta_i) = 1, so the subspace vanishing
     *  polynomials satisfy W_0(x) = x and W_{i+1}(x) = W_i(x)^2 + W_i(x).
     *  Thus Z_H(x) = W_m(x) where |H| = 2^m = Z.degree().
     *  Affine (nonzero) shifts are not handled here. */
    template <typename FieldT>
    FieldT lch_vanishing_polynomial_evaluation_at_point(const FieldT &evalpoint,
                                                        const vanishing_polynomial<FieldT> &Z)
    {
        if (!Z.is_cantor_basis())
        {
            throw std::invalid_argument(
                "lch_vanishing_polynomial_evaluation_at_point is only implemented for the Cantor basis");
        }

        const size_t deg = Z.degree();
        if (deg == 0)
        {
            return FieldT::zero();
        }
        if ((deg & (deg - 1)) != 0)
        {
            throw std::invalid_argument(
                "lch_vanishing_polynomial_evaluation_at_point: degree of the vanishing_polynomial must be a power of two");
        }

        const size_t dim = libff::log2(deg);
        FieldT result = evalpoint;
        for (size_t i = 0; i < dim; ++i)
        {
            result = result.squared() + result;
        }
        return result;
    }

} // namespace libiop
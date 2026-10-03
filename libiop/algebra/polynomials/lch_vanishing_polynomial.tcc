#include <cstddef>
#include <vector>
#include "libiop/algebra/polynomials/polynomial.hpp"
#include "libiop/algebra/polynomials/poly_basis.hpp"
#include "libiop/algebra/polynomials/lch_standard.hpp"

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
    const size_t vp_dim = libff::log2(chunk_size);
    const bool is_cantor = Z.is_cantor_basis();

    std::vector<FieldT> remainder_coeffs(P_coeffs.begin(),
                                         P_coeffs.begin() + chunk_size);

    std::vector<FieldT> quotient_coeffs(P_coeffs.begin() + chunk_size, P_coeffs.end());
    quotient_coeffs.resize(n_chunks * chunk_size);

    /** Round r folds W_m^{2^r} = W_{m+r} + (lower terms) back into the LCH basis.
     *  On Cantor z_i = 1 and the source-chunk scale is 1. On the standard basis
     *  the scale is (∏_{k<r} z_{m+k}) / z_{m+r} (fast_preon lch_poly_div_by_vp). */
    FieldT coeff_numinator = FieldT::one();
    for (size_t r = 0; r < m; ++r) {
        const FieldT z = lch_z_i<FieldT>(vp_dim + r, is_cantor);
        const FieldT coeff = is_cantor ? FieldT::one() : (coeff_numinator * z.inverse());
        const size_t s = size_t(1) << r;
        for (size_t offset = 0; offset < n_chunks; offset += 2 * s) {
            if (offset + s > n_chunks - 1) {
                break;
            }
            FieldT *src = &quotient_coeffs[(offset + s - 1) * chunk_size];
            if (!is_cantor)
            {
                for (size_t j = 0; j < chunk_size; ++j) {
                    src[j] *= coeff;
                }
            }
            for (size_t i = 0; i < s - 1; ++i) {
                FieldT *dst = &quotient_coeffs[(offset + i) * chunk_size];
                for (size_t j = 0; j < chunk_size; ++j) {
                    dst[j] += src[j];
                }
            }
        }
        if (!is_cantor)
        {
            coeff_numinator *= z;
        }
    }

    return std::make_pair(polynomial<FieldT>(std::move(quotient_coeffs)),
                          polynomial<FieldT>(std::move(remainder_coeffs)));
}

    /** Evaluate Z_H at a point for an m-dimensional linear subspace.
     *  Z_H = W_m, W_0(x)=x, W_{i+1}=W_i^2 + z_i W_i with z_i = W_i(β_i)
     *  (z_i = 1 on Cantor). Affine (nonzero) shifts of H itself are not handled here. */
    template <typename FieldT>
    FieldT lch_vanishing_polynomial_evaluation_at_point(const FieldT &evalpoint,
                                                        const vanishing_polynomial<FieldT> &Z)
    {
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

        return lch_linear_vanishing_at_point(evalpoint, libff::log2(deg), Z.is_cantor_basis());
    }

} // namespace libiop
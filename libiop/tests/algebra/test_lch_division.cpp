// Division by subspace vanishing polynomials in the LCH basis (Algorithm 3):
// checks A = Z_m * Q + R at random points and deg(R) < 2^m, on the standard and
// Cantor bases, for dividend lengths that are and are not multiples of 2^m.
#include <algorithm>
#include <cstddef>
#include <gtest/gtest.h>
#include <vector>

#include <libff/algebra/fields/binary/gf64.hpp>
#include <libff/algebra/fields/binary/gf128.hpp>
#include <libff/algebra/fields/binary/gf192.hpp>
#include <libff/algebra/fields/binary/gf256.hpp>
#include "libiop/algebra/field_subset/field_subset.hpp"
#include "libiop/algebra/polynomials/polynomial.hpp"
#include "libiop/algebra/polynomials/vanishing_polynomial.hpp"
#include "libiop/algebra/polynomials/poly_basis.hpp"
#include "libiop/algebra/polynomials/lch_standard.hpp"
#include "libiop/algebra/polynomials/lch_polynomial.hpp"
#include "libiop/algebra/polynomials/lch_vanishing_polynomial.hpp"
#include "libiop/algebra/utils.hpp"

namespace libiop {

// Stage-wise schedule (each source chunk added to every chunk below it in its
// group), used only as a reference: Algorithm 3 must return identical quotients.
template<typename FieldT>
std::vector<FieldT> stagewise_quotient(const std::vector<FieldT> &P, const std::size_t chunk, const bool is_cantor)
{
    const std::size_t n_chunks = (P.size() + chunk - 1) / chunk;
    const std::size_t m = (n_chunks <= 1) ? 0 : libff::log2(n_chunks - 1) + 1;
    const std::size_t vp_dim = libff::log2(chunk);
    std::vector<FieldT> q(P.begin() + chunk, P.end());
    q.resize(n_chunks * chunk);
    FieldT num = FieldT::one();
    for (std::size_t r = 0; r < m; ++r) {
        const FieldT z = lch_z_i<FieldT>(vp_dim + r, is_cantor);
        const FieldT c = is_cantor ? FieldT::one() : num * z.inverse();
        const std::size_t s = std::size_t(1) << r;
        for (std::size_t off = 0; off < n_chunks; off += 2 * s) {
            if (off + s > n_chunks - 1) break;
            FieldT *src = &q[(off + s - 1) * chunk];
            for (std::size_t j = 0; j < chunk; ++j) src[j] *= c;
            for (std::size_t i = 0; i + 1 < s; ++i)
                for (std::size_t j = 0; j < chunk; ++j) q[(off + i) * chunk + j] += src[j];
        }
        if (!is_cantor) num *= z;
    }
    return q;
}

template<typename FieldT>
void run_lch_division_test(const bool is_cantor)
{
    set_polynomial_basis_config(polynomial_basis_config::lch_poly_basis);
    set_lch_evaluation_is_cantor(is_cantor);
    for (std::size_t vp_dim = 2; vp_dim <= 6; ++vp_dim)  // dims 0,1 are ambiguous: {1} is both bases
    {
        const std::size_t chunk = std::size_t(1) << vp_dim;
        const field_subset<FieldT> H(chunk, is_cantor);
        const vanishing_polynomial<FieldT> Z(H);
        for (std::size_t len = chunk + 1; len <= 40 * chunk + 7; len += (len < 4 * chunk ? 1 : 3))
        {
            const polynomial<FieldT> A(random_vector<FieldT>(len));
            const auto QR = lch_basis_polynomial_over_vanishing_polynomial(A, Z);
            const std::vector<FieldT> &Q = QR.first.coefficients();
            const std::vector<FieldT> &R = QR.second.coefficients();
            ASSERT_LE(R.size(), chunk);
            std::vector<FieldT> expected = stagewise_quotient(A.coefficients(), chunk, is_cantor);
            std::vector<FieldT> got = Q;
            const std::size_t width = std::max(expected.size(), got.size());
            expected.resize(width, FieldT::zero());
            got.resize(width, FieldT::zero());
            ASSERT_TRUE(got == expected)
                << "vp_dim=" << vp_dim << " len=" << len << " cantor=" << is_cantor;
            for (std::size_t trial = 0; trial < 3; ++trial)
            {
                const FieldT x = FieldT::random_element();
                const FieldT lhs = lch_polynomial_evaluation_at_point(A.coefficients(), x);
                const FieldT rhs =
                    lch_vanishing_polynomial_evaluation_at_point(x, Z) * lch_polynomial_evaluation_at_point(Q, x)
                    + lch_polynomial_evaluation_at_point(R, x);
                ASSERT_TRUE(lhs == rhs) << "vp_dim=" << vp_dim << " len=" << len << " cantor=" << is_cantor;
            }
        }
    }
}

TEST(LCHDivisionTest, StandardBasis) {
    run_lch_division_test<libff::gf256>(false);
}

TEST(LCHDivisionTest, CantorBasis) {
    run_lch_division_test<libff::gf256>(true);
}

}

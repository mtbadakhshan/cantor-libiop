#include <cstddef>
#include <vector>

#include "depends/additive-fft/C++/LCH/fft.hpp"
#include "libiop/algebra/polynomials/lch_standard.hpp"


namespace libiop
{
    template<typename FieldT>
    FieldT lch_polynomial_evaluation_at_point(const std::vector<FieldT> &lch_coeffs,
                                            const FieldT &x)
    {
        if (lch_coeffs.empty()) {
            return FieldT::zero();
        }

        size_t n = 1;
        while (n < lch_coeffs.size()) {
            n <<= 1;
        }
        const size_t d = libff::log2(n);

        /* W_r(x) for r < d; d <= 64 since n fits in a size_t. */
        FieldT X[64];
        lch_normalized_W_chain(X, d, x, get_lch_evaluation_is_cantor());

        /** Fold f = f_lo + W_r(x) * f_hi level by level. Working copy padded to n
         *  (FieldT() is zero); the top level only touches the live prefix, so the
         *  padding never needs a separate zero pass. */
        std::vector<FieldT> c(lch_coeffs);
        c.resize(n);

        for (size_t r = d; r-- > 0;) {
            const size_t half = size_t(1) << r;
            const FieldT &Xr = X[r];
            for (size_t j = 0; j < half; ++j) {
                c[j] += Xr * c[j + half];
            }
        }
        return c[0];
    }

    template<typename FieldT>
    std::vector<FieldT> lch_fri_fold(const std::vector<FieldT> &lch_coeffs,
                                     const std::size_t eta,
                                     const FieldT &x_i,
                                     const FieldT &bit0_scale)
    {
        const size_t k = size_t(1) << eta;
        std::vector<FieldT> Xt(k, FieldT::one());
        if (eta > 0)
        {
            FieldT W[64];
            lch_normalized_W_chain(W, eta, x_i, get_lch_evaluation_is_cantor());
            W[0] *= bit0_scale;
            for (size_t t = 1; t < k; ++t)
            {
                Xt[t] = Xt[t & (t - 1)] * W[__builtin_ctzll(t)];
            }
        }

        const size_t out_len = (lch_coeffs.size() + k - 1) / k;
        std::vector<FieldT> out(out_len, FieldT::zero());
        for (size_t j = 0; j < out_len; ++j)
        {
            for (size_t t = 0; t < k && j * k + t < lch_coeffs.size(); ++t)
            {
                out[j] += Xt[t] * lch_coeffs[j * k + t];
            }
        }
        return out;
    }

    template<typename FieldT>
    std::vector<FieldT> lch_fri_fold(const std::vector<FieldT> &lch_coeffs,
                                     const std::size_t eta,
                                     const FieldT &x_i,
                                     const lch_basis_tables<FieldT> &tables,
                                     const FieldT &bit0_scale)
    {
        const size_t k = size_t(1) << eta;
        std::vector<FieldT> Xt(k, FieldT::one());
        if (eta > 0)
        {
            FieldT W[64];
            lch_normalized_W_chain_from_tables(W, eta, x_i, tables);
            W[0] *= bit0_scale;
            for (size_t t = 1; t < k; ++t)
            {
                Xt[t] = Xt[t & (t - 1)] * W[__builtin_ctzll(t)];
            }
        }

        const size_t out_len = (lch_coeffs.size() + k - 1) / k;
        std::vector<FieldT> out(out_len, FieldT::zero());
        for (size_t j = 0; j < out_len; ++j)
        {
            for (size_t t = 0; t < k && j * k + t < lch_coeffs.size(); ++t)
            {
                out[j] += Xt[t] * lch_coeffs[j * k + t];
            }
        }
        return out;
    }

    template<typename FieldT>
    std::vector<FieldT> lch_to_monomial(const std::vector<FieldT> &lch_coeffs)
    {
        if (lch_coeffs.empty())
        {
            return {};
        }
        size_t n = 1;
        while (n < lch_coeffs.size())
        {
            n <<= 1;
        }
        std::vector<FieldT> mono(lch_coeffs);
        mono.resize(n, FieldT::zero());
        lch::inv_basis_conversion(mono, n);
        return mono;
    }

} // namespace libiop
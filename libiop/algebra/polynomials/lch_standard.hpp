#ifndef LIBIOP_ALGEBRA_POLYNOMIALS_LCH_STANDARD_HPP_
#define LIBIOP_ALGEBRA_POLYNOMIALS_LCH_STANDARD_HPP_

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "libiop/algebra/polynomials/poly_basis.hpp"

namespace libiop
{

/** Normalized Ŵ_r(β_j) and unnormalized z_r = W_r(β_r) for the standard
 *  ordered basis β_j = 2^j. Matches fast_preon's (w_std, w_std_not_normalized).
 *  On Cantor every z_r is 1 and the W-recurrence needs no table. */
template<typename FieldT>
struct lch_standard_tables
{
    static constexpr std::size_t max_dim = 32;
    FieldT w[max_dim][max_dim];
    FieldT z[max_dim];
};

template<typename FieldT>
const lch_standard_tables<FieldT> &get_lch_standard_tables()
{
    thread_local lch_standard_tables<FieldT> tables;
    thread_local bool initialized = false;
    if (initialized)
    {
        return tables;
    }

    constexpr std::size_t D = lch_standard_tables<FieldT>::max_dim;
    std::vector<FieldT> w_nn_prev(D);
    for (std::size_t j = 0; j < D; ++j)
    {
        tables.w[0][j] = FieldT(1ull << j);
        w_nn_prev[j] = tables.w[0][j];
    }
    tables.z[0] = tables.w[0][0];

    std::vector<FieldT> w_nn_cur(D);
    for (std::size_t i = 1; i < D; ++i)
    {
        for (std::size_t j = 0; j < i; ++j)
        {
            tables.w[i][j] = FieldT::zero();
            w_nn_cur[j] = FieldT::zero();
        }
        for (std::size_t j = i; j < D; ++j)
        {
            w_nn_cur[j] = w_nn_prev[j].squared() + w_nn_prev[j] * w_nn_prev[i - 1];
        }
        tables.z[i] = w_nn_cur[i];
        const FieldT z_inv = tables.z[i].inverse();
        tables.w[i][i] = FieldT::one();
        for (std::size_t j = i + 1; j < D; ++j)
        {
            tables.w[i][j] = w_nn_cur[j] * z_inv;
        }
        w_nn_prev.swap(w_nn_cur);
    }

    initialized = true;
    return tables;
}

template<typename FieldT>
FieldT lch_z_i(const std::size_t i, const bool is_cantor_basis)
{
    if (is_cantor_basis)
    {
        return FieldT::one();
    }
    if (i >= lch_standard_tables<FieldT>::max_dim)
    {
        throw std::invalid_argument("lch_z_i: index exceeds the standard-basis table");
    }
    return get_lch_standard_tables<FieldT>().z[i];
}

/** Normalized W_r(x) for r < d of the active evaluation basis. */
template<typename FieldT>
void lch_normalized_W_chain(FieldT *X, const std::size_t d, const FieldT &x,
                            const bool is_cantor_basis)
{
    if (d == 0)
    {
        return;
    }
    X[0] = x;
    if (is_cantor_basis)
    {
        for (std::size_t r = 1; r < d; ++r)
        {
            X[r] = X[r - 1].squared() + X[r - 1];
        }
        return;
    }
    const auto &w = get_lch_standard_tables<FieldT>().w;
    if (d > lch_standard_tables<FieldT>::max_dim)
    {
        throw std::invalid_argument("lch_normalized_W_chain: degree exceeds the standard-basis table");
    }
    for (std::size_t r = 1; r < d; ++r)
    {
        const FieldT fe = w[r - 1][r].squared() + w[r - 1][r] * w[r - 1][r - 1];
        X[r] = (X[r - 1].squared() + w[r - 1][r - 1] * X[r - 1]) * fe.inverse();
    }
}

/** Unnormalized linearized vanishing polynomial of an m-dimensional linear
 *  subspace: W_0(x)=x, W_{i+1}=W_i^2 + z_i W_i. */
template<typename FieldT>
FieldT lch_linear_vanishing_at_point(const FieldT &evalpoint,
                                     const std::size_t dim,
                                     const bool is_cantor_basis)
{
    FieldT result = evalpoint;
    for (std::size_t i = 0; i < dim; ++i)
    {
        result = result.squared() + lch_z_i<FieldT>(i, is_cantor_basis) * result;
    }
    return result;
}

/** Twiddles at butterfly level r: tw[u] = Ŵ_r(shift + sum_k u_k β_{r+1+k}).
 *  shift_dim is a bitmask of standard-basis indices (β_i = 2^i). */
template<typename FieldT>
void lch_standard_twiddles(std::vector<FieldT> &tw,
                           const std::size_t r,
                           const std::size_t dimension,
                           const std::size_t shift_dim)
{
    constexpr std::size_t max_dim = lch_standard_tables<FieldT>::max_dim;
    const auto &w = get_lch_standard_tables<FieldT>().w;
    tw[0] = FieldT::zero();
    for (std::size_t k = 0, bits = shift_dim; bits != 0; ++k, bits >>= 1)
    {
        if ((bits & 1) && k < max_dim)
        {
            tw[0] += w[r][k];
        }
    }
    for (std::size_t k = 0; (r + 1) + k < dimension; ++k)
    {
        const FieldT shift = w[r][(r + 1) + k];
        const std::size_t half = std::size_t(1) << k;
        for (std::size_t u = 0; u < half; ++u)
        {
            tw[u + half] = tw[u] + shift;
        }
    }
}

/** LCH butterfly on a standard ordered basis. shift_dim is a bitmask of the shift. */
template<typename FieldT>
std::vector<FieldT> lch_standard_btfly(const std::vector<FieldT> &poly_coeffs,
                                       const std::size_t domain_dim,
                                       const std::size_t shift_dim)
{
    const std::size_t n = std::size_t(1) << domain_dim;
    if (poly_coeffs.size() > n)
    {
        throw std::invalid_argument("lch_standard_btfly: more coefficients than domain elements");
    }
    std::size_t d = 1;
    while (d < poly_coeffs.size())
    {
        d <<= 1;
    }
    const unsigned log_d = __builtin_ctzll(d);

    std::vector<FieldT> g(n);
    for (std::size_t block = 0; block < n; block += d)
    {
        std::copy(poly_coeffs.begin(), poly_coeffs.end(), g.begin() + block);
    }

    const std::size_t live_upper =
        poly_coeffs.size() > (d >> 1) ? poly_coeffs.size() - (d >> 1) : 0;
    std::vector<FieldT> tw(n);
    for (unsigned i = log_d; i > 0; i--)
    {
        const unsigned unit = (1u << i);
        const unsigned num = n / unit;
        const unsigned unit_2 = unit / 2;
        lch_standard_twiddles(tw, i - 1, domain_dim, shift_dim);
        const bool sparse_top = (i == log_d && live_upper < unit_2);
        for (unsigned j = 0; j < num; j++)
        {
            const unsigned n_live = sparse_top ? (unsigned)live_upper : unit_2;
            const unsigned off = j * unit;
            const FieldT &twu = tw[j];
            for (unsigned k = 0; k < n_live; k++)
            {
                g[off + k] += g[off + unit_2 + k] * twu;
                g[off + unit_2 + k] += g[off + k];
            }
            if (sparse_top)
            {
                for (unsigned k = n_live; k < unit_2; k++)
                {
                    g[off + unit_2 + k] = g[off + k];
                }
            }
        }
    }
    return g;
}

template<typename FieldT>
std::vector<FieldT> lch_standard_ibtfly(const std::vector<FieldT> &evals,
                                        const std::size_t domain_dim,
                                        const std::size_t shift_dim)
{
    const std::size_t n = std::size_t(1) << domain_dim;
    std::vector<FieldT> g(evals);
    g.resize(n, FieldT::zero());
    std::vector<FieldT> tw(n);
    for (unsigned i = 1; i <= domain_dim; i++)
    {
        const unsigned unit = (1u << i);
        const unsigned num = n / unit;
        const unsigned unit_2 = unit / 2;
        lch_standard_twiddles(tw, i - 1, domain_dim, shift_dim);
        for (unsigned j = 0; j < num; j++)
        {
            const unsigned off = j * unit;
            const FieldT &twu = tw[j];
            for (unsigned k = 0; k < unit_2; k++)
            {
                g[off + unit_2 + k] += g[off + k];
                g[off + k] += g[off + unit_2 + k] * twu;
            }
        }
    }
    return g;
}

/** Ŵ_r(γ_j) and z_r = W_r(γ_r) for an arbitrary ordered basis γ.
 *  Same recurrence as lch_standard_tables; the first row is γ_j, not 2^j. */
template<typename FieldT>
struct lch_basis_tables
{
    static constexpr std::size_t max_dim = lch_standard_tables<FieldT>::max_dim;
    std::size_t dim = 0;
    FieldT w[max_dim][max_dim];
    FieldT z[max_dim];
};

template<typename FieldT>
void lch_fill_basis_tables(lch_basis_tables<FieldT> &tables,
                           const std::vector<FieldT> &basis)
{
    const std::size_t D = basis.size();
    if (D > lch_basis_tables<FieldT>::max_dim)
    {
        throw std::invalid_argument("lch_fill_basis_tables: basis exceeds max_dim");
    }
    tables.dim = D;
    if (D == 0)
    {
        return;
    }
    std::vector<FieldT> w_nn_prev(D);
    for (std::size_t j = 0; j < D; ++j)
    {
        tables.w[0][j] = basis[j];
        w_nn_prev[j] = basis[j];
    }
    tables.z[0] = basis[0];

    std::vector<FieldT> w_nn_cur(D);
    for (std::size_t i = 1; i < D; ++i)
    {
        for (std::size_t j = 0; j < i; ++j)
        {
            tables.w[i][j] = FieldT::zero();
            w_nn_cur[j] = FieldT::zero();
        }
        for (std::size_t j = i; j < D; ++j)
        {
            w_nn_cur[j] = w_nn_prev[j].squared() + w_nn_prev[j] * w_nn_prev[i - 1];
        }
        tables.z[i] = w_nn_cur[i];
        const FieldT z_inv = tables.z[i].inverse();
        tables.w[i][i] = FieldT::one();
        for (std::size_t j = i + 1; j < D; ++j)
        {
            tables.w[i][j] = w_nn_cur[j] * z_inv;
        }
        w_nn_prev.swap(w_nn_cur);
    }
}

template<typename FieldT>
void lch_normalized_W_chain_from_tables(FieldT *X, const std::size_t d, const FieldT &x,
                                        const lch_basis_tables<FieldT> &tables)
{
    if (d == 0)
    {
        return;
    }
    if (d > tables.dim)
    {
        throw std::invalid_argument("lch_normalized_W_chain_from_tables: degree exceeds the table");
    }
    X[0] = x;
    for (std::size_t r = 1; r < d; ++r)
    {
        const FieldT fe = tables.w[r - 1][r].squared() + tables.w[r - 1][r] * tables.w[r - 1][r - 1];
        X[r] = (X[r - 1].squared() + tables.w[r - 1][r - 1] * X[r - 1]) * fe.inverse();
    }
}

/** Twiddles Ŵ_r(shift + sum_k u_k γ_{r+1+k}). Wr_shift is Ŵ_r(shift). */
template<typename FieldT>
void lch_basis_twiddles(std::vector<FieldT> &tw,
                        const std::size_t r,
                        const std::size_t dimension,
                        const FieldT &Wr_shift,
                        const lch_basis_tables<FieldT> &tables)
{
    tw[0] = Wr_shift;
    for (std::size_t k = 0; (r + 1) + k < dimension; ++k)
    {
        const FieldT add = tables.w[r][(r + 1) + k];
        const std::size_t half = std::size_t(1) << k;
        for (std::size_t u = 0; u < half; ++u)
        {
            tw[u + half] = tw[u] + add;
        }
    }
}

/** LCH butterfly for a general ordered basis. The affine shift is a field
 *  element (later FRI images are not a standard-basis bitmask). */
template<typename FieldT>
std::vector<FieldT> lch_basis_btfly(const std::vector<FieldT> &poly_coeffs,
                                    const std::size_t domain_dim,
                                    const FieldT &shift,
                                    const lch_basis_tables<FieldT> &tables)
{
    const std::size_t n = std::size_t(1) << domain_dim;
    if (poly_coeffs.size() > n)
    {
        throw std::invalid_argument("lch_basis_btfly: more coefficients than domain elements");
    }
    if (domain_dim > tables.dim)
    {
        throw std::invalid_argument("lch_basis_btfly: domain_dim exceeds the table");
    }
    std::size_t d = 1;
    while (d < poly_coeffs.size())
    {
        d <<= 1;
    }
    const unsigned log_d = __builtin_ctzll(d);

    std::vector<FieldT> g(n);
    for (std::size_t block = 0; block < n; block += d)
    {
        std::copy(poly_coeffs.begin(), poly_coeffs.end(), g.begin() + block);
    }

    FieldT Wshift[lch_basis_tables<FieldT>::max_dim];
    lch_normalized_W_chain_from_tables(Wshift, domain_dim, shift, tables);

    const std::size_t live_upper =
        poly_coeffs.size() > (d >> 1) ? poly_coeffs.size() - (d >> 1) : 0;
    std::vector<FieldT> tw(n);
    for (unsigned i = log_d; i > 0; i--)
    {
        const unsigned unit = (1u << i);
        const unsigned num = n / unit;
        const unsigned unit_2 = unit / 2;
        lch_basis_twiddles(tw, i - 1, domain_dim, Wshift[i - 1], tables);
        /* Ŵ_r(γ_r)=1 for r≥1; Ŵ_0(γ_0)=γ_0, which is not 1 on a Z-image. */
        const FieldT wr = tables.w[i - 1][i - 1];
        const bool sparse_top = (i == log_d && live_upper < unit_2);
        for (unsigned j = 0; j < num; j++)
        {
            const unsigned n_live = sparse_top ? (unsigned)live_upper : unit_2;
            const unsigned off = j * unit;
            const FieldT &twu = tw[j];
            for (unsigned k = 0; k < n_live; k++)
            {
                g[off + k] += g[off + unit_2 + k] * twu;
                g[off + unit_2 + k] = g[off + k] + g[off + unit_2 + k] * wr;
            }
            if (sparse_top)
            {
                for (unsigned k = n_live; k < unit_2; k++)
                {
                    g[off + unit_2 + k] = g[off + k];
                }
            }
        }
    }
    return g;
}

template<typename FieldT>
std::vector<FieldT> lch_basis_ibtfly(const std::vector<FieldT> &evals,
                                     const std::size_t domain_dim,
                                     const FieldT &shift,
                                     const lch_basis_tables<FieldT> &tables)
{
    const std::size_t n = std::size_t(1) << domain_dim;
    if (domain_dim > tables.dim)
    {
        throw std::invalid_argument("lch_basis_ibtfly: domain_dim exceeds the table");
    }
    std::vector<FieldT> g(evals);
    g.resize(n, FieldT::zero());
    FieldT Wshift[lch_basis_tables<FieldT>::max_dim];
    lch_normalized_W_chain_from_tables(Wshift, domain_dim, shift, tables);
    std::vector<FieldT> tw(n);
    for (unsigned i = 1; i <= domain_dim; i++)
    {
        const unsigned unit = (1u << i);
        const unsigned num = n / unit;
        const unsigned unit_2 = unit / 2;
        lch_basis_twiddles(tw, i - 1, domain_dim, Wshift[i - 1], tables);
        const FieldT wr_inv = tables.w[i - 1][i - 1].inverse();
        for (unsigned j = 0; j < num; j++)
        {
            const unsigned off = j * unit;
            const FieldT &twu = tw[j];
            for (unsigned k = 0; k < unit_2; k++)
            {
                g[off + unit_2 + k] = (g[off + unit_2 + k] + g[off + k]) * wr_inv;
                g[off + k] += g[off + unit_2 + k] * twu;
            }
        }
    }
    return g;
}

template<typename FieldT>
void lch_copy_standard_prefix(lch_basis_tables<FieldT> &tables, const std::size_t dim)
{
    if (dim > lch_basis_tables<FieldT>::max_dim)
    {
        throw std::invalid_argument("lch_copy_standard_prefix: dim exceeds max_dim");
    }
    const auto &src = get_lch_standard_tables<FieldT>();
    tables.dim = dim;
    for (std::size_t i = 0; i < dim; ++i)
    {
        tables.z[i] = src.z[i];
        for (std::size_t j = 0; j < dim; ++j)
        {
            tables.w[i][j] = src.w[i][j];
        }
    }
}

/** Tail of the current tower: Ŵ_{η+k}(β_{η+s}), for LCH in y = Ŵ_η. */
template<typename FieldT>
void lch_fill_tail_tables(lch_basis_tables<FieldT> &tail,
                          const lch_basis_tables<FieldT> &src,
                          const std::size_t eta)
{
    if (eta > src.dim)
    {
        throw std::invalid_argument("lch_fill_tail_tables: eta exceeds the source dimension");
    }
    tail.dim = src.dim - eta;
    for (std::size_t k = 0; k < tail.dim; ++k)
    {
        tail.z[k] = src.z[eta + k];
        for (std::size_t s = 0; s < tail.dim; ++s)
        {
            tail.w[k][s] = src.w[eta + k][eta + s];
        }
    }
}

/** Evaluate y-LCH of the fold at Ŵ_η(x) for x in the current domain, in
 *  L^{(i+1)} order: Butterfly on L^{(i+1)}/z_η with the tail table
 *  (first row is Ŵ_η(β_{η+s}), so the implicit basis is already scaled). */
template<typename FieldT>
std::vector<FieldT> lch_evals_of_y_lch_on_z_image(
    const std::vector<FieldT> &y_lch,
    const field_subset<FieldT> &next_domain,
    const std::size_t eta,
    const lch_basis_tables<FieldT> &current_tables)
{
    lch_basis_tables<FieldT> tail;
    lch_fill_tail_tables(tail, current_tables, eta);
    if (tail.dim != next_domain.dimension())
    {
        throw std::invalid_argument("lch_evals_of_y_lch_on_z_image: tail dim != next domain dim");
    }
    FieldT shift_y = next_domain.shift();
    if (eta > 0)
    {
        const FieldT z = current_tables.z[eta];
        if (z != FieldT::one())
        {
            shift_y *= z.inverse();
        }
    }
    return lch_basis_btfly<FieldT>(y_lch, next_domain.dimension(), shift_y, tail);
}

} // namespace libiop

#endif // LIBIOP_ALGEBRA_POLYNOMIALS_LCH_STANDARD_HPP_

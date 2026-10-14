#ifndef LIBIOP_PROTOCOLS_LDT_FRI_LOCAL_FOLD_HPP_
#define LIBIOP_PROTOCOLS_LDT_FRI_LOCAL_FOLD_HPP_

#include <limits>
#include <stdexcept>
#include <vector>
#include "libiop/algebra/field_subset/field_subset.hpp"

namespace libiop {

/** Prover arithmetic selection only; the folding map and verifier are unchanged.
 * configured preserves the existing polynomial-basis-dependent dispatch.
 * local_evaluations overrides it on additive domains. */
enum class fri_fold_implementation { configured, local_evaluations };

inline fri_fold_implementation &fri_fold_implementation_storage()
{
    thread_local fri_fold_implementation value = fri_fold_implementation::configured;
    return value;
}

inline void set_fri_fold_implementation(const fri_fold_implementation value)
{
    fri_fold_implementation_storage() = value;
}

inline fri_fold_implementation get_fri_fold_implementation()
{
    return fri_fold_implementation_storage();
}

/** Domain-only preprocessing for inversion-free local evaluation folding.
 * Each two-point fold is a + ((x+p)/gamma)*(a+b). For eta > 1,
 * both the domain and the SAME univariate challenge are mapped through
 * q_s(X) = X^2 + gamma_s*X at each stage. The result equals the
 * degree-<2^eta coset interpolant evaluated at the original challenge.
 * There is one inversion per binary stage during setup, none per coset.
 * Tables may be reused across codewords/challenges on the same domain. */
template<typename FieldT>
struct fri_local_fold_tables
{
    size_t eta = 0;
    std::vector<FieldT> pivot;      /* gamma_s, s < eta */
    std::vector<FieldT> pivot_inv;  /* 1/gamma_s */
    std::vector<FieldT> coset_c;    /* [j*eta + s] = c_j^{(s)} / gamma_s */
    std::vector<std::vector<FieldT>> intra_u; /* [s][t] = u_t^{(s)} / gamma_s, t < 2^{eta-s-1} */
};

/* Subset sums: out[j] = base + sum_{b} bit_b(j) gens[b]. Additions only. */
template<typename FieldT>
void fri_local_subset_sums(std::vector<FieldT> &out, const FieldT &base, const std::vector<FieldT> &gens)
{
    out.assign(size_t(1) << gens.size(), FieldT::zero());
    out[0] = base;
    for (size_t b = 0; b < gens.size(); ++b)
    {
        const size_t h = size_t(1) << b;
        for (size_t j = 0; j < h; ++j)
        {
            out[h + j] = out[j] + gens[b];
        }
    }
}

template<typename FieldT>
fri_local_fold_tables<FieldT> make_fri_local_fold_tables(const field_subset<FieldT> &dom, const size_t eta)
{
    if (dom.type() != affine_subspace_type || eta == 0 ||
        eta > dom.dimension() || dom.dimension() >= std::numeric_limits<size_t>::digits)
        throw std::invalid_argument("local FRI fold requires an additive domain and 1 <= eta <= dimension");
    fri_local_fold_tables<FieldT> T;
    T.eta = eta;
    const std::vector<FieldT> &basis = dom.basis();
    const size_t m = basis.size();
    const size_t num_cosets = size_t(1) << (m - eta);
    std::vector<FieldT> g(basis.begin(), basis.begin() + eta);   /* intra basis images */
    std::vector<FieldT> o(basis.begin() + eta, basis.end());     /* outer basis images */
    FieldT th = dom.shift();
    T.coset_c.assign(num_cosets * eta, FieldT::zero());
    T.intra_u.resize(eta);
    std::vector<FieldT> sums;
    for (size_t s = 0; s < eta; ++s)
    {
        const FieldT pv = g[s];
        const FieldT pinv = pv.inverse();
        T.pivot.push_back(pv);
        T.pivot_inv.push_back(pinv);
        std::vector<FieldT> o_scaled(o.size());
        for (size_t b = 0; b < o.size(); ++b) o_scaled[b] = o[b] * pinv;
        fri_local_subset_sums(sums, th * pinv, o_scaled);
        for (size_t j = 0; j < num_cosets; ++j) T.coset_c[j * eta + s] = sums[j];
        std::vector<FieldT> g_scaled;
        for (size_t b = s + 1; b < eta; ++b) g_scaled.push_back(g[b] * pinv);
        fri_local_subset_sums(T.intra_u[s], FieldT::zero(), g_scaled);
        /* apply q_s(X) = X^2 + gamma_s X */
        for (size_t b = s + 1; b < eta; ++b) g[b] = g[b].squared() + pv * g[b];
        for (FieldT &v : o) v = v.squared() + pv * v;
        th = th.squared() + pv * th;
    }
    return T;
}

template<typename FieldT>
std::vector<FieldT> evaluate_fri_local_fold(const std::vector<FieldT> &f, const fri_local_fold_tables<FieldT> &T, const FieldT &x)
{
    const size_t eta = T.eta;
    if (eta == 0 || eta >= std::numeric_limits<size_t>::digits)
        throw std::invalid_argument("invalid local FRI tables");
    const size_t k = size_t(1) << eta;
    const size_t num_cosets = f.size() / k;
    if (f.size() % k != 0 || num_cosets != T.coset_c.size() / eta)
        throw std::invalid_argument("local FRI input length does not match its tables");
    FieldT X[64];
    FieldT xs = x;
    for (size_t s = 0; s < eta; ++s)
    {
        X[s] = xs * T.pivot_inv[s];
        xs = xs.squared() + T.pivot[s] * xs;
    }
    std::vector<FieldT> out(num_cosets);
    std::vector<FieldT> g(k);
    for (size_t j = 0; j < num_cosets; ++j)
    {
        const FieldT *src = &f[j * k];
        const FieldT *cc = &T.coset_c[j * eta];
        /* step 0 reads the input directly */
        size_t half = k >> 1;
        {
            const FieldT base = X[0] + cc[0];
            const std::vector<FieldT> &U = T.intra_u[0];
            for (size_t t = 0; t < half; ++t)
            {
                const FieldT a = src[2 * t];
                g[t] = a + (base + U[t]) * (a + src[2 * t + 1]);
            }
        }
        for (size_t s = 1; s < eta; ++s)
        {
            half >>= 1;
            const FieldT base = X[s] + cc[s];
            const std::vector<FieldT> &U = T.intra_u[s];
            for (size_t t = 0; t < half; ++t)
            {
                const FieldT a = g[2 * t];
                g[t] = a + (base + U[t]) * (a + g[2 * t + 1]);
            }
        }
        out[j] = g[0];
    }
    return out;
}


} // namespace libiop

#endif

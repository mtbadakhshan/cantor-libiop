// FRI folding micro-benchmark over GF(2^256): the prover's fold loop only
// (no Merkle hashing, no query answering, no final interpolation).
//
// Methods
//   lagrange : libiop's evaluation-domain fold (per-coset Lagrange interpolation
//              with one batch inversion per coset), monomial configuration.
//   eval     : inversion-free evaluation-domain fold. Per coset of size 2^eta,
//              eta binary steps g <- g0 + (x_s + p)/gamma_s * (g0 + g1); the
//              per-coset constants p/gamma_s are domain-only preprocessing.
//   lch      : Algorithm 1 as in fri_ldt.tcc (lch_poly_basis): IButterfly on the
//              size-D prefix, coefficient mixing, Butterfly to the next codeword.
//
// One run = one call per round on a single codeword, as in FRI with one LDT
// instance and one interactive repetition.
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <libff/algebra/fields/binary/gf256.hpp>
#include <libff/common/profiling.hpp>

#include "libiop/algebra/fft.hpp"
#include "libiop/algebra/field_subset/field_subset.hpp"
#include "libiop/algebra/field_subset/subspace.hpp"
#include "libiop/algebra/polynomials/lch_polynomial.hpp"
#include "libiop/algebra/polynomials/lch_standard.hpp"
#include "libiop/algebra/polynomials/poly_basis.hpp"
#include "libiop/protocols/ldt/fri/fri_aux.hpp"
#include "libiop/protocols/ldt/fri/fri_local_fold.hpp"
#include "libiop/protocols/ldt/fri/fri_ldt.hpp"
#include "libiop/protocols/ldt/fri/localizer_polynomial.hpp"

using namespace libiop;
using FieldT = libff::gf256;
using clk = std::chrono::steady_clock;

static double ms_since(const clk::time_point &t0)
{
    return std::chrono::duration<double, std::milli>(clk::now() - t0).count();
}

/* Same as FRI_protocol::compute_domains (additive case). */
static std::vector<field_subset<FieldT>> compute_domains(const field_subset<FieldT> &L,
                                                         const std::vector<size_t> &etas)
{
    std::vector<field_subset<FieldT>> domains{L};
    for (size_t i = 0; i < etas.size(); ++i)
    {
        const std::vector<FieldT> &basis = domains[i].basis();
        const affine_subspace<FieldT> loc_space(
            std::vector<FieldT>(basis.begin(), basis.begin() + etas[i]), FieldT::zero());
        const localizer_polynomial<FieldT> loc(loc_space);
        std::vector<FieldT> next_basis(basis.begin() + etas[i], basis.end());
        for (FieldT &b : next_basis)
        {
            b = loc.evaluation_at_point(b);
        }
        domains.emplace_back(affine_subspace<FieldT>(next_basis, loc.evaluation_at_point(domains[i].shift())));
    }
    return domains;
}

/* ---------------- inversion-free evaluation fold ---------------- */

using eval_round_tables = fri_local_fold_tables<FieldT>;

static eval_round_tables eval_preprocess(const field_subset<FieldT> &dom, const size_t eta)
{
    return make_fri_local_fold_tables(dom, eta);
}

static std::vector<FieldT> eval_fold(const std::vector<FieldT> &f, const eval_round_tables &T, const FieldT &x)
{
    return evaluate_fri_local_fold(f, T, x);
}

/* ---------------- configuration ---------------- */

struct config
{
    bool cantor = false;
    size_t log_D = 12;
    size_t rs_extra = 5;      /* rho = 2^{-rs_extra} */
    size_t eta = 1;
    std::string schedule = "libiop";  /* libiop: [1, eta, eta, ...]; uniform: [eta, eta, ...] */
    int reps = 5;
    int warmups = 1;
    bool validate = false;
    std::vector<std::string> methods{"lagrange", "eval", "lch"};
    std::string csv;
};

static std::vector<size_t> schedule_for(const config &c)
{
    const size_t logn = c.log_D + c.rs_extra;
    if (c.schedule == "libiop")
    {
        return FRI_protocol_parameters<FieldT>::localization_parameter_to_array(c.eta, logn, c.rs_extra);
    }
    if (c.schedule == "uniform")
    {
        return std::vector<size_t>(c.log_D / c.eta, c.eta);
    }
    throw std::invalid_argument("schedule must be libiop or uniform");
}

struct run_result
{
    double total_ms = 0, pre_ms = 0, init_ms = 0;
    std::vector<double> round_ms;
    std::vector<std::vector<FieldT>> outputs; /* per round, only when validating */
};

static run_result run_lagrange(const std::vector<FieldT> &f0, const std::vector<field_subset<FieldT>> &doms,
                               const std::vector<size_t> &etas, const std::vector<FieldT> &xs, bool keep)
{
    set_polynomial_basis_config(polynomial_basis_config::monomial_poly_basis);
    run_result r;
    auto cur = std::make_shared<std::vector<FieldT>>(f0);
    const auto t_all = clk::now();
    for (size_t i = 0; i < etas.size(); ++i)
    {
        const auto t0 = clk::now();
        cur = evaluate_next_f_i_over_entire_domain(cur, doms[i], size_t(1) << etas[i], xs[i]);
        r.round_ms.push_back(ms_since(t0));
        if (keep) r.outputs.push_back(*cur);
    }
    r.total_ms = ms_since(t_all);
    return r;
}

static run_result run_eval(const std::vector<FieldT> &f0, const std::vector<field_subset<FieldT>> &doms,
                           const std::vector<size_t> &etas, const std::vector<FieldT> &xs, bool keep)
{
    run_result r;
    const auto tp = clk::now();
    std::vector<eval_round_tables> tables;
    for (size_t i = 0; i < etas.size(); ++i) tables.push_back(eval_preprocess(doms[i], etas[i]));
    r.pre_ms = ms_since(tp);
    std::vector<FieldT> cur = f0;
    const auto t_all = clk::now();
    for (size_t i = 0; i < etas.size(); ++i)
    {
        const auto t0 = clk::now();
        cur = eval_fold(cur, tables[i], xs[i]);
        r.round_ms.push_back(ms_since(t0));
        if (keep) r.outputs.push_back(cur);
    }
    r.total_ms = ms_since(t_all);
    return r;
}

/* Mirrors the lch_fold branch of FRI_protocol::calculate_and_submit_proof. */
static run_result run_lch(const std::vector<FieldT> &f0, const std::vector<field_subset<FieldT>> &doms,
                          const std::vector<size_t> &etas, const std::vector<FieldT> &xs,
                          const size_t D, bool keep)
{
    set_polynomial_basis_config(polynomial_basis_config::lch_poly_basis);
    run_result r;
    /* Domain-only tables (lch_copy_standard_prefix / lch_fill_basis_tables). */
    const auto tp = clk::now();
    std::vector<lch_basis_tables<FieldT>> tables(etas.size());
    std::vector<bool> is_cantor(etas.size()), is_standard(etas.size());
    for (size_t i = 0; i < etas.size(); ++i)
    {
        is_cantor[i] = doms[i].is_cantor_basis();
        is_standard[i] = doms[i].subspace().is_standard_basis();
        if (!is_cantor[i])
        {
            if (is_standard[i]) lch_copy_standard_prefix(tables[i], doms[i].dimension());
            else lch_fill_basis_tables(tables[i], doms[i].basis());
        }
    }
    r.pre_ms = ms_since(tp);

    const auto t_all = clk::now();
    const auto ti = clk::now();
    std::vector<FieldT> coeffs = IFFT_of_known_degree_over_field_subset<FieldT>(f0, D, doms[0]);
    r.init_ms = ms_since(ti);
    FieldT bit0_scale = FieldT::one();
    std::vector<FieldT> evals;
    for (size_t i = 0; i < etas.size(); ++i)
    {
        const auto t0 = clk::now();
        if (is_cantor[i] || is_standard[i])
            coeffs = lch_fri_fold(coeffs, etas[i], xs[i], bit0_scale);
        else
            coeffs = lch_fri_fold(coeffs, etas[i], xs[i], tables[i], bit0_scale);
        if (is_cantor[i])
            evals = FFT_over_field_subset<FieldT>(coeffs, doms[i + 1]);
        else
            evals = additive_FFT_wrapper<FieldT>(coeffs, doms[i + 1], etas[i], tables[i]);
        r.round_ms.push_back(ms_since(t0));
        if (keep) r.outputs.push_back(evals);
        bit0_scale = is_cantor[i] ? FieldT::one() : tables[i].z_inv[etas[i]];
    }
    r.total_ms = ms_since(t_all);
    return r;
}

static double median(std::vector<double> v)
{
    std::sort(v.begin(), v.end());
    const size_t n = v.size();
    return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

static std::vector<std::string> split(const std::string &s, char d)
{
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, d)) if (!item.empty()) out.push_back(item);
    return out;
}

int main(int argc, char **argv)
{
    config c;
    for (int a = 1; a < argc; ++a)
    {
        const std::string k = argv[a];
        auto next = [&]() -> std::string {
            if (a + 1 >= argc) throw std::invalid_argument("missing value for " + k);
            return argv[++a];
        };
        if (k == "--basis") { const std::string v = next(); c.cantor = (v == "cantor"); if (v != "cantor" && v != "standard") throw std::invalid_argument("basis"); }
        else if (k == "--log-D") c.log_D = std::stoul(next());
        else if (k == "--rs-extra") c.rs_extra = std::stoul(next());
        else if (k == "--eta") c.eta = std::stoul(next());
        else if (k == "--schedule") c.schedule = next();
        else if (k == "--reps") c.reps = std::stoi(next());
        else if (k == "--warmups") c.warmups = std::stoi(next());
        else if (k == "--methods") c.methods = split(next(), ',');
        else if (k == "--validate") c.validate = true;
        else if (k == "--csv") c.csv = next();
        else
        {
            std::cerr << "usage: bench_fri_fold --basis standard|cantor --log-D d --rs-extra r --eta e "
                         "[--schedule libiop|uniform] [--reps 5] [--warmups 1] "
                         "[--methods lagrange,eval,lch] [--validate] [--csv out.csv]\n";
            return 2;
        }
    }
    try
    {
        libff::inhibit_profiling_info = true;
        libff::inhibit_profiling_counters = true;
        set_additive_fft_cantor_implementation(additive_fft_cantor_implementation::lch_afft);
        set_lch_evaluation_is_cantor(c.cantor);

        const size_t D = size_t(1) << c.log_D;
        const size_t logn = c.log_D + c.rs_extra;
        const size_t n = size_t(1) << logn;
        const std::vector<size_t> etas = schedule_for(c);
        if (etas.empty()) throw std::invalid_argument("no rounds");

        /* Codeword domain as in aurora_iop: same-size subspace, shifted by an outside element. */
        const field_subset<FieldT> unshifted(n, c.cantor);
        const field_subset<FieldT> L(n, unshifted.element_outside_of_subset(), c.cantor);

        /* Honest codeword of degree < D (monomial coefficients, monomial FFT). */
        set_polynomial_basis_config(polynomial_basis_config::monomial_poly_basis);
        std::vector<FieldT> poly(D);
        for (FieldT &v : poly) v = FieldT::random_element();
        const std::vector<FieldT> f0 = FFT_over_field_subset<FieldT>(poly, L);
        std::vector<FieldT> xs;
        for (size_t i = 0; i < etas.size(); ++i) xs.push_back(FieldT::random_element());

        set_polynomial_basis_config(polynomial_basis_config::monomial_poly_basis);
        const auto doms_mono = compute_domains(L, etas);
        set_polynomial_basis_config(polynomial_basis_config::lch_poly_basis);
        const auto doms_lch = compute_domains(L, etas);

        if (c.validate)
        {
            const run_result ref = run_lagrange(f0, doms_mono, etas, xs, true);
            /* (ii) vs (i) on a random word (not a codeword). */
            std::vector<FieldT> rnd(n);
            for (FieldT &v : rnd) v = FieldT::random_element();
            const run_result ref_rnd = run_lagrange(rnd, doms_mono, etas, xs, true);
            const run_result ev_rnd = run_eval(rnd, doms_mono, etas, xs, true);
            const run_result ev = run_eval(f0, doms_mono, etas, xs, true);
            const run_result lc = run_lch(f0, doms_lch, etas, xs, D, true);
            bool ok_ev = true, ok_ev_rnd = true, ok_lch = true;
            for (size_t i = 0; i < etas.size(); ++i)
            {
                ok_ev_rnd &= ev_rnd.outputs[i] == ref_rnd.outputs[i];
                ok_ev &= ev.outputs[i] == ref.outputs[i];
                ok_lch &= lc.outputs[i] == ref.outputs[i];
            }
            std::cerr << "validate basis=" << (c.cantor ? "cantor" : "standard") << " logD=" << c.log_D
                      << " rs_extra=" << c.rs_extra << " eta=" << c.eta << " schedule=" << c.schedule
                      << " rounds=" << etas.size()
                      << " eval(random word)=" << (ok_ev_rnd ? "OK" : "FAIL")
                      << " eval(codeword)=" << (ok_ev ? "OK" : "FAIL")
                      << " lch(codeword)=" << (ok_lch ? "OK" : "FAIL") << "\n";
            if (!(ok_ev && ok_ev_rnd && ok_lch)) return 1;
        }

        std::ofstream out;
        bool header = false;
        if (!c.csv.empty())
        {
            std::ifstream probe(c.csv);
            header = !probe.good() || probe.peek() == std::ifstream::traits_type::eof();
            out.open(c.csv, std::ios::app);
            if (!out) throw std::runtime_error("cannot open csv");
        }
        if (header || c.csv.empty())
        {
            (c.csv.empty() ? std::cout : out)
                << "field,basis,method,schedule,log_D,log_n,rho,eta,num_rounds,reps,total_ms,preprocessing_ms,init_ms,round_ms\n";
        }
        std::string sched;
        for (size_t i = 0; i < etas.size(); ++i) sched += (i ? "-" : "") + std::to_string(etas[i]);

        for (const std::string &method : c.methods)
        {
            std::vector<double> tot, pre, init;
            std::vector<std::vector<double>> rounds(etas.size());
            for (int s = 0; s < c.warmups + c.reps; ++s)
            {
                run_result r;
                if (method == "lagrange") r = run_lagrange(f0, doms_mono, etas, xs, false);
                else if (method == "eval") r = run_eval(f0, doms_mono, etas, xs, false);
                else if (method == "lch") r = run_lch(f0, doms_lch, etas, xs, D, false);
                else throw std::invalid_argument("unknown method " + method);
                if (s < c.warmups) continue;
                tot.push_back(r.total_ms);
                pre.push_back(r.pre_ms);
                init.push_back(r.init_ms);
                for (size_t i = 0; i < etas.size(); ++i) rounds[i].push_back(r.round_ms[i]);
            }
            std::ostringstream row;
            row << "gf256," << (c.cantor ? "cantor" : "standard") << ',' << method << ',' << c.schedule
                << ':' << sched << ',' << c.log_D << ',' << logn << ",1/" << (size_t(1) << c.rs_extra) << ','
                << c.eta << ',' << etas.size() << ',' << c.reps << ',' << median(tot) << ',' << median(pre) << ','
                << median(init) << ',';
            for (size_t i = 0; i < etas.size(); ++i) row << (i ? ";" : "") << median(rounds[i]);
            row << '\n';
            (c.csv.empty() ? std::cout : out) << row.str();
            if (!c.csv.empty()) std::cerr << row.str();
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}

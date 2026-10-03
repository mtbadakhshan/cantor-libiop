#include <vector>
#include <benchmark/benchmark.h>
#include <chrono>


// #include <libff/algebra/fields/binary/gf265.hpp>
#include <libff/algebra/curves/bls12_381/bls12_381_pp.hpp>
#include "libiop/snark/aurora_snark.hpp"
#include "libiop/relations/examples/r1cs_examples.hpp"
#include "libiop/bcs/hashing/hashing.hpp"
#include "libiop/algebra/fft.hpp"

/* Cantor-basis rows: set_additive_fft_cantor_implementation(lch_afft|cantor_afft).
 * LCH polynomial basis is a separate switch, set_polynomial_basis_config(lch_poly_basis).
 * It works with Cantor (is_cantor_basis=true) or standard (false) evaluation domains.
 * See FFT_benchmark_column_mapping.txt. */

namespace libiop {

typedef binary_hash_digest hash_type;
const LDT_reducer_soundness_type ldt_reducer_soundness_type = LDT_reducer_soundness_type::proven;
const FRI_soundness_type fri_soundness_type = FRI_soundness_type::proven;
const std::size_t num_inputs = (1 << 5) - 1;
const size_t security_parameter = 128;
const size_t RS_extra_dimensions = 5;
const size_t FRI_localization_parameter = 1;
const bool make_zk = true;

// ######## PRIME FIELD ########

static void BM_AURORA_PROVER_PRIME_FIELD(benchmark::State &state)
{
    libff::bls12_381_pp::init_public_params();
    typedef libff::Fr<libff::bls12_381_pp> FieldT;
    const field_subset_type domain_type = multiplicative_coset_type;
    const size_t sz = state.range(0);
    const std::size_t num_constraints = 1 << sz;
    const std::size_t num_variables = (1 << sz) - 1;

    r1cs_example<FieldT> r1cs_params = generate_r1cs_example<FieldT>(
        num_constraints, num_inputs, num_variables);

    aurora_snark_parameters<FieldT, hash_type> params( security_parameter, ldt_reducer_soundness_type,
            fri_soundness_type, blake2b_type, FRI_localization_parameter, RS_extra_dimensions,
            make_zk, domain_type,  num_constraints, num_variables);
    aurora_snark_argument<FieldT, hash_type> argument;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(r1cs_params);
        benchmark::DoNotOptimize(params);
        benchmark::DoNotOptimize(argument = aurora_snark_prover<FieldT>(r1cs_params.constraint_system_, 
                                    r1cs_params.primary_input_, r1cs_params.auxiliary_input_, params));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
    state.counters["codeword_dim"] = params.iop_params_.codeword_domain_dim();
}


// ######## BINARY EXTENSION FIELD ########

static void BM_AURORA_PROVER_CANTOR_BASIS_LCH(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const field_subset_type domain_type = affine_subspace_type;
    const bool is_cantor_basis = true;
    const size_t sz = state.range(0);
    const std::size_t num_constraints = 1 << sz;
    const std::size_t num_variables = (1 << sz) - 1;

    r1cs_example<FieldT> r1cs_params = generate_r1cs_example<FieldT>(
        num_constraints, num_inputs, num_variables);

    aurora_snark_parameters<FieldT, hash_type> params( security_parameter, ldt_reducer_soundness_type,
            fri_soundness_type, blake2b_type, FRI_localization_parameter, RS_extra_dimensions,
            make_zk, domain_type, is_cantor_basis,  num_constraints, num_variables);
    aurora_snark_argument<FieldT, hash_type> argument;
    set_polynomial_basis_config(polynomial_basis_config::monomial_poly_basis);
    set_additive_fft_cantor_implementation(additive_fft_cantor_implementation::lch_afft);
    state.counters["poly_basis"] = static_cast<double>(polynomial_basis_config::monomial_poly_basis);
    state.counters["cantor_fft_impl"] = static_cast<double>(additive_fft_cantor_implementation::lch_afft);
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(r1cs_params);
        benchmark::DoNotOptimize(params);
        benchmark::DoNotOptimize(argument = aurora_snark_prover<FieldT>(r1cs_params.constraint_system_, 
                                    r1cs_params.primary_input_, r1cs_params.auxiliary_input_, params));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
    state.counters["codeword_dim"] = params.iop_params_.codeword_domain_dim();
}

static void BM_AURORA_PROVER_CANTOR_BASIS_CANTOR(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const field_subset_type domain_type = affine_subspace_type;
    const bool is_cantor_basis = true;
    const size_t sz = state.range(0);
    const std::size_t num_constraints = 1 << sz;
    const std::size_t num_variables = (1 << sz) - 1;

    r1cs_example<FieldT> r1cs_params = generate_r1cs_example<FieldT>(
        num_constraints, num_inputs, num_variables);

    aurora_snark_parameters<FieldT, hash_type> params( security_parameter, ldt_reducer_soundness_type,
            fri_soundness_type, blake2b_type, FRI_localization_parameter, RS_extra_dimensions,
            make_zk, domain_type, is_cantor_basis,  num_constraints, num_variables);
    aurora_snark_argument<FieldT, hash_type> argument;
    set_polynomial_basis_config(polynomial_basis_config::monomial_poly_basis);
    set_additive_fft_cantor_implementation(additive_fft_cantor_implementation::cantor_afft);
    state.counters["poly_basis"] = static_cast<double>(polynomial_basis_config::monomial_poly_basis);
    state.counters["cantor_fft_impl"] = static_cast<double>(additive_fft_cantor_implementation::cantor_afft);
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(r1cs_params);
        benchmark::DoNotOptimize(params);
        benchmark::DoNotOptimize(argument = aurora_snark_prover<FieldT>(r1cs_params.constraint_system_, 
                                    r1cs_params.primary_input_, r1cs_params.auxiliary_input_, params));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
    state.counters["codeword_dim"] = params.iop_params_.codeword_domain_dim();
}


static void BM_AURORA_VERIFIER_CANTOR_BASIS_LCH(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const field_subset_type domain_type = affine_subspace_type;
    const bool is_cantor_basis = true;
    const size_t sz = state.range(0);
    const std::size_t num_constraints = 1 << sz;
    const std::size_t num_variables = (1 << sz) - 1;

    r1cs_example<FieldT> r1cs_params = generate_r1cs_example<FieldT>(
        num_constraints, num_inputs, num_variables);

    aurora_snark_parameters<FieldT, hash_type> params( security_parameter, ldt_reducer_soundness_type,
            fri_soundness_type, blake2b_type, FRI_localization_parameter, RS_extra_dimensions,
            make_zk, domain_type, is_cantor_basis,  num_constraints, num_variables);
    set_polynomial_basis_config(polynomial_basis_config::monomial_poly_basis);
    set_additive_fft_cantor_implementation(additive_fft_cantor_implementation::lch_afft);
    aurora_snark_argument<FieldT, hash_type> argument = aurora_snark_prover<FieldT>(r1cs_params.constraint_system_, 
                                    r1cs_params.primary_input_, r1cs_params.auxiliary_input_, params);
    state.counters["poly_basis"] = static_cast<double>(polynomial_basis_config::monomial_poly_basis);
    state.counters["cantor_fft_impl"] = static_cast<double>(additive_fft_cantor_implementation::lch_afft);
    bool bit;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(r1cs_params);
        benchmark::DoNotOptimize(params);
        benchmark::DoNotOptimize(bit = aurora_snark_verifier<FieldT, hash_type>( r1cs_params.constraint_system_,
                                                            r1cs_params.primary_input_, argument, params));
        benchmark::ClobberMemory();
        
        if (!bit){
            throw std::invalid_argument("Verification Failed");
        }

    }
    state.SetItemsProcessed(state.iterations());
    state.counters["codeword_dim"] = params.iop_params_.codeword_domain_dim();
}

static void BM_AURORA_VERIFIER_CANTOR_BASIS_CANTOR(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const field_subset_type domain_type = affine_subspace_type;
    const bool is_cantor_basis = true;
    const size_t sz = state.range(0);
    const std::size_t num_constraints = 1 << sz;
    const std::size_t num_variables = (1 << sz) - 1;

    r1cs_example<FieldT> r1cs_params = generate_r1cs_example<FieldT>(
        num_constraints, num_inputs, num_variables);

    aurora_snark_parameters<FieldT, hash_type> params( security_parameter, ldt_reducer_soundness_type,
            fri_soundness_type, blake2b_type, FRI_localization_parameter, RS_extra_dimensions,
            make_zk, domain_type, is_cantor_basis,  num_constraints, num_variables);
    set_polynomial_basis_config(polynomial_basis_config::monomial_poly_basis);
    set_additive_fft_cantor_implementation(additive_fft_cantor_implementation::cantor_afft);
    aurora_snark_argument<FieldT, hash_type> argument = aurora_snark_prover<FieldT>(r1cs_params.constraint_system_, 
                                    r1cs_params.primary_input_, r1cs_params.auxiliary_input_, params);
    state.counters["poly_basis"] = static_cast<double>(polynomial_basis_config::monomial_poly_basis);
    state.counters["cantor_fft_impl"] = static_cast<double>(additive_fft_cantor_implementation::cantor_afft);
    bool bit;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(r1cs_params);
        benchmark::DoNotOptimize(params);
        benchmark::DoNotOptimize(bit = aurora_snark_verifier<FieldT, hash_type>( r1cs_params.constraint_system_,
                                                            r1cs_params.primary_input_, argument, params));
        benchmark::ClobberMemory();
        
        if (!bit){
            throw std::invalid_argument("Verification Failed");
        }

    }
    state.SetItemsProcessed(state.iterations());
    state.counters["codeword_dim"] = params.iop_params_.codeword_domain_dim();
}


static void BM_AURORA_PROVER_STANDARD_BASIS(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const field_subset_type domain_type = affine_subspace_type;
    const bool is_cantor_basis = false;
    const size_t sz = state.range(0);
    const std::size_t num_constraints = 1 << sz;
    const std::size_t num_variables = (1 << sz) - 1;

    r1cs_example<FieldT> r1cs_params = generate_r1cs_example<FieldT>(
        num_constraints, num_inputs, num_variables);

    aurora_snark_parameters<FieldT, hash_type> params( security_parameter, ldt_reducer_soundness_type,
            fri_soundness_type, blake2b_type, FRI_localization_parameter, RS_extra_dimensions,
            make_zk, domain_type, is_cantor_basis,  num_constraints, num_variables);
    set_polynomial_basis_config(polynomial_basis_config::monomial_poly_basis);
    state.counters["poly_basis"] = static_cast<double>(polynomial_basis_config::monomial_poly_basis);
    aurora_snark_argument<FieldT, hash_type> argument;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(r1cs_params);
        benchmark::DoNotOptimize(params);
        benchmark::DoNotOptimize(argument = aurora_snark_prover<FieldT>(r1cs_params.constraint_system_, 
                                    r1cs_params.primary_input_, r1cs_params.auxiliary_input_, params));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
    state.counters["codeword_dim"] = params.iop_params_.codeword_domain_dim();
}

static void BM_AURORA_VERIFIER_STANDARD_BASIS(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const field_subset_type domain_type = affine_subspace_type;
    const bool is_cantor_basis = false;
    const size_t sz = state.range(0);
    const std::size_t num_constraints = 1 << sz;
    const std::size_t num_variables = (1 << sz) - 1;

    r1cs_example<FieldT> r1cs_params = generate_r1cs_example<FieldT>(
        num_constraints, num_inputs, num_variables);

    aurora_snark_parameters<FieldT, hash_type> params( security_parameter, ldt_reducer_soundness_type,
            fri_soundness_type, blake2b_type, FRI_localization_parameter, RS_extra_dimensions,
            make_zk, domain_type, is_cantor_basis,  num_constraints, num_variables);
    set_polynomial_basis_config(polynomial_basis_config::monomial_poly_basis);
    state.counters["poly_basis"] = static_cast<double>(polynomial_basis_config::monomial_poly_basis);
    aurora_snark_argument<FieldT, hash_type> argument = aurora_snark_prover<FieldT>(r1cs_params.constraint_system_, 
                                    r1cs_params.primary_input_, r1cs_params.auxiliary_input_, params);
    bool bit;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(r1cs_params);
        benchmark::DoNotOptimize(params);
        benchmark::DoNotOptimize(bit = aurora_snark_verifier<FieldT, hash_type>( r1cs_params.constraint_system_,
                                                            r1cs_params.primary_input_, argument, params));
        benchmark::ClobberMemory();
        
        if (!bit){
            throw std::invalid_argument("Verification Failed");
        }

    }
    state.SetItemsProcessed(state.iterations());
    state.counters["codeword_dim"] = params.iop_params_.codeword_domain_dim();
}

/** Cantor evaluation basis, polynomials stored in the LCH basis.
 *  FFT/IFFT use only the butterfly (no monomial basis conversion).
 *  The additive-FFT implementation switch is not consulted on this path. */
static void BM_AURORA_PROVER_CANTOR_BASIS_LCH_POLY_BASIS(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const field_subset_type domain_type = affine_subspace_type;
    const bool is_cantor_basis = true;
    const size_t sz = state.range(0);
    const std::size_t num_constraints = 1 << sz;
    const std::size_t num_variables = (1 << sz) - 1;

    r1cs_example<FieldT> r1cs_params = generate_r1cs_example<FieldT>(
        num_constraints, num_inputs, num_variables);

    aurora_snark_parameters<FieldT, hash_type> params( security_parameter, ldt_reducer_soundness_type,
            fri_soundness_type, blake2b_type, FRI_localization_parameter, RS_extra_dimensions,
            make_zk, domain_type, is_cantor_basis,  num_constraints, num_variables);
    aurora_snark_argument<FieldT, hash_type> argument;
    set_polynomial_basis_config(polynomial_basis_config::lch_poly_basis);
    state.counters["poly_basis"] = static_cast<double>(polynomial_basis_config::lch_poly_basis);
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(r1cs_params);
        benchmark::DoNotOptimize(params);
        benchmark::DoNotOptimize(argument = aurora_snark_prover<FieldT>(r1cs_params.constraint_system_, 
                                    r1cs_params.primary_input_, r1cs_params.auxiliary_input_, params));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
    state.counters["codeword_dim"] = params.iop_params_.codeword_domain_dim();
}

static void BM_AURORA_VERIFIER_CANTOR_BASIS_LCH_POLY_BASIS(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const field_subset_type domain_type = affine_subspace_type;
    const bool is_cantor_basis = true;
    const size_t sz = state.range(0);
    const std::size_t num_constraints = 1 << sz;
    const std::size_t num_variables = (1 << sz) - 1;

    r1cs_example<FieldT> r1cs_params = generate_r1cs_example<FieldT>(
        num_constraints, num_inputs, num_variables);

    aurora_snark_parameters<FieldT, hash_type> params( security_parameter, ldt_reducer_soundness_type,
            fri_soundness_type, blake2b_type, FRI_localization_parameter, RS_extra_dimensions,
            make_zk, domain_type, is_cantor_basis,  num_constraints, num_variables);
    set_polynomial_basis_config(polynomial_basis_config::lch_poly_basis);
    aurora_snark_argument<FieldT, hash_type> argument = aurora_snark_prover<FieldT>(r1cs_params.constraint_system_, 
                                    r1cs_params.primary_input_, r1cs_params.auxiliary_input_, params);
    state.counters["poly_basis"] = static_cast<double>(polynomial_basis_config::lch_poly_basis);
    bool bit;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(r1cs_params);
        benchmark::DoNotOptimize(params);
        benchmark::DoNotOptimize(bit = aurora_snark_verifier<FieldT, hash_type>( r1cs_params.constraint_system_,
                                                            r1cs_params.primary_input_, argument, params));
        benchmark::ClobberMemory();
        
        if (!bit){
            throw std::invalid_argument("Verification Failed");
        }

    }
    state.SetItemsProcessed(state.iterations());
    state.counters["codeword_dim"] = params.iop_params_.codeword_domain_dim();
}

/** Standard evaluation basis, polynomials stored in the LCH basis.
 *  FFT/IFFT use the LCH butterfly with standard-basis twiddles. */
static void BM_AURORA_PROVER_STANDARD_BASIS_LCH_POLY_BASIS(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const field_subset_type domain_type = affine_subspace_type;
    const bool is_cantor_basis = false;
    const size_t sz = state.range(0);
    const std::size_t num_constraints = 1 << sz;
    const std::size_t num_variables = (1 << sz) - 1;

    r1cs_example<FieldT> r1cs_params = generate_r1cs_example<FieldT>(
        num_constraints, num_inputs, num_variables);

    aurora_snark_parameters<FieldT, hash_type> params( security_parameter, ldt_reducer_soundness_type,
            fri_soundness_type, blake2b_type, FRI_localization_parameter, RS_extra_dimensions,
            make_zk, domain_type, is_cantor_basis,  num_constraints, num_variables);
    aurora_snark_argument<FieldT, hash_type> argument;
    set_polynomial_basis_config(polynomial_basis_config::lch_poly_basis);
    state.counters["poly_basis"] = static_cast<double>(polynomial_basis_config::lch_poly_basis);
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(r1cs_params);
        benchmark::DoNotOptimize(params);
        benchmark::DoNotOptimize(argument = aurora_snark_prover<FieldT>(r1cs_params.constraint_system_,
                                    r1cs_params.primary_input_, r1cs_params.auxiliary_input_, params));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
    state.counters["codeword_dim"] = params.iop_params_.codeword_domain_dim();
}

static void BM_AURORA_VERIFIER_STANDARD_BASIS_LCH_POLY_BASIS(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const field_subset_type domain_type = affine_subspace_type;
    const bool is_cantor_basis = false;
    const size_t sz = state.range(0);
    const std::size_t num_constraints = 1 << sz;
    const std::size_t num_variables = (1 << sz) - 1;

    r1cs_example<FieldT> r1cs_params = generate_r1cs_example<FieldT>(
        num_constraints, num_inputs, num_variables);

    aurora_snark_parameters<FieldT, hash_type> params( security_parameter, ldt_reducer_soundness_type,
            fri_soundness_type, blake2b_type, FRI_localization_parameter, RS_extra_dimensions,
            make_zk, domain_type, is_cantor_basis,  num_constraints, num_variables);
    set_polynomial_basis_config(polynomial_basis_config::lch_poly_basis);
    aurora_snark_argument<FieldT, hash_type> argument = aurora_snark_prover<FieldT>(r1cs_params.constraint_system_,
                                    r1cs_params.primary_input_, r1cs_params.auxiliary_input_, params);
    state.counters["poly_basis"] = static_cast<double>(polynomial_basis_config::lch_poly_basis);
    bool bit;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(r1cs_params);
        benchmark::DoNotOptimize(params);
        benchmark::DoNotOptimize(bit = aurora_snark_verifier<FieldT, hash_type>( r1cs_params.constraint_system_,
                                                            r1cs_params.primary_input_, argument, params));
        benchmark::ClobberMemory();

        if (!bit){
            throw std::invalid_argument("Verification Failed");
        }

    }
    state.SetItemsProcessed(state.iterations());
    state.counters["codeword_dim"] = params.iop_params_.codeword_domain_dim();
}

// BENCHMARK(BM_AURORA_PROVER_PRIME_FIELD)->DenseRange(5, 18, 1)->Unit(benchmark::kMillisecond)->ReportAggregatesOnly(true);
BENCHMARK(BM_AURORA_PROVER_CANTOR_BASIS_LCH)->DenseRange(9, 19, 1)->Unit(benchmark::kMillisecond)->ReportAggregatesOnly(true);
BENCHMARK(BM_AURORA_PROVER_CANTOR_BASIS_CANTOR)->DenseRange(9, 19, 1)->Unit(benchmark::kMillisecond)->ReportAggregatesOnly(true);
BENCHMARK(BM_AURORA_VERIFIER_CANTOR_BASIS_LCH)->DenseRange(9, 19, 1)->Unit(benchmark::kMillisecond)->ReportAggregatesOnly(true);
BENCHMARK(BM_AURORA_VERIFIER_CANTOR_BASIS_CANTOR)->DenseRange(9, 19, 1)->Unit(benchmark::kMillisecond)->ReportAggregatesOnly(true);
BENCHMARK(BM_AURORA_PROVER_CANTOR_BASIS_LCH_POLY_BASIS)->DenseRange(9, 19, 1)->Unit(benchmark::kMillisecond)->ReportAggregatesOnly(true);
BENCHMARK(BM_AURORA_VERIFIER_CANTOR_BASIS_LCH_POLY_BASIS)->DenseRange(9, 19, 1)->Unit(benchmark::kMillisecond)->ReportAggregatesOnly(true);
BENCHMARK(BM_AURORA_PROVER_STANDARD_BASIS)->DenseRange(9, 16, 1)->Unit(benchmark::kMillisecond)->ReportAggregatesOnly(true);
BENCHMARK(BM_AURORA_VERIFIER_STANDARD_BASIS)->DenseRange(9, 16, 1)->Unit(benchmark::kMillisecond)->ReportAggregatesOnly(true);
BENCHMARK(BM_AURORA_PROVER_STANDARD_BASIS_LCH_POLY_BASIS)->DenseRange(9, 16, 1)->Unit(benchmark::kMillisecond)->ReportAggregatesOnly(true);
BENCHMARK(BM_AURORA_VERIFIER_STANDARD_BASIS_LCH_POLY_BASIS)->DenseRange(9, 16, 1)->Unit(benchmark::kMillisecond)->ReportAggregatesOnly(true);


}

BENCHMARK_MAIN();

// Cantor-basis rows (LCH FFT, cantor:: FFT, and LCH polynomial basis):
//   ./benchmark_aurora --benchmark_filter='BM_AURORA_(PROVER|VERIFIER)_CANTOR_BASIS_.*'
// Standard evaluation basis + LCH polynomial basis:
//   ./benchmark_aurora --benchmark_filter='BM_AURORA_(PROVER|VERIFIER)_STANDARD_BASIS_LCH_POLY_BASIS'
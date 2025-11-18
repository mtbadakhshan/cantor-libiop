#include <vector>
#include <benchmark/benchmark.h>
#include <chrono>


// #include <libff/algebra/fields/binary/gf265.hpp>
#include <libff/algebra/curves/bls12_381/bls12_381_pp.hpp>
#include "libiop/snark/ligero_snark.hpp"
#include "libiop/bcs/common_bcs_parameters.hpp"
#include "libiop/relations/examples/r1cs_examples.hpp"
#include "libiop/bcs/hashing/hashing.hpp"


namespace libiop {

typedef binary_hash_digest hash_type;
const LDT_reducer_soundness_type ldt_reducer_soundness_type = LDT_reducer_soundness_type::proven;
const std::size_t num_inputs = (1 << 5) - 1;
const size_t security_parameter = 128;
const size_t RS_extra_dimensions = 5;
const double height_width_ratio = 1.0;
const bool make_zk = true;

// ######## PRIME FIELD ########

static void BM_LIGERO_PROVER_PRIME_FIELD(benchmark::State &state)
{ //Todo: Not done yet
    libff::bls12_381_pp::init_public_params();
    typedef libff::Fr<libff::bls12_381_pp> FieldT;
    const field_subset_type domain_type = multiplicative_coset_type;
    const size_t sz = state.range(0);
    std::size_t constraint_dim = sz;
    const std::size_t num_constraints = 1 << sz;
    const std::size_t num_variables = (1 << sz) - 1;

    r1cs_example<FieldT> r1cs_params = generate_r1cs_example<FieldT>(
        num_constraints, num_inputs, num_variables);

        ligero_snark_parameters<FieldT, hash_type> parameters;
        parameters.security_level_ = security_parameter;
        parameters.height_width_ratio_ = height_width_ratio;
        parameters.RS_extra_dimensions_ = RS_extra_dimensions;
        parameters.make_zk_ = make_zk;
        parameters.domain_type_ = domain_type;
        parameters.LDT_reducer_soundness_type_ = LDT_reducer_soundness_type::proven;
        parameters.bcs_params_ = default_bcs_params<FieldT, binary_hash_digest>(
            blake2b_type, parameters.security_level_, constraint_dim);
    ligero_snark_argument<FieldT, hash_type> argument;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(r1cs_params);
        benchmark::DoNotOptimize(parameters);
        benchmark::DoNotOptimize(argument = ligero_snark_prover<FieldT>(r1cs_params.constraint_system_, 
                                    r1cs_params.primary_input_, r1cs_params.auxiliary_input_, parameters));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
    const bool bit = ligero_snark_verifier<FieldT, binary_hash_digest>(r1cs_params.constraint_system_, r1cs_params.primary_input_, argument, parameters);
    if (!bit){
        //throw exception
        throw std::invalid_argument("Verification Failed");
    }
    ligero_iop_parameters<FieldT> iop_params = 
                        obtain_iop_parameters_from_ligero_snark_params(parameters, num_constraints, num_variables);

    state.counters["systematic_domain_dim"] = iop_params.systematic_domain_dim();
    state.counters["codeword_domain_dim"] = iop_params.systematic_domain_dim() + iop_params.RS_extra_dimensions();
}


// ######## BINARY EXTENSION FIELD ########

static void BM_LIGERO_PROVER_CANTOR_BASIS(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const field_subset_type domain_type = affine_subspace_type;
    const bool is_cantor_basis = true;
    const size_t sz = state.range(0);
    std::size_t constraint_dim = sz;
    const std::size_t num_constraints = 1 << sz;
    const std::size_t num_variables = (1 << sz) - 1;

    r1cs_example<FieldT> r1cs_params = generate_r1cs_example<FieldT>(
        num_constraints, num_inputs, num_variables);
    
    ligero_snark_parameters<FieldT, hash_type> parameters;
        parameters.security_level_ = security_parameter;
        parameters.height_width_ratio_ = height_width_ratio;
        parameters.RS_extra_dimensions_ = RS_extra_dimensions;
        parameters.make_zk_ = make_zk;
        parameters.domain_type_ = affine_subspace_type;
        parameters.is_cantor_basis_ = is_cantor_basis;
        parameters.LDT_reducer_soundness_type_ = LDT_reducer_soundness_type::proven;
        parameters.bcs_params_ = default_bcs_params<FieldT, binary_hash_digest>(
            blake2b_type, parameters.security_level_, constraint_dim);

    ligero_snark_argument<FieldT, hash_type> argument;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(r1cs_params);
        benchmark::DoNotOptimize(parameters);
        benchmark::DoNotOptimize(argument = ligero_snark_prover<FieldT>(r1cs_params.constraint_system_, 
                                    r1cs_params.primary_input_, r1cs_params.auxiliary_input_, parameters));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
    const bool bit = ligero_snark_verifier<FieldT, binary_hash_digest>(r1cs_params.constraint_system_, r1cs_params.primary_input_, argument, parameters);
    if (!bit){
        //throw exception
        throw std::invalid_argument("Verification Failed");
    }
    ligero_iop_parameters<FieldT> iop_params = 
                        obtain_iop_parameters_from_ligero_snark_params(parameters, num_constraints, num_variables);

    state.counters["systematic_domain_dim"] = iop_params.systematic_domain_dim();
    state.counters["codeword_domain_dim"] = iop_params.systematic_domain_dim() + iop_params.RS_extra_dimensions();
}

static void BM_LIGERO_VERIFIER_CANTOR_BASIS(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const field_subset_type domain_type = affine_subspace_type;
    const bool is_cantor_basis = true;
    const size_t sz = state.range(0);
    std::size_t constraint_dim = sz;
    const std::size_t num_constraints = 1 << sz;
    const std::size_t num_variables = (1 << sz) - 1;

    r1cs_example<FieldT> r1cs_params = generate_r1cs_example<FieldT>(
        num_constraints, num_inputs, num_variables);
    
    ligero_snark_parameters<FieldT, hash_type> parameters;
        parameters.security_level_ = security_parameter;
        parameters.height_width_ratio_ = height_width_ratio;
        parameters.RS_extra_dimensions_ = RS_extra_dimensions;
        parameters.make_zk_ = make_zk;
        parameters.domain_type_ = affine_subspace_type;
        parameters.is_cantor_basis_ = is_cantor_basis;
        parameters.LDT_reducer_soundness_type_ = LDT_reducer_soundness_type::proven;
        parameters.bcs_params_ = default_bcs_params<FieldT, binary_hash_digest>(
            blake2b_type, parameters.security_level_, constraint_dim);

    ligero_snark_argument<FieldT, hash_type> argument = ligero_snark_prover<FieldT>(r1cs_params.constraint_system_, 
                                    r1cs_params.primary_input_, r1cs_params.auxiliary_input_, parameters);
    bool bit;
                                    
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(r1cs_params);
        benchmark::DoNotOptimize(parameters);
        benchmark::DoNotOptimize(bit = ligero_snark_verifier<FieldT, binary_hash_digest>(r1cs_params.constraint_system_, r1cs_params.primary_input_, argument, parameters));
        benchmark::ClobberMemory();
        if (!bit){
            throw std::invalid_argument("Verification Failed");
        }
    }
    state.SetItemsProcessed(state.iterations());
    ligero_iop_parameters<FieldT> iop_params = 
                        obtain_iop_parameters_from_ligero_snark_params(parameters, num_constraints, num_variables);
    state.counters["systematic_domain_dim"] = iop_params.systematic_domain_dim();
    state.counters["codeword_domain_dim"] = iop_params.systematic_domain_dim() + iop_params.RS_extra_dimensions();
}


static void BM_LIGERO_PROVER_STANDARD_BASIS(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const field_subset_type domain_type = affine_subspace_type;
    const bool is_cantor_basis = false; //false for standard basis
    const size_t sz = state.range(0);
    std::size_t constraint_dim = sz;
    const std::size_t num_constraints = 1 << sz;
    const std::size_t num_variables = (1 << sz) - 1;

    r1cs_example<FieldT> r1cs_params = generate_r1cs_example<FieldT>(
        num_constraints, num_inputs, num_variables);
    
    ligero_snark_parameters<FieldT, hash_type> parameters;
        parameters.security_level_ = security_parameter;
        parameters.height_width_ratio_ = height_width_ratio;
        parameters.RS_extra_dimensions_ = RS_extra_dimensions;
        parameters.make_zk_ = make_zk;
        parameters.domain_type_ = affine_subspace_type;
        parameters.is_cantor_basis_ = is_cantor_basis;
        parameters.LDT_reducer_soundness_type_ = LDT_reducer_soundness_type::proven;
        parameters.bcs_params_ = default_bcs_params<FieldT, binary_hash_digest>(
            blake2b_type, parameters.security_level_, constraint_dim);

    ligero_snark_argument<FieldT, hash_type> argument;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(r1cs_params);
        benchmark::DoNotOptimize(parameters);
        benchmark::DoNotOptimize(argument = ligero_snark_prover<FieldT>(r1cs_params.constraint_system_, 
                                    r1cs_params.primary_input_, r1cs_params.auxiliary_input_, parameters));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
    const bool bit = ligero_snark_verifier<FieldT, binary_hash_digest>(r1cs_params.constraint_system_, r1cs_params.primary_input_, argument, parameters);
    if (!bit){
        //throw exception
        throw std::invalid_argument("Verification Failed");
    }
    ligero_iop_parameters<FieldT> iop_params = 
                        obtain_iop_parameters_from_ligero_snark_params(parameters, num_constraints, num_variables);

    state.counters["systematic_domain_dim"] = iop_params.systematic_domain_dim();
    state.counters["codeword_domain_dim"] = iop_params.systematic_domain_dim() + iop_params.RS_extra_dimensions();
}

static void BM_LIGERO_VERIFIER_STANDARD_BASIS(benchmark::State &state)
{
    typedef libff::gf256 FieldT;
    const field_subset_type domain_type = affine_subspace_type;
    const bool is_cantor_basis = false; //false for standard basis
    const size_t sz = state.range(0);
    std::size_t constraint_dim = sz;
    const std::size_t num_constraints = 1 << sz;
    const std::size_t num_variables = (1 << sz) - 1;

    r1cs_example<FieldT> r1cs_params = generate_r1cs_example<FieldT>(
        num_constraints, num_inputs, num_variables);
    
    ligero_snark_parameters<FieldT, hash_type> parameters;
        parameters.security_level_ = security_parameter;
        parameters.height_width_ratio_ = height_width_ratio;
        parameters.RS_extra_dimensions_ = RS_extra_dimensions;
        parameters.make_zk_ = make_zk;
        parameters.domain_type_ = affine_subspace_type;
        parameters.is_cantor_basis_ = is_cantor_basis;
        parameters.LDT_reducer_soundness_type_ = LDT_reducer_soundness_type::proven;
        parameters.bcs_params_ = default_bcs_params<FieldT, binary_hash_digest>(
            blake2b_type, parameters.security_level_, constraint_dim);

    ligero_snark_argument<FieldT, hash_type> argument = ligero_snark_prover<FieldT>(r1cs_params.constraint_system_, 
                                    r1cs_params.primary_input_, r1cs_params.auxiliary_input_, parameters);
    bool bit;
                                    
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(r1cs_params);
        benchmark::DoNotOptimize(parameters);
        benchmark::DoNotOptimize(bit = ligero_snark_verifier<FieldT, binary_hash_digest>(r1cs_params.constraint_system_, r1cs_params.primary_input_, argument, parameters));
        benchmark::ClobberMemory();
        if (!bit){
            throw std::invalid_argument("Verification Failed");
        }
    }
    state.SetItemsProcessed(state.iterations());
    ligero_iop_parameters<FieldT> iop_params = 
                        obtain_iop_parameters_from_ligero_snark_params(parameters, num_constraints, num_variables);
    state.counters["systematic_domain_dim"] = iop_params.systematic_domain_dim();
    state.counters["codeword_domain_dim"] = iop_params.systematic_domain_dim() + iop_params.RS_extra_dimensions();
}


BENCHMARK(BM_LIGERO_PROVER_PRIME_FIELD)->DenseRange(9, 20, 1)->Unit(benchmark::kMillisecond)->ReportAggregatesOnly(true);
BENCHMARK(BM_LIGERO_PROVER_CANTOR_BASIS)->DenseRange(9, 20, 1)->Unit(benchmark::kMillisecond)->ReportAggregatesOnly(true);
BENCHMARK(BM_LIGERO_VERIFIER_CANTOR_BASIS)->DenseRange(9, 20, 1)->Unit(benchmark::kMillisecond)->ReportAggregatesOnly(true);
BENCHMARK(BM_LIGERO_PROVER_STANDARD_BASIS)->DenseRange(9, 20, 1)->Unit(benchmark::kMillisecond)->ReportAggregatesOnly(true);
BENCHMARK(BM_LIGERO_VERIFIER_STANDARD_BASIS)->DenseRange(9, 20, 1)->Unit(benchmark::kMillisecond)->ReportAggregatesOnly(true);


}

BENCHMARK_MAIN();

//./libiop/benchmark_ligero --benchmark_repetitions=100 --benchmark_out=cantor_gm.json --benchmark_out_format=json
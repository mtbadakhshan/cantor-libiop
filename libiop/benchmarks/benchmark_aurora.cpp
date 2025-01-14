#include <vector>
#include <benchmark/benchmark.h>


// #include <libff/algebra/fields/binary/gf265.hpp>
#include "libiop/snark/aurora_snark.hpp"
#include "libiop/relations/examples/r1cs_examples.hpp"
#include "libiop/bcs/hashing/hashing.hpp"


namespace libiop {
typedef libff::gf256 FieldT;
typedef binary_hash_digest hash_type;
const LDT_reducer_soundness_type ldt_reducer_soundness_type = LDT_reducer_soundness_type::optimistic_heuristic;
const FRI_soundness_type fri_soundness_type = FRI_soundness_type::heuristic;
const field_subset_type domain_type = affine_subspace_type;
const std::size_t num_inputs = (1 << 5) - 1;
const size_t security_parameter = 128;
const size_t RS_extra_dimensions = 2;
const size_t FRI_localization_parameter = 3;
const bool make_zk = true;


static void BM_AURORA_PROVER_CANTOR_BASIS(benchmark::State &state)
{

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
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(r1cs_params);
        benchmark::DoNotOptimize(params);
        benchmark::DoNotOptimize(argument = aurora_snark_prover<FieldT>(r1cs_params.constraint_system_, 
                                    r1cs_params.primary_input_, r1cs_params.auxiliary_input_, params));
        benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations());
}

BENCHMARK(BM_AURORA_PROVER_CANTOR_BASIS)->DenseRange(17, 20, 1)->Unit(benchmark::kMilisecond)->ReportAggregatesOnly(true);

}

BENCHMARK_MAIN();

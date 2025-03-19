#include <cstdint>

#include <gtest/gtest.h>

#include <libff/algebra/fields/binary/gf64.hpp>
#include <libff/algebra/curves/edwards/edwards_pp.hpp>
#include <libff/algebra/curves/bn128/bn128_pp.hpp>
#include <libff/algebra/curves/bls12_381/bls12_381_pp.hpp>
#include "libiop/snark/aurora_snark.hpp"
#include "libiop/relations/examples/r1cs_examples.hpp"

namespace libiop {

TEST(AuroraSnarkTest, SimpleTest) {
    /* Set up R1CS */
    typedef libff::gf256 FieldT;
    typedef binary_hash_digest hash_type;

    const std::size_t num_constraints = 1 << 15;
    const std::size_t num_inputs = (1 << 5) - 1;
    const std::size_t num_variables = (1 << 15) - 1;
    const size_t security_parameter = 128;
    const size_t RS_extra_dimensions = 5;
    const size_t FRI_localization_parameter = 1;
    const LDT_reducer_soundness_type ldt_reducer_soundness_type = LDT_reducer_soundness_type::proven;
    const FRI_soundness_type fri_soundness_type = FRI_soundness_type::proven;
    const field_subset_type domain_type = affine_subspace_type;
    const bool is_cantor_basis = true;

    r1cs_example<FieldT> r1cs_params = generate_r1cs_example<FieldT>(
        num_constraints, num_inputs, num_variables);
    EXPECT_TRUE(r1cs_params.constraint_system_.is_satisfied(
        r1cs_params.primary_input_, r1cs_params.auxiliary_input_));

    /* Actual SNARK test */
    for (std::size_t i = 1; i < 2; i++) {
        const bool make_zk = (i == 0) ? false : true;
        aurora_snark_parameters<FieldT, hash_type> params(
            security_parameter,
            ldt_reducer_soundness_type,
            fri_soundness_type,
            blake2b_type,
            FRI_localization_parameter,
            RS_extra_dimensions,
            make_zk,
            domain_type,
            is_cantor_basis, 
            num_constraints,
            num_variables);
        const aurora_snark_argument<FieldT, hash_type> argument = aurora_snark_prover<FieldT>(
            r1cs_params.constraint_system_,
            r1cs_params.primary_input_,
            r1cs_params.auxiliary_input_,
            params);

        printf("iop size in bytes %lu\n", argument.IOP_size_in_bytes());
        printf("bcs size in bytes %lu\n", argument.BCS_size_in_bytes());
        printf("argument size in bytes %lu\n", argument.size_in_bytes());

        const bool bit = aurora_snark_verifier<FieldT, hash_type>(
            r1cs_params.constraint_system_,
            r1cs_params.primary_input_,
            argument,
            params);

        EXPECT_TRUE(bit) << "failed on make_zk = " << i << " test";
    }
}

TEST(AuroraSnarkMultiplicativeTest, SimpleTest) {
    /* Set up R1CS */
    libff::bls12_381_pp::init_public_params();
    typedef libff::Fr<libff::bls12_381_pp> FieldT;
    typedef binary_hash_digest hash_type;

    const size_t num_constraints = 1 << 15;
    const size_t num_inputs = (1 << 5) - 1;
    const size_t num_variables = (1 << 15) - 1;
    const size_t security_parameter = 128;
    const size_t RS_extra_dimensions = 5;
    const size_t FRI_localization_parameter = 1;
    const LDT_reducer_soundness_type ldt_reducer_soundness_type = LDT_reducer_soundness_type::optimistic_heuristic;
    const FRI_soundness_type fri_soundness_type = FRI_soundness_type::heuristic;
    const field_subset_type domain_type = multiplicative_coset_type;

    r1cs_example<FieldT> r1cs_params = generate_r1cs_example<FieldT>(
        num_constraints, num_inputs, num_variables);
    EXPECT_TRUE(r1cs_params.constraint_system_.is_satisfied(
        r1cs_params.primary_input_, r1cs_params.auxiliary_input_));

    /* Actual SNARK test */
    for (std::size_t i = 1; i < 2; i++) {
        const bool make_zk = (i == 0) ? false : true;
        aurora_snark_parameters<FieldT, hash_type> params(
            security_parameter,
            ldt_reducer_soundness_type,
            fri_soundness_type,
            blake2b_type,
            FRI_localization_parameter,
            RS_extra_dimensions,
            make_zk,
            domain_type,
            num_constraints,
            num_variables);
        const aurora_snark_argument<FieldT, hash_type> argument = aurora_snark_prover<FieldT>(
            r1cs_params.constraint_system_,
            r1cs_params.primary_input_,
            r1cs_params.auxiliary_input_,
            params);

        printf("iop size in bytes %lu\n", argument.IOP_size_in_bytes());
        printf("bcs size in bytes %lu\n", argument.BCS_size_in_bytes());
        printf("argument size in bytes %lu\n", argument.size_in_bytes());
        const bool bit = aurora_snark_verifier<FieldT>(
            r1cs_params.constraint_system_,
            r1cs_params.primary_input_,
            argument,
            params);

        EXPECT_TRUE(bit) << "failed on make_zk = " << i << " test";
    }
}

}

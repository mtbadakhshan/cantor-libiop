#include <algorithm>
#include <gtest/gtest.h>
#include <libff/algebra/fields/binary/gf256.hpp>
#include "libiop/protocols/ldt/fri/fri_aux.hpp"
#include "libiop/protocols/ldt/fri/fri_local_fold.hpp"
#include "libiop/snark/aurora_snark.hpp"
#include "libiop/relations/examples/r1cs_examples.hpp"

namespace libiop {

TEST(FRILocalFold, MatchesLagrangeOnArbitraryWords)
{
    using FieldT = libff::gf256;
    // No degree assumption: compare the full local linear map, including a
    // challenge inside the domain, on standard, Cantor, and non-native bases.
    for (bool cantor : {false, true})
    {
        const field_subset<FieldT> native(256, cantor);
        for (bool change_basis : {false, true})
        {
            auto basis = native.basis();
            if (change_basis)
            {
                std::reverse(basis.begin(), basis.end());
                for (auto &b : basis) b *= FieldT(7);
            }
            const field_subset<FieldT> domain(affine_subspace<FieldT>(basis, FieldT::random_element()));
            auto word = std::make_shared<std::vector<FieldT>>(domain.num_elements());
            for (auto &v : *word) v = FieldT::random_element();
            for (size_t eta : {1, 2, 3, 4, 8})
            {
                const auto tables = make_fri_local_fold_tables(domain, eta);
                for (const auto &challenge : {FieldT::random_element(), domain.element_by_index(7)})
                {
                    SCOPED_TRACE(::testing::Message() << "cantor=" << cantor
                        << " changed=" << change_basis << " eta=" << eta);
                    const auto expected = evaluate_next_f_i_over_entire_domain(word, domain, size_t(1) << eta, challenge);
                    EXPECT_EQ(evaluate_fri_local_fold(*word, tables, challenge), *expected);
                }
            }
        }
    }
}

TEST(FRILocalFold, AuroraProofsVerifyWithOriginalArithmeticSettings)
{
    using FieldT = libff::gf256;
    const auto previous_basis = get_polynomial_basis_config();
    const auto previous_fold = get_fri_fold_implementation();
    const auto example = generate_r1cs_example<FieldT>(512, 31, 511);
    for (bool cantor : {false, true})
    {
        for (size_t eta : {1, 2, 3, 4})
        {
            SCOPED_TRACE(::testing::Message() << "cantor=" << cantor << " eta=" << eta);
            set_polynomial_basis_config(polynomial_basis_config::lch_poly_basis);
            set_fri_fold_implementation(fri_fold_implementation::local_evaluations);
            const aurora_snark_parameters<FieldT, binary_hash_digest> params(
                128, LDT_reducer_soundness_type::proven, FRI_soundness_type::proven,
                blake2b_type, eta, 3, true, affine_subspace_type, cantor, 512, 511);
            const auto proof = aurora_snark_prover<FieldT>(example.constraint_system_,
                example.primary_input_, example.auxiliary_input_, params);
            // The new selector and LCH representation must not be needed to
            // interpret the transcript, including the final polynomial.
            set_polynomial_basis_config(polynomial_basis_config::monomial_poly_basis);
            set_fri_fold_implementation(fri_fold_implementation::configured);
            EXPECT_TRUE(aurora_snark_verifier<FieldT>(example.constraint_system_,
                example.primary_input_, proof, params));
        }
    }
    set_polynomial_basis_config(previous_basis);
    set_fri_fold_implementation(previous_fold);
}

TEST(FRILocalFold, RejectsInvalidDimensionsAndInputLength)
{
    using FieldT = libff::gf256;
    const field_subset<FieldT> domain(256, false);
    EXPECT_THROW(make_fri_local_fold_tables(domain, 0), std::invalid_argument);
    EXPECT_THROW(make_fri_local_fold_tables(domain, 9), std::invalid_argument);
    const auto tables = make_fri_local_fold_tables(domain, 4);
    EXPECT_THROW(evaluate_fri_local_fold(std::vector<FieldT>(255), tables, FieldT::one()), std::invalid_argument);
    EXPECT_THROW(evaluate_fri_local_fold(std::vector<FieldT>(128), tables, FieldT::one()), std::invalid_argument);
}

} // namespace libiop

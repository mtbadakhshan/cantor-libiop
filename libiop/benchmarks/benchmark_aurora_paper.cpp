// One proof per sample; export libff counters directly, never parse console logs.
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <libff/common/profiling.hpp>
#include "libiop/snark/aurora_snark.hpp"
#include "libiop/relations/examples/r1cs_examples.hpp"
#include "libiop/algebra/fft.hpp"

int main(int argc, char **argv)
{
    using namespace libiop;
    using FieldT = libff::gf256;
    try {
        if (argc != 7) throw std::invalid_argument("Usage: benchmark_aurora_paper LOG_SIZE REPS WARMUPS ETA RS_EXTRA OUTPUT.csv");
#ifdef MULTICORE
        throw std::invalid_argument("Use MULTICORE=OFF: libff stage counters are shared mutable state.");
#endif
        const int size = std::stoi(argv[1]), reps = std::stoi(argv[2]), warmups = std::stoi(argv[3]);
        const int eta = std::stoi(argv[4]), rs_extra = std::stoi(argv[5]);
        if (size < 9 || size > 19 || reps < 1 || warmups < 1 || eta < 1 || eta > 4 || rs_extra < 1 || rs_extra > 8)
            throw std::invalid_argument("Invalid size/repetitions/warmups/eta/rate parameters");
        std::ofstream out(argv[6]);
        if (!out) throw std::runtime_error("Cannot open output CSV");
        out << "log_constraints,config,sample,warmup,codeword_dim,stage,seconds,calls\n" << std::setprecision(12);
        const size_t n = size_t(1) << size;
        // All four configurations use exactly the same R1CS instance and witness.
        const auto example = generate_r1cs_example<FieldT>(n, 31, n - 1);
        const char *names[] = {"general_gm", "general_lch", "cantor_bsg", "cantor_lch"};
        libff::inhibit_profiling_info = true;
        libff::inhibit_profiling_counters = false;
        for (int sample = 0; sample < warmups + reps; ++sample) {
            // Rotate order to reduce systematic drift between the four columns.
            for (int position = 0; position < 4; ++position) {
                const int config = (position + sample) % 4;
                const bool lch = config == 1 || config == 3;
                set_polynomial_basis_config(lch ? polynomial_basis_config::lch_poly_basis : polynomial_basis_config::monomial_poly_basis);
                set_additive_fft_cantor_implementation(additive_fft_cantor_implementation::lch_afft);
                const aurora_snark_parameters<FieldT, binary_hash_digest> params(
                    128, LDT_reducer_soundness_type::proven, FRI_soundness_type::proven,
                    blake2b_type, eta, rs_extra, true, affine_subspace_type, config >= 2, n, n - 1);
                libff::clear_profiling_counters();
                const auto start = std::chrono::steady_clock::now();
                const auto proof = aurora_snark_prover<FieldT>(example.constraint_system_, example.primary_input_, example.auxiliary_input_, params);
                const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
                const auto times = libff::cumulative_times;
                const auto counts = libff::invocation_counts;
                const auto count = [&](const std::string &name) { auto i = counts.find(name); return i == counts.end() ? size_t(0) : i->second; };
                if ((count("FRI LCH fold") > 0) != lch || (count("FRI evaluation fold") > 0) == lch)
                    throw std::runtime_error("Actual FRI branch does not match requested configuration");
                // Verification and CSV writing are outside the prover timing.
                if (!aurora_snark_verifier<FieldT>(example.constraint_system_, example.primary_input_, proof, params))
                    throw std::runtime_error(std::string("Proof verification failed: ") + names[config]);
                const auto row = [&](const std::string &stage, double seconds, size_t calls) {
                    out << size << ',' << names[config] << ',' << sample << ',' << (sample < warmups) << ','
                        << params.iop_params_.codeword_domain_dim() << ",\"";
                    for (char c : stage) { if (c == '"') out << '"'; out << c; }
                    out << "\"," << seconds << ',' << calls << '\n';
                };
                row("prover_total", elapsed, 1);
                for (const auto &entry : times) row(entry.first, entry.second * 1e-9, counts.at(entry.first));
                out.flush();
                if (!out) throw std::runtime_error("Failed writing CSV");
                std::cerr << "2^" << size << ' ' << names[config] << ' ' << (sample < warmups ? "warmup" : "measured") << ' ' << sample << ": " << elapsed << " s; verified\n";
            }
        }
    } catch (const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
}

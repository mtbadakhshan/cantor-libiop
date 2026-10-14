// Fractal prover benchmark: one proof per sample; export libff counters directly.
// Fractal in libiop supports the standard (non-Cantor) additive domain only, so the
// configurations are the monomial baseline and LCH-basis arithmetic on that basis.
// Indexing (preprocessing) is outside the timed interval; every proof is verified.
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <libff/common/profiling.hpp>
#include <libff/algebra/fields/binary/gf256.hpp>
#include "libiop/snark/fractal_snark.hpp"
#include "libiop/relations/examples/r1cs_examples.hpp"
#include "libiop/algebra/fft.hpp"
#include "libiop/protocols/ldt/fri/fri_local_fold.hpp"

int main(int argc, char **argv)
{
    using namespace libiop;
    using FieldT = libff::gf256;
    using hash_type = binary_hash_digest;
    try {
        if (argc < 7 || argc > 9) throw std::invalid_argument("Usage: benchmark_fractal_paper LOG_SIZE REPS WARMUPS ETA RS_EXTRA OUTPUT.csv [configured|local] [CONFIGS_COMMA_SEPARATED]");
#ifdef MULTICORE
        throw std::invalid_argument("Use MULTICORE=OFF: libff stage counters are shared mutable state.");
#endif
        const int size = std::stoi(argv[1]), reps = std::stoi(argv[2]), warmups = std::stoi(argv[3]);
        const int eta = std::stoi(argv[4]), rs_extra = std::stoi(argv[5]);
        if (size < 9 || size > 19 || reps < 1 || warmups < 1 || eta < 1 || eta > 4 || rs_extra < 1 || rs_extra > 8)
            throw std::invalid_argument("Invalid size/repetitions/warmups/eta/rate parameters");
        const std::string fri_mode = argc >= 8 ? argv[7] : "configured";
        if (fri_mode != "configured" && fri_mode != "local")
            throw std::invalid_argument("FRI mode must be configured or local");
        const bool local_fri = fri_mode == "local";
        set_fri_fold_implementation(local_fri ? fri_fold_implementation::local_evaluations
                                             : fri_fold_implementation::configured);
        const char *names[] = {"general_gm", "general_lch"};
        std::vector<size_t> configs;
        if (argc == 9) {
            std::istringstream input(argv[8]);
            std::string name;
            bool selected[2] = {false, false};
            while (std::getline(input, name, ',')) {
                size_t index = 0;
                while (index < 2 && name != names[index]) ++index;
                if (index == 2 || selected[index])
                    throw std::invalid_argument("Unknown or duplicate configuration: " + name);
                selected[index] = true;
                configs.push_back(index);
            }
            if (configs.empty()) throw std::invalid_argument("Empty configuration list");
        } else {
            configs = {0, 1};
        }
        std::ofstream out(argv[6]);
        if (!out) throw std::runtime_error("Cannot open output CSV");
        out << "log_constraints,config,sample,warmup,stage,seconds,calls,fri_fold,proof_bytes,verifier_seconds\n" << std::setprecision(12);
        const size_t n = size_t(1) << size;
        // All selected configurations use exactly the same R1CS instance and witness.
        const auto example = generate_r1cs_example<FieldT>(n, 31, n - 1);
        const auto cs = std::make_shared<r1cs_constraint_system<FieldT>>(example.constraint_system_);
        libff::inhibit_profiling_info = true;
        libff::inhibit_profiling_counters = false;
        for (int sample = 0; sample < warmups + reps; ++sample) {
            for (size_t position = 0; position < configs.size(); ++position) {
                const size_t config = configs[(position + sample) % configs.size()];
                const bool lch = config == 1;
                set_polynomial_basis_config(lch ? polynomial_basis_config::lch_poly_basis : polynomial_basis_config::monomial_poly_basis);
                const fractal_snark_parameters<FieldT, hash_type> params(
                    128, LDT_reducer_soundness_type::proven, FRI_soundness_type::proven,
                    blake2b_type, eta, rs_extra, true, affine_subspace_type, cs);
                // Indexing is preprocessing: outside the timed interval.
                auto index = fractal_snark_indexer(params);
                bcs_prover_index<FieldT, hash_type> prover_index = index.first;
                libff::clear_profiling_counters();
                const auto start = std::chrono::steady_clock::now();
                const auto proof = fractal_snark_prover<FieldT, hash_type>(prover_index, example.primary_input_, example.auxiliary_input_, params);
                const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
                const auto times = libff::cumulative_times;
                const auto counts = libff::invocation_counts;
                const auto count = [&](const std::string &stage) { auto i = counts.find(stage); return i == counts.end() ? size_t(0) : i->second; };
                if ((count("FRI local evaluation fold") > 0) != local_fri ||
                    (count("FRI LCH fold") > 0) != (lch && !local_fri) ||
                    (count("FRI evaluation fold") > 0) != (!lch && !local_fri))
                    throw std::runtime_error("Actual FRI branch does not match requested configuration");
                const auto vstart = std::chrono::steady_clock::now();
                const bool ok = fractal_snark_verifier<FieldT, hash_type>(index.second, example.primary_input_, proof, params);
                const double velapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - vstart).count();
                if (!ok) throw std::runtime_error(std::string("Proof verification failed: ") + names[config]);
                const size_t proof_bytes = proof.size_in_bytes();
                const auto row = [&](const std::string &stage, double seconds, size_t calls) {
                    out << size << ',' << names[config] << ',' << sample << ',' << (sample < warmups) << ",\"";
                    for (char c : stage) { if (c == '"') out << '"'; out << c; }
                    out << "\"," << seconds << ',' << calls << ',' << fri_mode << ',' << proof_bytes << ',' << velapsed << '\n';
                };
                row("prover_total", elapsed, 1);
                for (const auto &entry : times) row(entry.first, entry.second * 1e-9, counts.at(entry.first));
                out.flush();
                if (!out) throw std::runtime_error("Failed writing CSV");
                std::cerr << "2^" << size << ' ' << names[config] << " fri=" << fri_mode << ' ' << (sample < warmups ? "warmup" : "measured") << ' ' << sample
                          << ": " << elapsed << " s; verified (" << velapsed << " s); proof " << proof_bytes << " B\n";
            }
        }
    } catch (const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
}

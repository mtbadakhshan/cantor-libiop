# Four-configuration Aurora benchmark

Run commands from the repository root. This is a dedicated driver, independent
of Google Benchmark's automatic iteration calibration. Each sample is one proof.

| CSV configuration | Evaluation basis | Polynomial representation and transforms | FRI |
|---|---|---|---|
| `general_gm` | Standard, non-Cantor | Monomial; Gao–Mateer FFT | Evaluation-domain fold |
| `general_lch` | Standard, non-Cantor | LCH; butterflies and existing LCH polynomial optimizations | LCH coefficient fold |
| `cantor_bsg` | Cantor | Monomial; full LCH additive FFT including basis conversion | Evaluation-domain fold |
| `cantor_lch` | Cantor | LCH; butterflies and existing LCH polynomial optimizations | LCH coefficient fold |

“General” benchmarks the implementation's standard non-Cantor basis, not a
randomly sampled basis. `cantor_bsg` denotes the requested BSG-style algorithm
configuration in this checkout, not a reproduction using a historical BSG
revision. Monomial mode disables the LCH polynomial optimizations while still
allowing the full LCH FFT. LCH polynomial mode activates the existing optimizations
and automatically selects the LCH FRI branch. Every proof checks actual FRI
branch invocation counters and fails if they disagree with the configuration.

## Build and run

Use the same source checkout and submodules on the measurement machine, including
the new benchmark files. Build there because the project's default optimization
flags include `-march=native`.

```sh
cmake -S . -B build-paper -DCMAKE_BUILD_TYPE=Release -DMULTICORE=OFF -DUSE_ASM=ON -DDEBUG=OFF -DCPPDEBUG=OFF -DPROFILE_OP_COUNTS=OFF -DUSE_ASAN=OFF -DUSE_UBSAN=OFF
cmake --build build-paper --target benchmark_aurora_paper -j2

# Smoke test: four configurations, one warmup and one measured proof each.
python3 libiop/benchmarks/run_aurora_paper.py --build-dir build-paper --min-log 9 --max-log 9 --reps 1 --warmups 1 --out benchmark_results/paper-smoke

# Full experiment: 11 sizes, four configurations, 20 measured proofs + 1 warmup.
python3 libiop/benchmarks/run_aurora_paper.py --build-dir build-paper --min-log 9 --max-log 19 --reps 20 --warmups 1 --out benchmark_results/paper-full
```

Output directories must be new, to prevent accidental overwriting. The full
command runs 924 proofs, of which 880 are measured. `--eta 2` or `--eta 3` selects
larger localization parameters; keep each experiment in a separate directory.
`--rs-extra` defaults to 5. Default eta is 1. Other fixed settings match the
existing binary-field Aurora benchmark: GF(2^256), 128-bit security, proven
soundness settings, zero knowledge, BLAKE2b, 31 public inputs, and
2^k constraints with 2^k - 1 variables.

Sizes 9–19 mean **log2 R1CS constraint count**, not log2 codeword length. The
actual `codeword_dim` is in the raw CSV. For example, the default size-9 smoke
run has codeword dimension 16. Larger sizes need substantial memory; avoid
swapping when collecting performance measurements.

## Measurement and output

One random R1CS instance and witness per size is shared by all four configurations.
The configuration order rotates between repetitions. Warmups are additional
proofs and are retained in raw data with `warmup=1`, but excluded from summaries.
Each proof is verified; verification, instance generation, parameter construction,
and CSV writing are outside the prover interval. The timed interval includes
protocol initialization, the complete prover call and its internal cleanup,
commitments, proof of work, and transcript extraction. Proof destruction is outside
that interval. Library stage instrumentation remains enabled, with profiling
console messages suppressed; these are instrumented wall-clock measurements.
Remaining library diagnostics go to a log file, not the terminal.

* `samples_K.csv`: every stage's seconds and invocation count per proof, including
  warmups, plus the external `prover_total` wall time.
* `summary.csv`: sample count, mean, median, sample standard deviation, coefficient
  of variation, minimum and maximum, excluding warmups.
* `comparison.csv`: one row per size and stage, with the four configuration columns
  containing median seconds. Empty cells mean a stage was not executed.
* `metadata.json`, `cmake_cache.txt`, `cpuinfo.txt`, `source.diff`, and copies of
  the driver sources record the environment and experiment configuration.
* `prover_K.log`: library diagnostics for troubleshooting; **not used to extract
  timings**. Verification failure or an unexpected folding branch aborts the run.

Stage times come directly from libff's in-memory counters. They are **inclusive**
and summed across invocations within one proof before taking sample statistics.
Do not sum a parent stage with its children. In particular:

* `Submit witness oracles` includes witness randomization and A/B/Cz work.
* `Aurora encoded protocol proof` covers the encoded protocol's proof calculation.
* `LDT Reducer: Calculate and submit proof` includes FRI and its commitments.
* `evaluating next FRI codeword` compares the folding loop across all four
  configurations. LCH's initial interpolation is outside that loop and is reported
  separately; include it when comparing the full arithmetic cost of the paths.
* `FRI LCH fold` contains coefficient mixing, output evaluation, and (on general
  domains) next-round basis conversion. Cantor does not execute the conversion.
* `Construct Merkle tree` and FFT-wrapper rows are cross-cutting subtotals already
  included in protocol stages. FFT-wrapper totals do not include direct butterfly
  calls in the LCH FRI implementation; use the FRI-specific rows for those.

For a compact paper table, select total prover time, witness construction,
encoded protocol proof, LDT/FRI proof, and Merkle construction, clearly labeling
the latter as an included substage. A separate arithmetic table can show initial
FRI interpolation, fold-loop cost, and its LCH substages. These are timing
comparisons, not field-operation counts or a disjoint accounting of total time.

Twenty measured samples are a starting point, not an accuracy guarantee. Inspect
the variation, run on an otherwise idle machine, and increase repetitions when
needed. One hundred runs are generally more statistically informative than twenty
under comparable conditions; neither count corrects systematic thermal, frequency,
or scheduling effects. The smoke samples are correctness checks, not paper results.

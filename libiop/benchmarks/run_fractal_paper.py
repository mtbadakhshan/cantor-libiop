#!/usr/bin/env python3
"""Run selected Fractal configurations and export sample statistics and stage tables.

Fractal in libiop supports only standard-basis (general) domains, so the configurations are
general_gm (libiop) and general_lch (LCH-basis arithmetic). The indexer runs outside the timed
interval. --fri-fold local uses local evaluation folding in each selected column; configured
keeps libiop's Lagrange-based fold for general_gm. The default, both, measures
general_gm:configured, general_gm:local, and general_lch:local in one process, interleaved
proof by proof, and writes speedups.csv with the median prover times and the speedups of
general_lch over general_gm with each fold. Every proof is verified; the verifier time is
exported as the stage verifier_total, and the peak memory of each size is written to memory.csv.
"""
import argparse
import csv
import datetime
import hashlib
import shutil
import json
import os
from pathlib import Path
import platform
import resource
import statistics
import subprocess

CONFIGS = ['general_gm', 'general_lch']

def measure(args, root, binary, out, fri_fold, configs):
    out.mkdir(parents=True, exist_ok=False)  # Never overwrite earlier measurements.
    def git(*options):
        return subprocess.check_output(['git', *options], cwd=root, text=True)
    arguments = {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()}
    arguments.update(fri_fold=fri_fold, configs=configs)
    metadata = dict(arguments=arguments,
                    timestamp_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                    platform=platform.platform(), processor=platform.processor(),
                    cpu_affinity=sorted(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else None,
                    revision=git('rev-parse', 'HEAD').strip(), git_status=git('status', '--short'),
                    submodules=git('submodule', 'status'),
                    binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                    field='GF(2^256)', security_bits=128, zero_knowledge=True,
                    fri_fold=fri_fold,
                    soundness='proven', public_inputs=31,
                    stage_semantics='inclusive wall times; nested rows overlap; warmups excluded from summaries')
    for name, path in [('cmake_cache', args.build_dir / 'CMakeCache.txt'), ('cpuinfo', Path('/proc/cpuinfo'))]:
        if path.exists():
            (out / (name + '.txt')).write_text(path.read_text())
    (out / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    (out / 'source.diff').write_text(git('diff', 'HEAD'))
    for name in ['benchmark_fractal_paper.cpp', 'run_fractal_paper.py']:
        shutil.copyfile(Path(__file__).parent / name, out / name)
    for name in ['fri_local_fold.hpp', 'fri_ldt.tcc']:
        shutil.copyfile(root / 'libiop' / 'protocols' / 'ldt' / 'fri' / name, out / name)
    groups = {}
    memory = []
    for size in range(args.min_log, args.max_log + 1):
        raw = out / f'samples_{size}.csv'
        cmd = [str(binary), str(size), str(args.reps), str(args.warmups), str(args.eta), str(args.rs_extra), str(raw.resolve()), fri_fold, ','.join(configs)]
        print('Running:', ' '.join(cmd), flush=True)
        before = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        with (out / f'prover_{size}.log').open('w') as log:
            subprocess.run(cmd, stdout=log, check=True)
        memory.append((size, resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss, before))
        with (out / 'memory.csv').open('w') as f:
            w = csv.writer(f)
            w.writerow(['log_constraints', 'peak_rss_kb_so_far', 'peak_rss_kb_before_this_size'])
            w.writerows(memory)
        with raw.open() as f:
            rows = list(csv.DictReader(f))
        if not rows or any(row.get('fri_fold') != (row['config'].split(':')[1] if ':' in row['config'] else fri_fold) for row in rows):
            raise RuntimeError('FRI mode missing or inconsistent; rebuild benchmark_fractal_paper')
        if {row['config'] for row in rows} != set(configs):
            raise RuntimeError('Measured configurations do not match the requested selection')
        if len({row['proof_bytes'] for row in rows}) != 1:
            raise RuntimeError('Proof sizes differ between configurations')
        for config in configs:
            totals = [r for r in rows if r['config'] == config and r['stage'] == 'prover_total']
            if len(totals) != args.reps + args.warmups or sum(r['warmup'] == '0' for r in totals) != args.reps:
                raise RuntimeError(f'Unexpected sample count for {config}, size {size}')
        for row in rows:
            if row['warmup'] == '0':
                key = (int(row['log_constraints']), row['stage'], row['config'])
                groups.setdefault(key, []).append(float(row['seconds']))
                if row['stage'] == 'prover_total':
                    groups.setdefault((key[0], 'verifier_total', key[2]), []).append(float(row['verifier_seconds']))
        # Update summaries after each completed size, retaining raw data on failure.
        with (out / 'summary.csv').open('w') as f:
            w = csv.writer(f)
            w.writerow(['log_constraints', 'stage', 'config', 'samples', 'mean_s', 'median_s', 'stdev_s', 'cv_percent', 'min_s', 'max_s', 'fri_fold'])
            for (sz, stage, config), values in sorted(groups.items()):
                mean = statistics.mean(values)
                sd = statistics.stdev(values) if len(values) > 1 else 0
                w.writerow([sz, stage, config, len(values), mean, statistics.median(values), sd, 100 * sd / mean if mean else 0, min(values), max(values), fri_fold])
        with (out / 'comparison.csv').open('w') as f:
            w = csv.writer(f)
            w.writerow(['log_constraints', 'stage', *[c + '_median_s' for c in configs], 'fri_fold'])
            for sz, stage in sorted({(k[0], k[1]) for k in groups}):
                w.writerow([sz, stage, *[statistics.median(groups[(sz, stage, c)]) if (sz, stage, c) in groups else '' for c in configs], fri_fold])
    print(f'Results: {out}/comparison.csv ({len(configs)} configuration columns, inclusive stage medians)')
    return groups

def main():
    root = Path(__file__).resolve().parents[2]
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--build-dir', type=Path, default=root / 'build')
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--min-log', type=int, default=9)
    p.add_argument('--max-log', type=int, default=19)
    p.add_argument('--reps', type=int, default=20, help='Measured proofs per configuration and size')
    p.add_argument('--warmups', type=int, default=1, help='Additional discarded proofs per configuration and size')
    p.add_argument('--eta', type=int, default=1, choices=[1, 2, 3, 4])
    p.add_argument('--rs-extra', type=int, default=5)
    p.add_argument('--fri-fold', choices=['both', 'configured', 'local'], default='both',
                   help='both: configured general_gm, then local with --configs; configured: original fold; local: local evaluation folding')
    p.add_argument('--configs', nargs='+', choices=CONFIGS, default=CONFIGS,
                   help='Configurations to measure (default: both)')
    args = p.parse_args()
    if len(set(args.configs)) != len(args.configs):
        p.error('Duplicate configurations are not allowed')
    if not (9 <= args.min_log <= args.max_log <= 19 and args.reps > 0 and args.warmups > 0 and 1 <= args.rs_extra <= 8):
        p.error('Invalid range, repetitions, warmups or RS extra dimensions')
    if args.fri_fold == 'both' and set(args.configs) != set(CONFIGS):
        p.error('--fri-fold both measures both configurations')
    binary = args.build_dir.resolve() / 'libiop/benchmark_fractal_paper'
    if not binary.is_file():
        p.error(f'Missing {binary}; build benchmark_fractal_paper first')
    if args.fri_fold != 'both':
        measure(args, root, binary, args.out, args.fri_fold, args.configs)
        return
    labels = ['general_gm:configured', 'general_gm:local', 'general_lch:local']
    groups = measure(args, root, binary, args.out, 'local', labels)
    with (args.out / 'speedups.csv').open('w') as f:
        w = csv.writer(f)
        w.writerow(['log_constraints', 'libiop_median_s', 'libiop_lch_fold_median_s', 'ours_median_s',
                    'speedup_over_libiop', 'speedup_over_libiop_lch_fold'])
        for size in range(args.min_log, args.max_log + 1):
            base = statistics.median(groups[(size, 'prover_total', labels[0])])
            base_fold = statistics.median(groups[(size, 'prover_total', labels[1])])
            ours = statistics.median(groups[(size, 'prover_total', labels[2])])
            w.writerow([size, base, base_fold, ours, base / ours, base_fold / ours])
    print(f'Speedups: {args.out}/speedups.csv')

if __name__ == '__main__':
    main()

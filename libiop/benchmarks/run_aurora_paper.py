#!/usr/bin/env python3
"""Run four Aurora configurations and export sample statistics and stage tables."""
import argparse
import csv
import datetime
import hashlib
import shutil
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess

CONFIGS = ['general_gm', 'general_lch', 'cantor_bsg', 'cantor_lch']

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
    args = p.parse_args()
    if not (9 <= args.min_log <= args.max_log <= 19 and args.reps > 0 and args.warmups > 0 and 1 <= args.rs_extra <= 8):
        p.error('Invalid range, repetitions, warmups or RS extra dimensions')
    binary = args.build_dir.resolve() / 'libiop/benchmark_aurora_paper'
    if not binary.is_file():
        p.error(f'Missing {binary}; build benchmark_aurora_paper first')
    args.out.mkdir(parents=True, exist_ok=False)  # Never overwrite earlier measurements.
    def git(*options):
        return subprocess.check_output(['git', *options], cwd=root, text=True)
    metadata = dict(arguments={k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
                    timestamp_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                    platform=platform.platform(), processor=platform.processor(),
                    cpu_affinity=sorted(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else None,
                    revision=git('rev-parse', 'HEAD').strip(), git_status=git('status', '--short'),
                    submodules=git('submodule', 'status'),
                    binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                    field='GF(2^256)', security_bits=128, zero_knowledge=True,
                    soundness='proven', public_inputs=31,
                    stage_semantics='inclusive wall times; nested rows overlap; warmups excluded from summaries')
    for name, path in [('cmake_cache', args.build_dir / 'CMakeCache.txt'), ('cpuinfo', Path('/proc/cpuinfo'))]:
        if path.exists():
            (args.out / (name + '.txt')).write_text(path.read_text())
    (args.out / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    (args.out / 'source.diff').write_text(git('diff', 'HEAD'))
    for name in ['benchmark_aurora_paper.cpp', 'run_aurora_paper.py']:
        shutil.copyfile(Path(__file__).parent / name, args.out / name)
    groups = {}
    for size in range(args.min_log, args.max_log + 1):
        raw = args.out / f'samples_{size}.csv'
        cmd = [str(binary), str(size), str(args.reps), str(args.warmups), str(args.eta), str(args.rs_extra), str(raw.resolve())]
        print('Running:', ' '.join(cmd), flush=True)
        with (args.out / f'prover_{size}.log').open('w') as log:
            subprocess.run(cmd, stdout=log, check=True)
        with raw.open() as f:
            rows = list(csv.DictReader(f))
        for config in CONFIGS:
            totals = [r for r in rows if r['config'] == config and r['stage'] == 'prover_total']
            if len(totals) != args.reps + args.warmups or sum(r['warmup'] == '0' for r in totals) != args.reps:
                raise RuntimeError(f'Unexpected sample count for {config}, size {size}')
        for row in rows:
            if row['warmup'] == '0':
                key = (int(row['log_constraints']), row['stage'], row['config'])
                groups.setdefault(key, []).append(float(row['seconds']))
        # Update summaries after each completed size, retaining raw data on failure.
        with (args.out / 'summary.csv').open('w') as f:
            w = csv.writer(f)
            w.writerow(['log_constraints', 'stage', 'config', 'samples', 'mean_s', 'median_s', 'stdev_s', 'cv_percent', 'min_s', 'max_s'])
            for (sz, stage, config), values in sorted(groups.items()):
                mean = statistics.mean(values)
                sd = statistics.stdev(values) if len(values) > 1 else 0
                w.writerow([sz, stage, config, len(values), mean, statistics.median(values), sd, 100 * sd / mean if mean else 0, min(values), max(values)])
        with (args.out / 'comparison.csv').open('w') as f:
            w = csv.writer(f)
            w.writerow(['log_constraints', 'stage', *[c + '_median_s' for c in CONFIGS]])
            for sz, stage in sorted({(k[0], k[1]) for k in groups}):
                w.writerow([sz, stage, *[statistics.median(groups[(sz, stage, c)]) if (sz, stage, c) in groups else '' for c in CONFIGS]])
    print(f'Results: {args.out}/comparison.csv (four columns, inclusive stage medians)')

if __name__ == '__main__':
    main()

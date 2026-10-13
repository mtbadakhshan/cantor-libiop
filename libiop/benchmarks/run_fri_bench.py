#!/usr/bin/env python3
"""Sweep the FRI-fold micro-benchmark (bench_fri_fold) over D, rate and eta.

Build (Release, single thread; on AArch64 add -DLIBFF_GF256_PMULL to CMAKE_CXX_FLAGS
for carry-less multiplication, on x86 use -DUSE_ASM=ON):

  cmake -S . -B build-bench -DCMAKE_BUILD_TYPE=Release -DMULTICORE=OFF
  cmake --build build-bench --target bench_fri_fold

Run from the repository root:

  python3 libiop/benchmarks/run_fri_bench.py --build-dir build-bench --out fri_bench_results

libiop's Lagrange fold does one inversion per coset and is slow at large codeword
domains: it runs with the full repetition count up to --lagrange-full-max-logn,
with a single measured run up to --lagrange-max-logn, and is skipped beyond.
"""
import argparse
import concurrent.futures as cf
import csv
import os
import subprocess
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build-dir", default="build-bench")
    ap.add_argument("--out", required=True, help="new output directory")
    ap.add_argument("--log-D", default="12,14,16,18,20")
    ap.add_argument("--rs-extra", default="1,2,3,4,5", help="rho = 2^-r")
    ap.add_argument("--eta", default="1,2,3,4")
    ap.add_argument("--basis", default="standard,cantor")
    ap.add_argument("--schedule", default="uniform,libiop",
                    help="uniform: [eta,...]; libiop: [1,eta,...] as in Aurora")
    ap.add_argument("--methods", default="lagrange,eval,lch")
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--warmups", type=int, default=1)
    ap.add_argument("--max-logn", type=int, default=25)
    ap.add_argument("--lagrange-full-max-logn", type=int, default=20)
    ap.add_argument("--lagrange-max-logn", type=int, default=23)
    ap.add_argument("--jobs", type=int, default=1, help="concurrent processes (1 = clean timings)")
    ap.add_argument("--validate", action="store_true", help="also check outputs against libiop's fold")
    a = ap.parse_args()

    exe = os.path.join(a.build_dir, "libiop", "bench_fri_fold")
    if not os.path.exists(exe):
        sys.exit(f"missing {exe}")
    os.makedirs(a.out, exist_ok=False)
    ints = lambda s: [int(x) for x in s.split(",") if x]

    jobs = []
    for basis in a.basis.split(","):
        for sched in a.schedule.split(","):
            for log_D in ints(a.log_D):
                for r in ints(a.rs_extra):
                    if log_D + r > a.max_logn:
                        continue
                    for eta in ints(a.eta):
                        if eta > log_D:
                            continue
                        for method in a.methods.split(","):
                            reps, warm = a.reps, a.warmups
                            if method == "lagrange":
                                if log_D + r > a.lagrange_max_logn:
                                    continue
                                if log_D + r > a.lagrange_full_max_logn:
                                    reps, warm = 1, 0
                            jobs.append((basis, sched, log_D, r, eta, method, reps, warm))

    def run(job):
        basis, sched, log_D, r, eta, method, reps, warm = job
        part = os.path.join(a.out, f"part_{basis}_{sched}_{log_D}_{r}_{eta}_{method}.csv")
        cmd = [exe, "--basis", basis, "--schedule", sched, "--log-D", str(log_D),
               "--rs-extra", str(r), "--eta", str(eta), "--methods", method,
               "--reps", str(reps), "--warmups", str(warm), "--csv", part]
        if a.validate:
            cmd.append("--validate")
        p = subprocess.run(cmd, capture_output=True, text=True)
        if p.returncode != 0:
            raise RuntimeError(" ".join(cmd) + "\n" + p.stderr)
        return part, p.stderr.strip()

    parts = []
    with cf.ThreadPoolExecutor(max_workers=a.jobs) as ex:
        for i, (part, msg) in enumerate(ex.map(run, jobs)):
            parts.append(part)
            print(f"[{i + 1}/{len(jobs)}] {msg.splitlines()[-1] if msg else part}", flush=True)

    rows, header = [], None
    for part in parts:
        with open(part) as f:
            rd = csv.reader(f)
            h = next(rd)
            header = header or h
            rows.extend(rd)
        os.remove(part)
    with open(os.path.join(a.out, "results.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)
    print(f"wrote {os.path.join(a.out, 'results.csv')} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
